####################################################################
## Upload packages
####################################################################

library(tidyverse)
library(fitdistrplus)
library(limma)
library(data.table)
library(janitor)

####################################################################
## Common functions
####################################################################

#' Left Right function to paste two strings
#'
#' @param left string
#' @param right string
#'
#' @return concatenated string
#' @export
#'
`%+%` <- function(left,right){
  return(paste0(left,right))
}

#' Left Right function to evaluate NOT in
#'
#' @param left vector
#' @param right vector
#'
#' @return boolean vector
#' @export
#'
`%nin%` <- function(left,right) {
  return(!(left %in% right))
}

#' function to evaluate NOT is.na
#'
#' @param target vector
#'
#' @return binary vector
#' @export
#'
nna <- function(target) {return(!is.na(target))}

####################################################################
## Functions for Score Calculation
####################################################################

#' Function to calculate aggregate reactivity scores
#'
#' Distributions are created by randomly selecting peptides
#' that were designed to viruses that are not the targets of
#' significant antibody responses
#' 
#' v 0.2.1 update - changed linear models to splines for interpolating
#' intermediate values of peptides
#' v 1.0 update - reduced length of 'representations' to expedite computation
#'  - fold changes are no longer log transformed
#'
#' @param norm_log dataframe containing intermediate reactivity metrics
#' @param all_peptide_fcs dataframe containing all peptide fold changes 
#' for an individual antibody profile
#' @param positives dataframe containing viruses determined to be 
#' significantly targeted by antibodies and therefore excluded from
#' random selection
#' @param exclusion_method string describing method for excluding likely reactive
#' peptides. Defaults to by genus unless exclusion_method = 'group' or 'species"
#'
#' @return dataframe containing aggregate reactivity scores and data on 
#' random distributions
#' 
#'
#' @import tidyverse fitdistrplus
calc_scores <- function(norm_log, all_peptide_fcs, positives, exclusion_method = "genus") {
  all_peptide_fcs <- all_peptide_fcs %>% mutate(fc = (value))
  
  if(exclusion_method == "species" | exclusion_method == "group"){
    all_peptide_fcs <- all_peptide_fcs %>% 
      filter(taxon_species %nin% positives$taxon_species)
  } else {
    all_peptide_fcs <- all_peptide_fcs %>% 
      filter(taxon_genus %nin% positives$taxon_genus)
  }
  print(nrow(all_peptide_fcs))
  # 1. 获取当前背景肽段的总数量（all_peptide_fcs的行数，即待抽样的总体大小）
  current_bg_size <- nrow(all_peptide_fcs)  # 关键：每轮迭代的背景肽段数量可能不同
  
  # 2. 获取实际数据中total_peps的最大值（需覆盖的目标范围）
  max_total_peps <- max(norm_log$total_peps, na.rm = TRUE)
  
  representations <- c()
  n <- 1
  t <- 2
  while(TRUE) {
    current_value <- 10 * (t^n)  # 生成10*t^n的值
    
    # 若当前值超过背景肽段数量，停止（避免抽样时超过总体）
    if(current_value > current_bg_size) {
      break
    }
    
    representations <- c(representations, current_value)
    
    # 若当前值已超过max_total_peps，且下一个值会超过背景量，可提前停止
    next_value <- 10 * (t^(n+1))
    if(current_value > max_total_peps && next_value > current_bg_size) {
      break
    }
    
    n <- n + 1
  }
  if(length(representations) == 0) {
    # 若背景肽段数量≥10，用10作为保底（最小可能的10*2^1=20不满足时，退到10）
    if(current_bg_size >= 10) {
      representations <- 10
    } else {
      # 若背景肽段数量<10，直接用当前背景量作为唯一值（允许有放回抽样）
      representations <- current_bg_size
    }
  }
  print(next_value)
  
  distributions <- data.frame(matrix(nrow = length(representations), ncol = 1000))
  row.names(distributions) <- representations
  
  #calculate random distributions
  set.seed(120)
  for(i in c(1:length(representations))){
    n = row.names(distributions)[i]
    for(j in c(1:1000)){
      distributions[i, j] = mean(sample(all_peptide_fcs$fc, as.numeric(n), replace = F))
      if(distributions[i, j] == 0) {distributions[i,j] = .001} 
    }
    distributions[i, 1000] = mean(all_peptide_fcs$fc)
    #print("calculating distributions " %+% i %+% " of " %+% length(representations))
  }
  
  #model random distributions as a gamma distribution
  dist_info <- distributions[,c(1, 2, 3)]
  names(dist_info) <- c("total_peps", "shape", "rate")
  for(i in c(1:length(representations))){
    
    dist_info[i, 1] <- row.names(dist_info)[i] %>% as.numeric()
    
    gammafit  <-  fitdist(distributions[i,] %>% t() %>% as.numeric(), "gamma")
    
    dist_info[i, 2] <- gammafit[["estimate"]][["shape"]]
    dist_info[i, 3] <- gammafit[["estimate"]][["rate"]]
  }
  
  dist_info <- dist_info %>% mutate(mean = shape / rate) %>%
    mutate(variance = shape / rate^2)

  # 数据检查与清洗
  dist_info_clean <- dist_info %>%
    filter(
      is.finite(total_peps),
      is.finite(shape),
      is.finite(rate),
      is.finite(mean)
    ) %>% distinct(total_peps, .keep_all = TRUE)
  
  # 检查有效数据点数量（至少需要4个不同的x值）
  x_values <- log10(dist_info_clean$total_peps)
  valid_x_count <- length(unique(x_values))
  if(valid_x_count < 4) {
    warning(paste0("当前迭代数据不足：有效x值仅", valid_x_count, "个（需至少4个），退出该轮迭代"))
    # 返回包含错误标记的列表，而非norm_log（方便上层函数识别）
    return(list(error = TRUE, message = "数据点不足，无法拟合样条", data = dist_info_clean))
  } else { dist_info <- dist_info_clean }
  
  # gamma distribution parameters vary with n
  #shape_lm <- lm(log10(shape) ~ I(log10(total_peps)), data = dist_info)
  shape_spline <- smooth.spline(x = log10(dist_info$total_peps), y = log10(dist_info$shape), spar = .5)
  #rate_lm <- lm(log10(rate) ~ I(log10(total_peps)), data = dist_info)
  rate_spline <- smooth.spline(x = log10(dist_info$total_peps), y = log10(dist_info$rate), spar = .5)
  #mean_lm <- lm(log10(mean) ~ I(log10(total_peps)), data = dist_info)
  mean_spline <- smooth.spline(x = log10(dist_info$total_peps), y = log10(dist_info$mean), spar = .5)
  
  # add distribution info to norm_log
  norm_log <- norm_log %>%
    mutate(shape = 10^(predict(shape_spline, log10(total_peps))$y),
           rate = 10^(predict(rate_spline, log10(total_peps))$y)) %>%
    mutate(mean = shape / rate) %>%
    mutate(variance = shape / rate^2)
  
  # calculate scores and p values
  norm_log <- norm_log %>%
    mutate(ARscore = zscoreGamma(score_norm, shape = shape / sqrt(total_peps), rate = rate / sqrt(total_peps)),
           virus_fc = score_norm / mean) %>%
    mutate(ARscore = ifelse(ARscore < -10, -10, ARscore)) %>%
    mutate(p_val = -log10((1-pnorm(abs(zscoreGamma(score_norm, shape = shape, rate = rate)))))) %>%
    mutate(p_val = ifelse(p_val>15, 15, p_val))
  
  return(list(error = FALSE, result = norm_log))  # 正常返回结果
}


#' Shell function to iterate calculations of aggregate reactivity scores
#'
#' v 1.0 update - two internal hit criteria are now used
#'
#' @param norm_log_1 dataframe containing intermediate reactivity metrics
#' @param all_peptide_fcs_1 dataframe containing all peptide fold changes 
#' for an individual antibody profile
#' @param max_iterations integer limiting the number of itererations to run
#' calc_scores. defaults to 10
#' @param p_cutoff numeric -log10 p value cutoff to determine significant  
#' targets of antibody reactivity. defaults to 4
#' @param score_cutoff numeric aggregate reactivity score cutoff to determine
#' signifiant targets of antibody reactivity. defaults to 2 the empircally 
#' determined threshold for VRC VirScan data. Recommended to adjust
#' @param exclusion_method_1 string describing method for excluding likely reactive
#' peptides. Defaults to by genus unless exclusion_method = 'group' or 'species"
#'
#' @return returns aggregate reactivity scores as well as a list denoting the
#' significant targets of reactivity detected during each iteration
#' 
#' @import tidyverse limma
iterative_scores <- function(norm_log_1, all_peptide_fcs_1, max_iterations = 10,
                             p_cutoff = -log10(.0001), score_cutoff = 2,
                             exclusion_method_1 = "genus") {
  iterations = 0
  positives_1 = norm_log_1 %>% filter(F)
  positives_output = list()
  library(limma)
  
  while(iterations < max_iterations) {
    print("iteration " %+% (iterations+1))
    # 调用calc_scores并检查是否返回错误
    calc_result <- calc_scores(
      norm_log = norm_log_1, 
      all_peptide_fcs = all_peptide_fcs_1, 
      positives = positives_1,
      exclusion_method = exclusion_method_1
    )
    
    # 若calc_scores返回错误，退出当前迭代
    if(calc_result$error) {
      warning("calc_scores出错：", calc_result$message, "，停止当前样本的迭代")
      # 可返回截至上一轮的结果（若有）或空结果
      final_scores <- if(iterations > 0) positives_output[[iterations]] else NULL
      return(list(final_scores, positives_output))
    }
    
    # 数据正常时，继续迭代逻辑
    scores <- calc_result$result  # 提取正常结果
    
    #update variables
    positives_2 <- scores %>% filter(p_val == 15) %>% filter(ARscore > 0) #12/02/23 added second internal criteria
    positives_1 <- scores %>% filter(p_val > p_cutoff) %>% filter(ARscore > score_cutoff) %>% full_join(positives_2, 
                            join_by(taxon_species, taxon_genus, sample_id, total_peps, score, score_norm, shape, 
                            rate, mean, variance, ARscore, virus_fc, p_val))
    positives_output[[iterations+1]] <- positives_1 
    
    #add break if no new positives
    if(iterations > 0) {
      if((positives_1$taxon_species %>% length()) <= (positives_output[[iterations]]$taxon_species %>% length())) {
        iterations <- iterations + 10
      }
    }
    
    iterations <- iterations + 1
  }
  
  outputs <- list(scores, positives_output)
  return(outputs) 
}

#' Calculate Aggregate Reactivity Scores
#'
#'v 1.0 update - reduced length of 'representations' to expedite computation
#'  - fold changes are no longer log transformed
#'  - additional internal positive criteria added
#'  - default thresholds changed to reflect shift to linear fcs
#'  - option to leave hfc as null
#'  - added option to choose peptide exclusion method
#'  
#'  If custom peptide groupings are desired, replace the taxon_species column
#'  with desired peptide grouping variable
#'
#' @param hfc a dataframe containing hits fold change data from a PhIP-Seq
#' experiment. Used to exclude peptides that are reactive in beads only controls.
#' If left NULL, no peptides will be excluded
#' @param fc a dataframe containing fold change data from a PhIP-Seq experiment.
#' Used to generate intermediate reactivity metrics, null distributions, and 
#' ARscores
#' @param set_max_iterations integer limiting the number of itererations to run
#' calc_scores. defaults to 10
#' @param set_p_cutoff numeric -log10 p value cutoff to determine significant  
#' targets of antibody reactivity. defaults to 4
#' @param set_score_cutoff numeric aggregate reactivity score cutoff to 
#' determine signifiant targets of antibody reactivity. defaults to 2 the 
#' empircally determined threshold for VRC VirScan data. Recommended to adjust
#' @param required_number_of_peptides integer that represents the minimum number
#' of peptides attributed to an antigen/pathogen required to calculate an 
#' ARscore
#' @param exclusion_method string describing method for excluding likely reactive
#' peptides. Defaults to by genus unless exclusion_method = 'group' or 'species"
#'
#' @return a dataframe containing the aggregate reactivity scores for a PhIP run
#' @export
#'
#' @import tidyverse fitdistrplus limma
ARscore_algorithm <- function(hfc = NULL, fc, set_max_iterations = 10,
                               set_p_cutoff = -log10(0.01), set_score_cutoff = 2,
                               required_number_of_peptides = 40,
                              exclusion_method = "genus") {
  
  ## peptides that have hits in beads
  if(!is.null(hfc)){
    bad_beads <- hfc %>% dplyr::select(pep_aa,contains('BEADS'),contains('Beads'),contains('beads')) %>% pivot_longer(-pep_aa) %>%
      filter(value > 1) %>% dplyr::select(pep_aa) %>% pull
    
    ### fold change ###
    d1 <- fc
    
    ## make longform dataframe out of fc, no longer floored fc at 1
    d2 <- d1 %>% filter(pep_aa %nin% bad_beads) %>% 
      dplyr::select(-contains('PBS')) %>%
      #common pulldowns. Has to be edited in some screens
      pivot_longer(cols = c(contains('-')),names_to = 'sample_id')
    
  } else {
    
    ### fold change ###
    d1 <- fc
    
    ## make longform dataframe out of fc, no longer floored fc at 1
    d2 <- d1 %>% 
      dplyr::select(-contains('PBS')) %>%
      #common pulldowns. Has to be edited in some screens
      pivot_longer(cols = c(contains('-')),names_to = 'sample_id')
    
  }
  
  ids <- d2 %>% dplyr::select(sample_id) %>% unique()
  
  # need to group by taxon_species, sample_id and count number of reactivities. Then add back sample annotations 
  d4 <- d2 %>% group_by(taxon_species, taxon_genus, sample_id) %>% summarise(score = sum((value))) %>% left_join(ids)
  d4_2 <- d2 %>% group_by(taxon_species, taxon_genus, sample_id) %>% summarise(total_peps = n()) %>%
    left_join(ids) %>% left_join(d4) %>% mutate(score = ifelse(is.na(score), 0, score))
  d4_2 <- d4_2 %>% mutate(score_norm = score / total_peps) %>% filter(total_peps >= required_number_of_peptides)
  
  library(tidyverse)
  library(fitdistrplus)
  library(limma)
  
  algorithm_output <- list()
  
  for(R in c(1:length(ids$sample_id))) {
    print(ids$sample_id[R])
    algorithm_output[[R]] <- iterative_scores(norm_log_1 = d4_2 %>% filter(sample_id == ids$sample_id[R]), 
                                              all_peptide_fcs_1 = d2 %>% filter(sample_id == ids$sample_id[R]),
                                              max_iterations = set_max_iterations,
                                              p_cutoff = set_p_cutoff, 
                                              score_cutoff = set_score_cutoff,
                                              exclusion_method_1 = exclusion_method)
    
  }
  
  #algorithm output has the final output as well as a structure that tracks positive 
  scores <- algorithm_output[[1]][[1]]
  if(length(algorithm_output)>1){
  for(i in c(2:length(algorithm_output))){
    scores <- rbind(scores, algorithm_output[[i]][[1]])
  }
  }
  
  scores <- scores %>% dplyr::select(taxon_species, taxon_genus, sample_id,
                                     total_peps, score_norm, 
                                     ARscore, p_val)
  
  return(scores)
}












