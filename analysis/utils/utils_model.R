#' Functions needed for modeling

#' Extract LOCF outcomes for a specific visit week
#'
#' @param data A data frame with preprocessed data from all trials
#' @param outcome_base_code String that identifies the outcome of interest, e.g. 
#'   "MDS_UPDRS23_total_diff" no week or locf should be included here
#' @return A string ID that identifies NCT02337738
extract_visit_locf <- function(data, outcome_base_code, week_number) {
  
  data_week <- data |>
    select(trial, all_of(arm_fac), all_of(covariate),
           all_of(paste0(outcome_base_code, "_w", week_number, "_locf"))) |>
    mutate(week = week_number) |>
    rename(outcome_locf = all_of(
      paste0(outcome_base_code, "_w", week_number, "_locf"))
    )
  return(data_week)
}

#' Extract needed LOCF outcomes and associated visit weeks from the data
#'
#' @param data A data frame with preprocessed data from all trials
#' @param outcome_base_code String that identifies the outcome of interest, e.g. 
#'   "MDS_UPDRS23_total_diff" no week or locf should be included here
#' @param weeks_of_interest A list of numbers indicating weeks of interest
#' @return A data frame LOCF outcomes, visit weeks and all other information 
#'   needed to calculate average treatment effects
extract_all_visits_locf <- function(data, outcome_base_code, weeks_of_interest) {
  
  data_locf_list <- vector("list", length(weeks_of_interest))
  
  # difference from baseline at baseline is zero
  data_locf_list[[1]] <- data |>
    select(trial, all_of(arm_fac), all_of(covariate)) |>
    mutate(week = 0) |>
    mutate(outcome_locf = 0)
  
  for (i_week in 2:length(weeks_of_interest)) {
    
    # outcome difference from baseline at other weeks
    data_locf_list[[i_week]] <- extract_visit_locf(
      data = data, 
      outcome_base_code = outcome_base_code, 
      week = weeks_of_interest[i_week]
    )
  }
  data_locf <- do.call(rbind, data_locf_list)
  return(data_locf)
}


#' Calculate average observed treatment effects per visit and per trial by using 
#' a last observation carried forward principle
#'
#' @param data A data frame with preprocessed data from all trials
#' @param outcome_base_code String that identifies the outcome of interest, e.g. 
#'   "MDS_UPDRS23_total_diff" no week or locf should be included here
#' @param covariate_code A string indicating the covariate used to correct 
#'   average treatment effect
#' @param arm_code A string indicating the treatment arm
#' @param weeks_of_interest A list of numbers indicating weeks of interest
#' @return A data frame with average treatment effects
calculate_average_effects_locf <- function(data, 
                                           outcome_base_code, 
                                           covariate_code, 
                                           arm_code, 
                                           weeks_of_interest) {
  
  data_locf <- extract_all_visits_locf(data, outcome_base_code, weeks_of_interest)
  
  weeks <- unique(data_locf$week)
  trials <- unique(data_locf$trial)
  
  average_effect_list <- vector("list", length(weeks))
  
  for (i_week in 1:length(weeks)){
    
    average_effect_list[[i_week]] <- vector("list", length(trials))
    
    for (i_trial in 1:length(trials)){
      
      current_data <- data_locf |>
        filter(week == weeks[i_week]) |>
        filter(trial == trials[i_trial])
      
      if(!(sum(is.na(current_data)) > 200)){
        
        current_average_effect <- calculate_average_effect(
          data = current_data, 
          outcome_code = "outcome_locf",
          covariate_code = covariate_code,
          arm_code = arm_code)
        
        to_add = data.frame(week = weeks[i_week], trial = trials[i_trial])
        
        current_average_effect <- cbind(current_average_effect, to_add)
        average_effect_list[[i_week]][[i_trial]] <- current_average_effect
      }
    }
  }
  average_effects <- do.call(rbind, unlist(average_effect_list, recursive = FALSE))
  return(average_effects)
}

#' Determine if the average treatment effect, expressed as a difference in 
#' change from baseline in MDS-UPDRS score is sufficiently similar
#'
#' @param ids A list of IDs that determine a training-validation split
#' @param data A data frame with the current trial data to split
#' @param outcome_code A string indicating the outcome
#' @param covariate_code A string indicating the covariate used to correct 
#'   average treatment effect
#' @param arm_code A string indicating the treatment arm
#' @param threshold A number indicating the threshold for similarity, a smaller 
#'   number will require a higher similarity in average treatment effect
#' @return TRUE if sufficiently similar, FALSE otherwise
is_similar_difference <- function(ids, 
                                  data, 
                                  outcome_code, 
                                  covariate_code, 
                                  arm_code, 
                                  threshold = 0.02) {
  
  formula <- as.formula(paste(outcome_code, " ~ ", covariate_code, " + ", arm_code))
  
  is_similar <- TRUE
  
  # calculate average treatment effect in the training set
  ancova_fit_t <- lm(formula, data = data[ids[[1]],])
  results_t <- confint(pairs(emmeans(ancova_fit_t, arm_code), reverse = TRUE))
  difference_t <- subset(summary(results_t), select = "estimate")[[1]]
  
  # calculate average treatment effect in the validation set
  ancova_fit_v <- lm(formula, data = data[ids[[2]],])
  results_v <- confint(pairs(emmeans(ancova_fit_v, arm_code), reverse = TRUE))
  difference_v <- subset(summary(results_v), select = "estimate")[[1]]
  
  # check if the average effects are similar
  if(abs((difference_t - difference_v) / difference_v) > threshold){
    is_similar <- FALSE
  }
  return(is_similar)
}

#' Search for a seed for the current trial that result in a split that has a 
#' similar average treatment effect in training and validation set
#'
#' @param data A data frame with the current trial data to split
#' @param outcome_code A string indicating the outcome
#' @param covariate_code A string indicating the covariate used to correct 
#'   average treatment effect
#' @param arm_code A string indicating the treatment arm
#' @param threshold A number indicating the threshold for similarity, a smaller 
#'   number will require a higher similarity in average treatment effect
#' @param nb_splits A number to indicate how many times the data should be split
#' @return A data frame with split seeds for each trial
search_cv_seeds_one_trial <- function(data_trial, 
                                      trial_code,
                                      outcome_code, 
                                      covariate_code, 
                                      arm_code, 
                                      threshold = 0.05,
                                      nb_splits = 100) {
  
  seeds_loop <- list()
  
  while(length(seeds_loop) < nb_splits) {
    while(TRUE){
      
      rs <- sample(1:1000000000, 1)
      
      # split the data in a training and validation set
      split_ids <- splitTools::partition(
        y = pull(data_trial, arm_code), 
        p = c(0.6, 0.4), 
        seed = rs
      )
      
      valid_split <- TRUE
      
      # check if the average effects in training and validation set are similar
      if(!is_similar_difference(
        ids = split_ids, 
        data = data_trial, 
        outcome_code = outcome_code, 
        covariate_code = covariate_code,
        arm_code =  arm_code, 
        threshold = threshold
      )){
        valid_split <- FALSE
      }
      
      #print(rs)
      
      # search until a split with similar average effects is found
      if(valid_split){
        #print(rs)
        break
      }
    }
    
    # add the seed to the list
    if (!rs %in% seeds_loop){
      seeds_loop <- append(seeds_loop, rs)
      cat(trial_code, "appended seed", rs, "current length", length(seeds_loop),"\n")
    }
  }
  
  new_seeds <- unlist(seeds_loop)
  
  seeds_trial <- data.frame(setNames(list(new_seeds), trial_code))
  return(seeds_trial)
}


#' Search for seeds that result in a split that has a similar average treatment 
#' effect in training and validation set
#'
#' @param data A data frame with all trial data to split
#' @param outcome_code A string indicating the outcome
#' @param covariate_code A string indicating the covariate used to correct 
#'   average treatment effect
#' @param arm_code A string indicating the treatment arm
#' @param threshold A number indicating the threshold for similarity, a smaller 
#'   number will require a higher similarity in average treatment effect
#' @param nb_splits A number to indicate how many times the data should be split
#' @return A data frame with split seeds for each trial
search_cv_seeds_all_trials <- function(data,
                                       outcome_code, 
                                       covariate_code, 
                                       arm_code, 
                                       threshold = 0.05,
                                       nb_splits = 100) {
  trials <- unique(data$trial)
  seeds_list <- vector("list", length(trials))
  
  for (i_trial in 1:length(trials)) {
    
    data_current <- data |> filter(trial == trials[i_trial])
    
    seeds_list[[i_trial]] <- search_cv_seeds_one_trial(
      data_trial = data_current, 
      trial_code = trials[i_trial],
      outcome_code = outcome_code, 
      covariate_code = covariate_code, 
      arm_code = arm_code, 
      threshold = threshold, 
      nb_splits = nb_splits
    )
  }
  seeds <- do.call(cbind, seeds_list)
  
  return(seeds)
}


#' Create a grid with values of features of interest at which partial 
#' dependence will be calculated
#'
#' @param data_pd A data frame with trial data, should contain columns with the 
#'   features of interest
#' @param pd_feature_codes A vector with the names of features of interest
#' @param nb_steps The number of steps in the range of a feature at which 
#'   partial dependence will be calculated
#' @return A data frame that contains a grid with values of the features at 
#'   which partial dependence will be calculated
create_partial_dependence_grid <- function(data_pd, 
                                           pd_feature_codes, 
                                           nb_steps) {
  
  # values for the features of interest measured during the trial
  original_features <- na.omit(subset(data_pd, select=pd_feature_codes))
  
  # create an empty grid
  pd_grid_list <- vector("list", length(pd_feature_codes))
  
  for (i in 1:ncol(original_features)) {
    # get a min-max range for the current feature, exclude strong outliers
    current_feature <- deframe(original_features[i])
    current_feature <- Winsorize(
      x = current_feature, 
      val = quantile(current_feature, probs = c(0, 0.99))
    )
    min_value <- min(current_feature)
    max_value <- max(current_feature)
    
    # grid with feature values at which partial dependence will be calculated
    pd_grid_feature <- data.frame(
      step = seq(1, nb_steps+1, 1),
      feature_code = rep(colnames(original_features)[i], nb_steps+1),
      grid_value = seq(min_value, max_value, (max_value-min_value)/nb_steps)
    )
    pd_grid_list[[i]] <- pd_grid_feature
  }
  pd_grid <- do.call(rbind, pd_grid_list)
  
  return(pd_grid)
}

#' Split the data into a train and validation set based on predetermined seeds
#' This split is stratified by treatment arm and by trial
#' The predetermined seeds should ensure that average treatment effect is equal 
#' in train and validation set
#'
#' @param data_to_split A data frame with data to split, should contain a 
#'   "trial" column that identifies each trial
#' @param split_seeds A data frame with seeds for splitting, one column per 
#'   trial, one row per cross-validation fold
#' @param current_split A number to indicate which row of seeds should be used 
#'   in the current cross-validation fold
#' @return A list that contains a train and validation data frame
split_on_seed <- function(data_to_split, split_seeds, current_split) {
  
  # trials present in the data
  trials <- unique(data_to_split$trial)
  
  # create empty lists for train and validation sets
  train_list <- list()
  validation_list <- list()
  
  # split per trial into a stratified train and validation set
  for (i_trial in 1:length(trials)) {
    data_trial <- data_to_split |> filter(trial == trials[i_trial])
    idx <- splitTools::partition(
      y = pull(data_trial, arm_fac), 
      p = c(0.6, 0.4), 
      seed = pull(split_seeds, trials[i_trial])[[current_split]]
    )
    train_list[[i_trial]] <- data_trial[idx[[1]], ]
    validation_list[[i_trial]] <- data_trial[idx[[2]], ]
  }
  
  # combine train and validation sets into one large train and validation set
  train <- do.call(rbind, train_list)
  validation <- do.call(rbind, validation_list)
  
  return(list("train"=train, "validation"=validation))
}

#' Save estimates of individualized treatment effect for an evaluation set
#' Format results so that this is easily saved during cross-validation
#'
#' @param evaluation_ids A vector with IDs for the evaluation set
#' @param model_estimates A vector with model estimates
#' @param current_split A number to indicate the current cross-validation fold
#' @param current_mice A number to indicate the current imputation
#' @return A data frame with IDs and corresponding model estimates
save_model_estimates <- function(model_estimates, 
                                 evaluation_ids, 
                                 current_split, 
                                 current_mice) {
  
  # format the result
  df_estimates <- data.frame(
    split = current_split,
    mice = current_mice,
    USUBJID = evaluation_ids,
    tau_hat = model_estimates
  )
  return(df_estimates)
}

#' Calculate feature importance from the causal forest 
#' Format results so that this is easily saved during cross-validation
#'
#' @param causal_forest A trained causal forest to evaluate
#' @param feature_names A vector with names of features used for training
#' @param current_split A number to indicate the current cross-validation fold
#' @param current_mice A number to indicate the current imputation
#' @return A data frame with feature importances
calculate_feature_importance_cf <- function(causal_forest, 
                                            feature_codes, 
                                            current_split, 
                                            current_mice) {
  
  # extract feature importance from the trained causal forest
  importances_cf <- variable_importance(causal_forest)
  
  # format the result
  values <- as.list(c(current_split, current_mice, importances_cf))
  col_names <- c("split", "mice", feature_codes)
  df_importances <- as.data.frame(setNames(values, col_names))
  
  return(df_importances)
}

#' Calculate performance on an evaluation set 
#' With the area under the operator characteristic curve (AUTOC)
#' Format results so that this is easily saved during cross-validation
#'
#' @param prioritization_rule A vector with prioritizations to evaluate
#' @param evaluation_features A data frame with features for the evaluation set
#' @param evaluation_outcome A vector with the outcome for the evaluation set
#' @param evaluation_arm  A vector with the treatment arm for the evaluation set
#' @param current_split A number to indicate the current cross-validation fold
#' @param current_mice A number to indicate the current imputation
#' @param skip_calculation Skip calculation if TRUE, default is FALSE
#' @return A data frame with AUTOC performance values for the evaluation set
calculate_autoc_performance <- function(prioritization_rule, 
                                        evaluation_features, 
                                        evaluation_outcome, 
                                        evaluation_arm, 
                                        current_split, 
                                        current_mice,
                                        skip_calculation = FALSE) {
  
  # bypass AUTOC performance calculation
  # to save calculation time when working on other parts of the workflow
  if (skip_calculation) {
    # format the result
    df_performance <- data.frame(
      split = current_split,
      mice = current_mice,
      autoc_est = NA,
      autoc_std = NA,
      autoc_lci = NA,
      autoc_uci = NA,
      autoc_p = NA
    )
  }
  
  else {
    # evaluation set causal forest needed for AUTOC calculation
    cf_ev <- grf::causal_forest(
      X = evaluation_features, 
      Y = evaluation_outcome, 
      W = evaluation_arm, 
      seed = 0
    )
    
    # calculate AUTOC
    set.seed(0)
    autoc <- rank_average_treatment_effect(
      forest = cf_ev, 
      priorities = prioritization_rule, 
      q = seq(0.01, 1, by = 0.01)
    )
    autoc_est <- autoc$estimate
    autoc_std <- autoc$std.err
    
    # format the result
    df_performance <- data.frame(
      split = current_split,
      mice = current_mice,
      autoc_est = round(autoc_est, 4),
      autoc_std = round(autoc_std, 4),
      autoc_lci = round(autoc_est - 1.96*autoc_std, 4),
      autoc_uci = round(autoc_est + 1.96*autoc_std, 4),
      autoc_p = 2*pnorm(-abs(autoc_est/autoc_std))
    )
  }
  return(df_performance)
}

#' Calculate partial dependence for a causal forest
#' For specified model features and feature values
#' Format results so that this is easily saved during cross-validation
#'
#' @param causal_forest A trained causal forest to evaluate
#' @param evaluation_features A data frame with features for the evaluation set
#' @param pd_grid A data frame with features and feature values to evaluate at
#' @param current_split A number to indicate the current cross-validation fold
#' @param current_mice A number to indicate the current imputation
#' @param skip_calculation Skip calculation if TRUE, default is FALSE
#' @return A data frame with partial dependence results
calculate_partial_dependence_cf <- function(causal_forest, 
                                            evaluation_features, 
                                            pd_grid, 
                                            current_split, 
                                            current_mice,
                                            skip_calculation = FALSE) {
  
  # bypass partial dependence calculation
  # to save calculation time when working on other parts of the workflow
  if (skip_calculation) {
    # format the result
    pd_grid <- pd_grid |>
      mutate(mean_value = NA) |>
      mutate(mice = current_mice) |>
      relocate(mice) |>
      mutate(split = current_split) |>
      relocate(split)
  }
  
  else {
    # all features for which partial dependence should be calculated
    pd_feature_codes <- unique(pd_grid$feature_code)
    
    # add an empty row for results
    pd_grid <- pd_grid |> mutate(mean_value = NA)
    
    # for all features of interest
    for (i in 1:length(pd_feature_codes)) {
      
      # for all feature values of interest
      for (j in 1:max(pd_grid$step)) {
        
        # for all patients, set the value of the current feature to the current 
        # value of interest
        current_grid_value <- pd_grid |> 
          filter(feature_code == pd_feature_codes[i] & step == j) |>
          pull(grid_value)
        current_X <- mutate(evaluation_features, !!pd_feature_codes[i]:=current_grid_value)
        
        # calculate mean estimated treatment effect over all patients
        current_mean = mean(predict(causal_forest, newdata = current_X)$predictions)
        
        # save that mean for the current feature and current feature value of interest
        pd_grid <- pd_grid |>
          mutate(mean_value = ifelse(
            feature_code == pd_feature_codes[i] & step == j, 
            current_mean, 
            mean_value)
          )
      }
    }
    pd_grid <- pd_grid |> 
      mutate(mice = current_mice) |>
      relocate(mice) |>
      mutate(split = current_split) |> 
      relocate(split)
  }
  return(pd_grid)
}

#' Calculate median performance from cross-validation results
#'
#' @param cv_performance A data frame with performance from cross-validation
#' @return A data frame with median performance
calculate_median_performance <- function(cv_performance) {
  
  cols_performance <- colnames(cv_performance)
  cols_metrics <- cols_performance[!cols_performance %in% c("split", "mice")]
  
  df_median_performance <- cv_performance |>
    # calculate median performance over all cross-validation folds
    group_by(mice) |>
    summarise(across(all_of(cols_metrics), median)) |>
    # calculate mean performance over all imputations
    ungroup() |>
    summarise(across(all_of(cols_metrics), mean))
  
  return(df_median_performance)
}


#' Calculate median feature importance from cross-validation results
#' Format results for easy visualization
#'
#' @param cv_feature_importances A data frame with feature importances from 
#'   cross-validation
#' @param feature_map A data frame that maps feature codes to feature texts
#' @return A data frame with median importances formatted for easy visualization
calculate_median_feature_importance <- function(cv_feature_importances, 
                                                feature_map) {
  
  df_median_importances <- cv_feature_importances |>
    pivot_longer(
      cols = all_of(model_features), 
      names_to = "feature_code", 
      values_to = "value") |>
    # calculate median importance over all cross-validation folds
    group_by(feature_code, mice) |>
    summarise(
      median = median(value), 
      lower_quartile = quantile(value, probs = 0.25), 
      upper_quartile = quantile(value, probs = 0.75), 
      .groups = "drop_last") |>
    # calculate mean importance over all imputations
    group_by(feature_code) |>
    summarise(across(everything(), mean)) |>
    # map variable codes to full text definitions
    mutate(feature_text = recode_values(
      feature_code, 
      from = feature_map$feature_code, 
      to = feature_map$feature_text)) |>
    arrange(desc(median))
  
  return(df_median_importances)
}

#' Calculate average observed treatment effect with an analysis of covariance 
#' (ANCOVA) model
#'
#' @param data A data frame with trial data, should contain columns with the 
#'   outcome, covariate, and treatment arm of interest
#' @param outcome_code A string indicating the column name of the outcome
#' @param covariate_code A string indicating the column name of the covariate
#' @param arm_code A string indicating the column name of the treatment arm
#' @return A data frame with the estimated average observed treatment effect
calculate_average_effect <- function(data, outcome_code, covariate_code, arm_code) {
  
  formula <- as.formula(paste(outcome_code, " ~ ", covariate_code, " + ", arm_code))
  
  # analysis of covariance model to calculate the average observed effect
  ancova_fit   <- aov(formula, data = data)
  effect_means <- emmeans::emmeans(ancova_fit, arm_code)
  effect_means_ci <- summary(confint(pairs(effect_means, reverse = TRUE)))
  
  average_effect <- data.frame(
    estimate = effect_means_ci$estimate,
    lower_ci = effect_means_ci$lower.CL,
    upper_ci = effect_means_ci$upper.CL,
    std = summary(pairs(effect_means, reverse = TRUE))$SE
  )
}

#' Add quantile details to an estimate of average observed treatment effect
#' To be used later on for visualizing results
#'
#' @param average_treatment_estimate A data frame containing the estimate
#' @param prioritization_rule A vector with a value of the prioritization rule 
#'  for each person
#' @param quantile A number indicating the current quantile
#' @param quantile_text A string indicating the current quantile
#' @param current_split A number indicating the current cross-validation split
#' @param current_mice A number indicating the current MICE iteration
#' @return A data frame with estimates of average observed treatment effect
add_quantile_details <- function(average_effect_estimate, 
                                 prioritization_rule, 
                                 quantile, 
                                 quantile_text, 
                                 current_split = 0,
                                 current_mice = 0) {
  
  to_add <- data.frame(
    quantile = quantile,
    quantile_text = quantile_text,
    prior_mean = mean(prioritization_rule, na.rm = TRUE),
    prior_std = sd(prioritization_rule, na.rm = TRUE),
    nb_people = length(prioritization_rule),
    split = current_split,
    mice = current_mice
  )
  
  quantile_effect <- cbind(average_effect_estimate, to_add)
}

#' Calculate average observed treatment effect per quantile of a given 
#' prioritization rule
#'
#' @param data A data frame with trial data, should contain columns with the 
#'   prioritization rule, outcome, covariate, and treatment arm of interest
#' @param prior_code A string indicating the prioritization rule
#' @param outcome_code A string indicating the outcome
#' @param covariate_code A string indicating the covariate
#' @param arm_code A string indicating the treatment arm
#' @param nb_quantiles A number indicating the number of quantiles to use
#' @param current_split A number indicating the current cross-validation split
#' @param current_mice A number indicating the current MICE iteration
#' @param reverse A boolean indicating if the order of quantiles should be reversed
#' @return A data frame with estimates of average observed treatment effect
calculate_quantile_effects <- function(data, 
                                       prior_code, 
                                       outcome_code, 
                                       covariate_code, 
                                       arm_code, 
                                       nb_quantiles, 
                                       current_split = 0,
                                       current_mice = 0,
                                       reverse = FALSE) {
  
  # split into quantiles based on a prioritization rule
  data_split_quantiles <- data |>
    mutate(quantile = ntile(pull(data, prior_code), nb_quantiles))
  
  # reverse the order of quantiles if needed
  if (reverse) {
    reverse_quantiles <- tibble(
      before = sort(unique(data_split_quantiles$quantile)), 
      after = sort(unique(data_split_quantiles$quantile), decreasing = TRUE)
    )
    data_split_quantiles <- data_split_quantiles |>
      mutate(quantile = recode_values(
        quantile, from = reverse_quantiles$before, 
        to = reverse_quantiles$after)
      )
  }
  
  quantile_effects_list <- vector("list", nb_quantiles + 1)
  
  # overall average observed treatment effect
  average_effect <- calculate_average_effect(
    data = data_split_quantiles, 
    outcome_code = outcome_code, 
    covariate_code = covariate_code, 
    arm_code = arm_code
  )
  quantile_effects_list[[1]] <- add_quantile_details(
    average_effect_estimate = average_effect, 
    prioritization_rule = pull(data_split_quantiles, prior_code), 
    quantile = 0, 
    quantile_text = "Overall",
    current_split = current_split,
    current_mice = current_mice
  )
  
  # average observed treatment effect per quantile
  for (i_quantile in 1:nb_quantiles) {
    data_quantile <- data_split_quantiles |> filter(quantile == i_quantile)
    average_effect <- calculate_average_effect(
      data = data_quantile, 
      outcome_code = outcome_code, 
      covariate_code = covariate_code, 
      arm_code = arm_code
    )
    quantile_effects_list[[i_quantile + 1]] <- add_quantile_details(
      average_effect_estimate = average_effect, 
      prioritization_rule = pull(data_quantile, prior_code), 
      quantile = i_quantile, 
      quantile_text = paste("Q", i_quantile, sep=""),
      current_split = current_split,
      current_mice = current_mice
    )
  }
  quantile_effects <- do.call(rbind, quantile_effects_list)
  
  return(quantile_effects)
}

#' Calculate median observed average treatment effect (ATE) per quantile of a 
#' given prioritization rule over all cross-validation results
#'
#' @param all_trial_data A data frame with all trial data, should contain 
#'   columns with IDs, outcome, covariate, and treatment arm of interest that 
#'   are needed for calculating average treatment effects
#' @param cv_prior_estimates A data frame with IDs and the prioritization rule 
#'   values that were estimated in each cross-validation validation set
#' @param cv_prior_code A string indicating the prioritization rule of interest
#' @param outcome_code A string indicating the outcome of interest
#' @param covariate_code A string indicating the covariate for ATE correction
#' @param arm_code A string indicating the treatment arm in factor format
#' @param nb_quantiles A number indicating the number of quantiles to use
#' @return A data frame with a median ATE estimates per quantile
calculate_median_quantile_effects <- function(all_trial_data, 
                                              cv_prior_estimates,
                                              prior_code,
                                              outcome_code, 
                                              covariate_code, 
                                              arm_code, 
                                              nb_quantiles) {
  
  # additional information needed for ATE calculation
  to_add <- all_trial_data |> 
    select("USUBJID", outcome_code, covariate_code, arm_code)
  
  estimates_extended <- cv_prior_estimates |>
    merge(to_add, by="USUBJID")
  
  # calculate quantile effects in each validation set
  cv_q_effects_list <- vector("list", max(estimates_extended$split))
  
  for (i_split in 1:max(estimates_extended$split)) {
    
    cat("split seed", i_split, "\n")
    
    # calculate quantile effects in each imputation
    cv_q_effects_list[[i_split]] <- vector("list", max(estimates_extended$mice))
    
    for (i_mice in 1:max(estimates_extended$mice)) {
      
      evaluation_data <- estimates_extended |> 
        filter(split == i_split & mice == i_mice)
      
      cv_q_effects_list[[i_split]][[i_mice]] <- calculate_quantile_effects(
        data = evaluation_data, 
        prior_code = prior_code, 
        outcome_code = outcome_code,
        covariate_code = covariate_code,
        arm_code = arm_code,
        nb_quantiles = nb_quantiles,
        current_split = i_split,
        current_mice = i_mice
      )
    }
  }
  
  cv_q_effects <- do.call(rbind, unlist(cv_q_effects_list, recursive = FALSE))
  
  # calculate median quantile effects over all validation sets
  cols_q <- colnames(cv_q_effects)
  cols_median <- cols_q[!cols_q %in% c("split", "mice", "quantile_text")]
  
  df_median_quantile_effects <- cv_q_effects |>
    # calculate median effects over all cross-validation folds
    group_by(quantile_text, mice) |>
    summarise(across(all_of(cols_median), median), .groups = "drop_last") |>
    # calculate mean importance over all imputations
    select(-mice) |>
    group_by(quantile_text) |>
    summarise(across(everything(), mean)) |>
    # show a real number for the number of people
    mutate(nb_people = floor(nb_people))
  
  return(df_median_quantile_effects)
}

#' Calculate median partial dependence from cross-validation results
#'
#' @param cv_partial_dependence A data frame with partial dependence from 
#'   cross-validation
#' @return A data frame with median partial dependence
calculate_median_partial_dependence <- function(cv_partial_dependence) {
  
  df_median_partial_dependence <- cv_partial_dependence |>
    # calculate median partial dependence over all cross-validation folds
    select(-split) |>
    group_by(mice, step, grid_value, feature_code) |>
    summarise(pd_median = median(mean_value), 
              pd_lower_quartile = quantile(mean_value, probs = 0.25), 
              pd_upper_quartile = quantile(mean_value, probs = 0.75), 
              .groups = "drop_last") |>
    # calculate mean partial dependence over all imputations
    ungroup() |>
    select(-mice) |>
    group_by(step, grid_value, feature_code) |>
    summarise(across(everything(), mean), .groups = "drop_last")
  
  return(df_median_partial_dependence)
}

#' Calculate baseline characteristics per quantile for all cross-validation 
#' results
#'
#' @param all_trial_data A data frame data from all trials
#' @param cv_prior_estimates A data frame with IDs and the prioritization rule 
#'   values that were estimated in each cross-validation validation set
#' @param arm_code A string indicating the treatment arm
#' @param prior_code A string indicating the prioritization rule of interest
#' @param con_feature_codes A list of strings indicating the baseline 
#'   characteristics that should be evaluated. Currenly this function only works 
#'   for continuous characterstics
#' @param nb_quantiles The number of quantiles to use
#' @param feature_map A data frame that maps feature codes to feature texts
#' @return A data frame with baseline characteristics per quantile
calculate_cv_quantile_features <- function(all_trial_data, 
                                           cv_prior_estimates,
                                           arm_code,
                                           prior_code,
                                           con_feature_codes,
                                           nb_quantiles, 
                                           feature_map) {
  
  # additional features needed for tableone calculation
  to_add <- all_trial_data |>
    select("USUBJID", arm_code, all_of(con_feature_codes))
  
  estimates_extended <- cv_prior_estimates |>
    merge(to_add, by="USUBJID") |>
    group_by(split, mice) |>
    mutate(quantile =  ntile(tau_hat, nb_quantiles)) |>
    relocate(quantile) |>
    ungroup()
  
  one_con_list <- vector("list", max(estimates_extended$split))
  
  for (i_split in 1:max(estimates_extended$split)) {
    
    one_con_list[[i_split]] <- vector("list", max(estimates_extended$mice))
    
    for (i_mice in 1:max(estimates_extended$mice)) {
      
      evaluation_data <- estimates_extended |> 
        filter(split == i_split & mice == i_mice)
      
      nb_people <- evaluation_data |>
        group_by(quantile) |>
        summarize(n = n(), .groups = "drop_last") |>
        rename(nb = n) |>
        mutate(feature_code = "n") |>
        mutate(feature_text = "Number of people") |>
        mutate(me = NA) |>
        mutate(lq = NA) |>
        mutate(uq = NA) |>
        mutate(split = i_split) |>
        mutate(mice = i_mice)
      
      df_one_con <- evaluation_data |>
        select(all_of(con_features), "quantile") |>
        group_by(quantile) |>
        summarise(across(everything(), list(
          me = ~ median(.x, na.rm = TRUE), 
          lq = ~ quantile(.x, probs = 0.25, na.rm = TRUE), 
          uq = ~ quantile(.x, probs = 0.75, na.rm = TRUE))),
          .groups = "drop_last"
        ) |>
        pivot_longer(
          cols = -c(quantile), 
          names_to = c("feature_code", ".value"), 
          names_pattern = "(.*)_(me|lq|uq)"
        ) |>
        mutate(split = i_split) |>
        mutate(mice = i_mice) |>
        mutate(feature_text = recode_values(
          feature_code, 
          from = feature_map$feature_code, 
          to = feature_map$feature_text)
        ) |>
        relocate(feature_text, .after = feature_code) |>
        mutate(nb = NA)
      
      df_one_con <- rbind(df_one_con, nb_people)
      
      one_con_list[[i_split]][[i_mice]] <- df_one_con
    }
  }
  one_con_all <- do.call(rbind, unlist(one_con_list, recursive = FALSE))
  
  return(one_con_all)
}

#' Calculate median baseline characteristics per quantile
#'
#' @param cv_quantile_features A data frame with baseline characteristics per 
#'   quantile from cross-validation
#' @return A data frame with median baseline characteristics per quantile
calculate_median_quantile_features <- function(cv_quantile_features) {
  
  median_features <- cv_quantile_features |>
    select(-split) |>
    # calculate median effects over all cross-validation folds
    group_by(feature_code, feature_text, mice, quantile) |>
    summarise(across(everything(), median), .groups = "drop_last") |>
    ungroup() |>
    # calculate mean importance over all imputations
    select(-mice) |>
    group_by(feature_code, feature_text,  quantile) |>
    summarise(across(everything(), mean), .groups = "drop_last") |>
    # show a real number for the number of people
    mutate(nb = floor(nb)) |>
    mutate(across(c(me, lq, uq), ~ sprintf("%.2f", .x))) |>
    # reformat the results
    mutate(summary = NA) |>
    mutate(summary = ifelse(
      !feature_code == "n", 
      paste0(me, " [", lq, ", ", uq, "]"), 
      summary
    )) |>
    mutate(summary = ifelse(
      feature_code == "n", 
      paste0(nb), 
      summary
    )) |>
    select(-c(me, lq, uq, nb))
  
  return(median_features)
}
