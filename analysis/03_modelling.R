#' Main modelling file
#' 
#' Loads and pools preprocessed data from 
#'  NCT00256204 - ADAGIO, trial in early Parkinson's disease
#'  NCT02337725 - Japanese trial in early Parkinson's disease
#'  NCT02337738 - Japanese trial in advanced Parkinson's disease
#'  
#' A causal forest is trained and evaluated by using a Monte Carlo 
#' cross-validation approach with 100 training and internal validation splits

# libraries --------------------------------------------------------------------

library(DescTools)
library(mice)
library(tidyverse)
library(grf)
library(tictoc)
library(splitTools)
library(ggplot2)
library(emmeans)


# load all data ----------------------------------------------------------------

data_NCT00256204_o <- read.csv("data/preprocessed/NCT00256204_preprocessed.csv")
data_NCT02337725_o <- read.csv("data/preprocessed/NCT02337725_preprocessed.csv")
data_NCT02337738_o <- read.csv("data/preprocessed/NCT02337738_preprocessed.csv")

feature_map <- read.csv("analysis/maps/feature_map_manuscript.csv")

source("analysis/utils/utils_model.R")
source("analysis/utils/utils_plot.R")


# directories for saving results -----------------------------------------------

if (!file.exists("analysis/seeds")) {
  dir.create("analysis/seeds")
}

if (!file.exists("results")) {
  dir.create("results")
}

if (!file.exists("results/tables")) {
  dir.create("results/tables")
}

if (!file.exists("results/plots")) {
  dir.create("results/plots")
}

results_tables_path <- "results/tables/"
results_plot_path <- "results/plots/"


# model setup details ----------------------------------------------------------

# determine if computationally heavy parts should be skipped
skip_performance_flag <- FALSE
skip_partial_dependence_flag <- FALSE

# reusing seeds only works if a set of seeds was calculated in a previous run
reuse_seeds <- FALSE

# determine the number and iterations for MICE
# the manuscript used 10 imputations and 5 iterations - note that this will take 
# hours to compute, results with less imputations are still very representative
mice_nb_imputations <- 10
mice_max_iterations <- 5


# model feature codes ----------------------------------------------------------

# all baseline characteristics
all_features <- feature_map |> 
  filter(all_flag == "all") |> 
  pull(feature_code)

# baseline characteristics for table one
table_one_features <- feature_map |> 
  filter(table_one_flag == "table_one") |> 
  pull(feature_code)

# categorical baseline characteristics for table one
cat_features <- feature_map |> 
  filter(table_one_type == "categorical") |> 
  pull(feature_code)

# continuous baseline characteristics for table one
con_features <- feature_map |> 
  filter(table_one_type == "continuous") |> 
  pull(feature_code)

# baseline characteristics used to train the causal forests
model_features <- feature_map |> 
  filter(model_flag == "model") |> 
  pull(feature_code)

# model features for which partial dependence should be calculated
pd_features <- c("QRSAG", "QTCBAG", "WEIGHT", "MDS_UPDRS2_total_w0_o", 
                 "MDS_UPDRS3_total_w0_o", "GGT")

# numeric treatment indicator (1 = treatment, 0 = control)
arm_num <- "arm_num"

# factor treatment indicator (1 = treatment, 0 = control)
arm_fac <- "arm_fac"

# the outcome of interest
outcome <- "MDS_UPDRS23_total_diff_w14_i"

# covariate for correction in average treatment effect calculation
covariate <- "MDS_UPDRS23_total_w0_o"


# select data needed for the analysis ------------------------------------------
# exclude people who do not have any follow-up

data_NCT00256204 <- data_NCT00256204_o |>
  select(c("USUBJID", "trial", all_of(all_features), all_of(arm_num), 
           all_of(arm_fac), all_of(outcome), all_of(covariate))) |>
  filter(!is.na(.data[[outcome]]))

data_NCT02337725 <- data_NCT02337725_o |>
  select(c("USUBJID", "trial", all_of(all_features), all_of(arm_num), 
           all_of(arm_fac), all_of(outcome), all_of(covariate))) |>
  filter(!is.na(.data[[outcome]]))

data_NCT02337738 <- data_NCT02337738_o |>
  select(c("USUBJID", "trial", all_of(all_features), all_of(arm_num), 
           all_of(arm_fac), all_of(outcome), all_of(covariate))) |>
  filter(!is.na(.data[[outcome]]))

data <- do.call(
  rbind, 
  list(data_NCT00256204, data_NCT02337725, data_NCT02337738)
  ) |>
  # sort because split seeds depend on the order of patients
  arrange(USUBJID)


# search balanced splits -------------------------------------------------------

if (!reuse_seeds) {
  new_seeds <- search_cv_seeds_all_trials(
    data = data,
    outcome_code = outcome,
    arm_code = arm_fac,
    covariate_code = covariate,
    threshold = 0.05,
    nb_splits = 100)
  
  write.csv(new_seeds, "analysis/seeds/03_seeds.csv")
  
  split_seeds <- new_seeds
  
} else {
  
  split_seeds <- read.csv("analysis/seeds/03_seeds.csv")
}

# cross-validation for causal forest -------------------------------------------


# empty lists for saving results
estimates_list <- vector("list", nrow(split_seeds))
importances_list <- vector("list", nrow(split_seeds))
performance_list <- vector("list", nrow(split_seeds))
partial_dependence_list <- vector("list", nrow(split_seeds))

# determine at which feature values partial dependence should be calculated
partial_dependence_grid <- create_partial_dependence_grid(data, pd_features, 20)

tic()

for(i_split in 1:nrow(split_seeds)){
  
  cat("split seed", i_split, "\n")
  
  # split the data for Monte-Carlo cross-validation, based on predetermined seeds
  data_split <- split_on_seed(data, split_seeds, i_split)
  
  # nested empty lists for saving results
  estimates_list[[i_split]] <- vector("list", mice_nb_imputations)
  importances_list[[i_split]] <- vector("list", mice_nb_imputations)
  performance_list[[i_split]] <- vector("list", mice_nb_imputations)
  partial_dependence_list[[i_split]] <- vector("list", mice_nb_imputations)
  
  # impute missing baseline characteristics with MICE
  X_tr_mice <- mice(
    data = subset(data_split$train, select=model_features), 
    m = mice_nb_imputations, 
    maxit = mice_max_iterations,
    seed = 1, 
    printFlag = FALSE
  )
  X_va_mice <- mice(
    data = subset(data_split$validation, select=model_features), 
    m = mice_nb_imputations, 
    maxit = mice_max_iterations,
    seed = 1, 
    printFlag = FALSE
  )
  
  for (i_mice in 1:mice_nb_imputations) {
    
    # data from the current imputation
    X_tr <- mice::complete(X_tr_mice, i_mice)
    X_va <- mice::complete(X_va_mice, i_mice)
    
    # train a causal forest on the current train set
    cf_tr <- grf::causal_forest(
      X = X_tr, 
      Y = pull(data_split$train, outcome), 
      W = pull(data_split$train, arm_num), 
      seed = 1)
    
    # estimate individualized treatment effect on the current validation set
    tau_hat_va <- predict(cf_tr, X_va)$predictions
    
    # save validation set estimates
    estimates_list[[i_split]][[i_mice]] <- save_model_estimates(
      evaluation_ids = pull(data_split$validation, USUBJID),
      model_estimates = tau_hat_va,
      current_split = i_split,
      current_mice = i_mice
    )
    
    # calculate and save performance on the validation set
    performance_list[[i_split]][[i_mice]] <- calculate_autoc_performance(
      prioritization_rule = tau_hat_va, 
      evaluation_features = X_va,
      evaluation_outcome = pull(data_split$validation, outcome), 
      evaluation_arm = pull(data_split$validation, arm_num), 
      current_split = i_split,
      current_mice = i_mice,
      skip_calculation = skip_performance_flag
    )
    
    # calculate and save feature importance
    importances_list[[i_split]][[i_mice]] <- calculate_feature_importance_cf(
      causal_forest = cf_tr,
      feature_codes = model_features,
      current_split = i_split,
      current_mice = i_mice
    )
    
    # calculate and save partial dependence
    partial_dependence_list[[i_split]][[i_mice]] <- calculate_partial_dependence_cf(
      causal_forest = cf_tr, 
      evaluation_features = X_va, 
      pd_grid = partial_dependence_grid, 
      current_split = i_split,
      current_mice = i_mice,
      skip_calculation = skip_partial_dependence_flag
    )
  }
}

# save all cross-validation results in data frames
estimates <- do.call(rbind, unlist(estimates_list, recursive = FALSE))
importances <- do.call(rbind, unlist(importances_list, recursive = FALSE))
performance <- do.call(rbind, unlist(performance_list, recursive = FALSE))
partial_dependence <- do.call(rbind, unlist(partial_dependence_list, recursive = FALSE))

toc()

write.csv(estimates, paste0(results_tables_path, "03_all_estimates.csv"))
write.csv(importances, paste0(results_tables_path, "03_all_importances.csv"))
write.csv(performance, paste0(results_tables_path, "03_all_performance.csv"))
write.csv(partial_dependence, paste0(results_tables_path, "03_all_partial_dependence.csv"))


# evaluate performance with AUTOC ----------------------------------------------

median_performance <- calculate_median_performance(performance)


# evaluate feature importance --------------------------------------------------

median_importance <- calculate_median_feature_importance(importances, feature_map)

feature_importance_plot <- plot_median_feature_importance(median_importance, y_max = 0.11)

ggsave(
  filename = paste0(results_plot_path, "feature_importance", ".tiff"), 
  plot = feature_importance_plot, 
  width = 6, 
  height = 4,
  dpi = 200
)


# evaluate effect per quintile -------------------------------------------------

median_quantile_effects <- calculate_median_quantile_effects(
  all_trial_data = data, 
  cv_prior_estimates = estimates, 
  prior_code = "tau_hat",
  outcome_code = outcome,
  covariate_code = covariate,
  arm_code = arm_fac,
  nb_quantiles = 5
)

quantile_effects_plot <- plot_quantile_effects(
  median_quantile_effects, 
  title = "MDS-UPDRS part 2+3 at week 14", 
  x_label = "Estimated treatment effect"
)

ggsave(
  filename = paste0(results_plot_path, "quantile_effects", ".tiff"), 
  plot = quantile_effects_plot, 
  width = 4.5, 
  height = 4.5,
  dpi = 200
)

quantile_calibration_plot <- plot_quantile_calibration_curve(
  median_quantile_effects, 
  title = "Binned calibration curve"
)

ggsave(
  filename = paste0(results_plot_path, "quantile_calibration", ".tiff"), 
  plot = quantile_calibration_plot, 
  width = 4.5, 
  height = 4.5,
  dpi = 200
)


# evaluate features per quintile -----------------------------------------------

cv_quantiles <- calculate_cv_quantile_features(
  all_trial_data = data, 
  cv_prior_estimates = estimates, 
  arm_code = arm_fac,
  prior_code = "tau_hat",
  con_feature_codes = con_features, 
  nb_quantiles = 5, 
  feature_map = feature_map)

median_quantile_features <- calculate_median_quantile_features(cv_quantiles)

write.csv(median_quantile_features, paste0(results_tables_path, "03_median_quantile_features.csv"))

# evaluate partial dependence --------------------------------------------------

median_partial_dependence <- calculate_median_partial_dependence(partial_dependence)

pd_plot <- plot_partial_dependence(median_partial_dependence, median_importance, feature_map)

ggsave(
  filename = paste0(results_plot_path, "partial_dependence", ".tiff"), 
  plot = pd_plot, 
  width = 9, 
  height = 5,
  dpi = 200
)



