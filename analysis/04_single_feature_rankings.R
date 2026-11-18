#' Single feature rankings
#' 
#' Loads and pools preprocessed data from 
#'  NCT00256204 - ADAGIO, trial in early Parkinson's disease
#'  NCT02337725 - Japanese trial in early Parkinson's disease
#'  NCT02337738 - Japanese trial in advanced Parkinson's disease
#'  
#' Investigates quintiles of observed average treatment effect for rankings 
#' based on single important features

# libraries --------------------------------------------------------------------

library(tidyverse)
library(ggplot2)

source("analysis/utils/utils_model.R")
source("analysis/utils/utils_plot.R")


# load all data ----------------------------------------------------------------

data_NCT00256204_o <- read.csv("data/preprocessed/NCT00256204_preprocessed.csv")
data_NCT02337725_o <- read.csv("data/preprocessed/NCT02337725_preprocessed.csv")
data_NCT02337738_o <- read.csv("data/preprocessed/NCT02337738_preprocessed.csv")

feature_map <- read.csv("analysis/maps/feature_map_manuscript.csv")


# directories for saving results -----------------------------------------------

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


# feature codes ----------------------------------------------------------

# single features to investigate
single_features <- c("MDS_UPDRS2_total_w0_o", "MDS_UPDRS3_total_w0_o", "GGT", 
                     "WEIGHT", "QTCBAG", "QRSAG")

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
  select(c("USUBJID", "trial", all_of(single_features), all_of(arm_num), 
           all_of(arm_fac), all_of(outcome), all_of(covariate))) |>
  filter(!is.na(.data[[outcome]]))

data_NCT02337725 <- data_NCT02337725_o |>
  select(c("USUBJID", "trial", all_of(single_features), all_of(arm_num), 
           all_of(arm_fac), all_of(outcome), all_of(covariate))) |>
  filter(!is.na(.data[[outcome]]))

data_NCT02337738 <- data_NCT02337738_o |>
  select(c("USUBJID", "trial", all_of(single_features), all_of(arm_num), 
           all_of(arm_fac), all_of(outcome), all_of(covariate))) |>
  filter(!is.na(.data[[outcome]]))

data <- do.call(rbind, list(data_NCT00256204, data_NCT02337725, data_NCT02337738))


# create supplementary figures 1-6 ---------------------------------------------

for (i_feature in 1:length(single_features)) {
  
  feature <- single_features[i_feature]
  
  x_label <- feature_map |> 
    filter(feature_code == feature) |> 
    pull(feature_text_unit)
  
  # ranking in all patients
  trial = "All patients"
  save_name = paste0("quintiles_", feature, "_all")
  
  all_quintiles <- calculate_quantile_effects(
    data = data,
    prior_code = feature,
    outcome_code = outcome,
    covariate_code = covariate,
    arm_code = arm_fac,
    nb_quantiles = 5,
    reverse = FALSE)
  
  # create and save plot
  plot <- plot_quantile_effects_single(
    all_quintiles, 
    title = paste0(trial, "\nMDS-UPDRS part 2+3 at week 14"), 
    x_label = x_label
  )
  ggsave(
    filename = paste0(results_plot_path, save_name, ".tiff"), 
    plot = plot, 
    width = 4.5, 
    height = 4.5,
    dpi = 100
  )
  
  trials <- unique(data$trial)
  
  for (i_trial in 1: length(trials)) {
    
    # ranking in separate trials
    current_trial = trials[i_trial]
    save_name = paste0("quintiles_", feature, "_", current_trial)
    
    trial_data <- data |>
      filter(trial == current_trial)
    
    trial_quintiles <- calculate_quantile_effects(
      data = trial_data,
      prior_code = feature,
      outcome_code = outcome,
      covariate_code = covariate,
      arm_code = arm_fac,
      nb_quantiles = 5,
      reverse = FALSE)
    
    # create and save plot
    plot <- plot_quantile_effects_single(
      trial_quintiles, 
      title = paste0(current_trial, "\nMDS-UPDRS part 2+3 at week 14"), 
      x_label = x_label
    )
    ggsave(
      filename = paste0(results_plot_path, save_name, ".tiff"), 
      plot = plot, 
      width = 4.5, 
      height = 4.5,
      dpi = 100
    )
  }
}

