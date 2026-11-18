#' Exploratory data analysis
#' 
#' Loads and pools preprocessed data from 
#'  NCT00256204 - ADAGIO, trial in early Parkinson's disease
#'  NCT02337725 - Japanese trial in early Parkinson's disease
#'  NCT02337738 - Japanese trial in advanced Parkinson's disease
#'  
#' Investigates summary statistics of baseline characteristics and visualizes 
#' average treatment effect per trial

# libraries --------------------------------------------------------------------

library(tidyverse)
library(ggplot2)
library(tableone)

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


# feature codes ----------------------------------------------------------------

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

# all last observation carried forward outcomes
all_locf <- c("MDS_UPDRS23_total_diff_w6_locf", "MDS_UPDRS23_total_diff_w10_locf",
              "MDS_UPDRS23_total_diff_w12_locf", "MDS_UPDRS23_total_diff_w14_locf",
              "MDS_UPDRS23_total_diff_w20_locf", "MDS_UPDRS23_total_diff_w24_locf",
              "MDS_UPDRS23_total_diff_w26_locf", "MDS_UPDRS23_total_diff_w36_locf")

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
           all_of(arm_fac), all_of(outcome), all_of(covariate),
           all_of(all_locf))) |>
  filter(!is.na(.data[[outcome]]))

data_NCT02337725 <- data_NCT02337725_o |>
  select(c("USUBJID", "trial", all_of(all_features), all_of(arm_num), 
           all_of(arm_fac), all_of(outcome), all_of(covariate),
           all_of(all_locf))) |>
  filter(!is.na(.data[[outcome]]))

data_NCT02337738 <- data_NCT02337738_o |>
  select(c("USUBJID", "trial", all_of(all_features), all_of(arm_num), 
           all_of(arm_fac), all_of(outcome), all_of(covariate),
           all_of(all_locf))) |>
  filter(!is.na(.data[[outcome]]))

data <- do.call(rbind, list(data_NCT00256204, data_NCT02337725, data_NCT02337738))


# test normal distribution -----------------------------------------------------

# none of the continuous characteristics have an overall normal distribution
normal <- data |>
  select(all_of(con_features)) |>
  summarise(across(everything(), \(x) stats::shapiro.test(x)$p.value))


# create table 1 ---------------------------------------------------------------
# table of baseline characteristics of patients in the trials

# table 1, per treatment arm
one_arm <- CreateTableOne(
  data = data, 
  strata = arm_fac, 
  vars = table_one_features, 
  factorVars = cat_features)

one_arm_print <- as.data.frame(print(
  x = one_arm, 
  quote = FALSE, 
  noSpaces = TRUE, 
  printToggle = FALSE, 
  nonnormal = con_features, 
  catDigits = 2, 
  contDigits = 2, 
  showAllLevels = TRUE)
)

# table 1, per trial
one_trial <- CreateTableOne(
  data = data, 
  strata = c("trial"), 
  vars = table_one_features, 
  factorVars = cat_features)

one_trial_print <- as.data.frame(print(
  x = one_trial, 
  quote = FALSE, 
  noSpaces = TRUE, 
  printToggle = FALSE, 
  nonnormal = con_features, 
  catDigits = 2, 
  contDigits = 2, 
  showAllLevels = TRUE)
)

# count the number of missing values per feature
missing <- data |>
  select(all_of(table_one_features)) |>
  summarise(across(everything(), ~ sum(is.na(.))))

# save results
write.csv(one_trial_print, paste0(results_tables_path, "02_one_trial.csv"), row.names=TRUE)
write.csv(one_arm_print, paste0(results_tables_path, "02_one_arm.csv"), row.names = TRUE)

# create figure 2---------------------------------------------------------------
# figure of average LOCF treatment effect at different visits in the trials


# calculate average effects at weeks of interest
weeks <- c(0, 6, 10, 12, 14, 20, 24, 26, 36)

average_effects_locf <- calculate_average_effects_locf(
  data = data, 
  outcome_base_code = "MDS_UPDRS23_total_diff",
  covariate_code = covariate,
  arm_code = arm_fac,
  weeks_of_interest = weeks)

# visualize average effects
plot <- plot_average_effects_locf(average_effects_locf)

ggsave(
  filename = paste0(results_plot_path, "average_effects_locf", ".tiff"), 
  plot = plot, 
  width = 6.5, 
  height = 4,
  dpi = 200
)

