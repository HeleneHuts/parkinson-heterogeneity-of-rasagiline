#' Preprocessing of NCT02337725, a Japanese trial in patients with early 
#' Parkinson's disease
#' 
#' Input: Data from the Integrated Parkinson’s database of the Critical Path 
#'   Institute in SDTM format
#'
#' Output: A table with one row per patient and one column per baseline 
#'   characteristic and per outcome

# Needed libraries and functions -----------------------------------------------

library(tidyverse)

source("analysis/utils/utils_preprocessing_NCT02337725.R")
source("analysis/utils/utils_preprocessing_Japanese_trials.R")

# Load files -------------------------------------------------------------------

# raw data from the trial
lb_all <- read.csv("data/raw/lb.csv")
vs_all <- read.csv("data/raw/vs.csv")
xd_all <- read.csv("data/raw/xd.csv")
qs_all <- read.csv("data/raw/qs.csv")
dm_all <- read.csv("data/raw/dm.csv")
eg_all <- read.csv("data/raw/eg.csv")
cm_all <- read.csv("data/raw/cm.csv")

# maps for concomitant medication and UPDRS questions
med_map <- read.csv("analysis/maps/medication_map_manuscript.csv")
updrs_map <- read.csv("analysis/maps/updrs_map_manuscript.csv")

# Extract NCT02337725 ----------------------------------------------------------

trial_id <- get_trial_id_NCT02337725(dm_all)

lb_raw <- lb_all |> filter(STUDYID == trial_id)
vs_raw <- vs_all |> filter(STUDYID == trial_id)
xd_raw <- xd_all |> filter(STUDYID == trial_id)
qs_raw <- qs_all |> filter(STUDYID == trial_id)
dm_raw <- dm_all |> filter(STUDYID == trial_id)
eg_raw <- eg_all |> filter(STUDYID == trial_id)
cm_raw <- cm_all |> filter(STUDYID == trial_id)


# Preprocessing of each separate file ------------------------------------------

dm <- preprocess_dm_NCT02337725(dm_raw)
lb <- preprocess_lb_NCT02337725(lb_raw)
vs <- preprocess_vs_Japanese_trial(vs_raw)
xd <- preprocess_xd_Japanese_trial(xd_raw)
eg <- preprocess_eg_Japanese_trial(eg_raw)
cm <- preprocess_cm_NCT02337725(cm_raw, dm_raw, med_map)
qs_other <- preprocess_qs_other_NCT02337725(qs_raw, updrs_map)
qs_updrs <- preprocess_qs_updrs_NCT02337725(qs_raw, updrs_map)


# Combine the preprocessed files -----------------------------------------------

preprocessed_data <- dm |>
  left_join(qs_updrs, by = "USUBJID") |>
  left_join(qs_other, by = "USUBJID") |>
  left_join(cm, by = "USUBJID") |>
  left_join(lb, by = "USUBJID") |>
  left_join(vs, by = "USUBJID") |>
  left_join(xd, by = "USUBJID") |>
  left_join(eg, by = "USUBJID") |>
  mutate(trial = "NCT02337725")

# Save results -----------------------------------------------------------------

if (!file.exists("data/preprocessed")) {
  dir.create("data/preprocessed")
}

write.csv(preprocessed_data, "data/preprocessed/NCT02337725_preprocessed.csv")
