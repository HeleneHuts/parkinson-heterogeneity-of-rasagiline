#' Functions needed for preprocessing of NCT02337738, a Japanese trial in 
#' patients with advanced Parkinson's disease
#'
#' The functions in this file are specifically written for NCT02337738, due to 
#' differences in the structure of this trail compared to the other trials

source("analysis/utils/utils_preprocessing_Japanese_trials.R")

#' Get the trial ID of NCT02337738
#'
#' @param dm_all A data frame with raw demographics for all trials
#' @return A string ID that identifies NCT02337738
get_trial_id_NCT02337738 <- function(dm_all) {
  
  trial_arms <- c("ACTIVE ARM 1", "ACTIVE ARM 2", "Placebo")
  
  dm_id <- dm_all |>
    group_by(STUDYID) |>
    filter(all(trial_arms %in% ACTARM)) |>
    filter(sum(!is.na(ACTARM)) == 403) |>
    filter(sum(ACTARM == "ACTIVE ARM 1") == 129) |>
    filter(sum(ACTARM == "ACTIVE ARM 2") == 133) |>
    filter(sum(ACTARM == "Placebo") == 141) |>
    pull(STUDYID)
  
  trial_id <- dm_id[1]
  
  return(trial_id)
}

#' Preprocessing of features from demographics
#'
#' @param demograhics_raw A data frame with raw demographics for this trial
#' @return A data frame with preprocessed demographic features
preprocess_dm_NCT02337738 <- function(demographics_raw) {
  
  dm <- demographics_raw |>
    select(USUBJID, AGE, SEX, ACTARM, COUNTRY) |>
    mutate(arm = recode_values(
      ACTARM, 
      "Placebo" ~ "placebo", 
      "ACTIVE ARM 1" ~ "rasagiline", 
      "ACTIVE ARM 2" ~ "rasagiline"
    )) |>
    mutate(arm_d = recode_values(
      ACTARM, 
      "Placebo" ~ "placebo", 
      "ACTIVE ARM 1" ~ "rasagiline_1mg", 
      "ACTIVE ARM 2" ~ "rasagiline_05mg"
    )) |>
    mutate(arm_num = recode_values(arm, "placebo" ~ 0, "rasagiline" ~ 1)) |>
    mutate(arm_fac = as.factor(arm)) |>
    mutate(male_num = recode_values(SEX, "M" ~ 1, "F" ~ 0)) |>
    mutate(sex_fac = factor(SEX, levels = c("F", "M"))) |>
    # group countries into regions:  
    # 1 = Asia: JPN, 
    # 2 = Europe: AUT DEU ESP FRA GBR HUN ITA NLD PRT ROU, 
    # 3 = North America: CAN USA, 
    # 4 = Middle East: ISR, 
    # 5 = South America: ARG 
    mutate(region = recode_values(COUNTRY, "JPN" ~ "Asia")) |>
    mutate(region_num = recode_values(COUNTRY, "JPN" ~ 1)) |>
    mutate(region_fac = factor(
      region, 
      levels = c("Asia", "Europe", "North America", "Middle East", "South America")
    )) |>
    # one-hot encoding for region, 
    # largest group Europe will not be used for the model to avoid colinearity
    mutate(
      region_asia = ifelse(COUNTRY == "JPN", 1, 0),
      region_north_america = ifelse(COUNTRY %in% c("CAN", "USA"), 1, 0),
      region_middle_east = ifelse(COUNTRY == "ISR", 1, 0),
      region_south_america = ifelse(COUNTRY == "ARG", 1, 0)
    )
  return(dm)
}

#' Preprocessing of features from laboratory results
#'
#' @param laboratory_results_raw A data frame with raw laboratory results for 
#'   this trial
#' @return A data frame with preprocessed laboratory results features
preprocess_lb_NCT02337738 <- function(laboratory_results_raw) {
  
  # lab features used in the current analysis
  lb_features <- c("ALB", "ALT", "AST", "BASOLE", "BILI", "CA", "CHOL", "CL", 
                   "CREAT", "EOSLE", "GGT", "GLUC", "HGB", "K", "LDH", "LYMLE", 
                   "MONOLE", "NEUTLE", "PHOS", "PLAT", "PROT", "UREAN", "WBC")
  
  # values at screening
  lb_sc <- laboratory_results_raw |>
    select(USUBJID, LBTESTCD, LBSTRESC, LBCAT, VISIT) |>
    filter(VISIT=="VISIT 1(Week -2)") |>
    filter(LBTESTCD %in% lb_features) |>
    filter(!(LBTESTCD == "PROT" & LBCAT == "URINALYSIS"))  |>
    select(-c(LBCAT, VISIT)) |>
    rename(LBSTRESC_sc = LBSTRESC)
  
  # values at baseline
  lb_bl <- laboratory_results_raw |>
    select(USUBJID, LBTESTCD, LBSTRESC, LBCAT, VISIT) |>
    filter(VISIT=="VISIT 2(Week 0)" & LBSTRESC != "COAGULATION") |>
    filter(LBTESTCD %in% lb_features) |>
    filter(!(LBTESTCD == "PROT" & LBCAT == "URINALYSIS"))  |>
    select(-c(LBCAT, VISIT)) |>
    rename(LBSTRESC_bl = LBSTRESC)
  
  # use the values at screening to complete baseline if missing
  lb <- join_lb_Japanese_trial(lb_sc, lb_bl)
  
  return(lb)
}

#' Preprocessing of features from concomitant medication
#'
#' @param concomitant_medication_raw A data frame with raw concomitant 
#'   medication for this trial
#' @param demograhics_raw A data frame with raw demographics for this trial
#'   needed for a full list of USUBJIDs
#' @return A data frame with preprocessed concomitant medication features
preprocess_cm_NCT02337738 <- function(concomitant_medication_raw, 
                                      demographics_raw, 
                                      medication_map) {
  
  # get medication names per class
  doag <- medication_map |> filter(class == "Dopamine agonists") |> pull(medication)
  COMT <- medication_map |> filter(class == "COMT inhibitors") |> pull(medication)
  aman <- medication_map |> filter(class == "Amantadine") |> pull(medication)
  anti <- medication_map |> filter(class == "Anticholinergics") |> pull(medication)
  drox <- medication_map |> filter(class == "Droxidopa") |> pull(medication)
  istr <- medication_map |> filter(class == "Istradefylline") |> pull(medication)
  zoni <- medication_map |> filter(class == "Zonisamide") |> pull(medication)
  carb <- medication_map |> filter(class == "Carbidopa") |> pull(medication)
  levo <- medication_map |> filter(class == "Levodopa") |> pull(medication)
  
  
  # concomitant Parkinson medications
  cm_p <- concomitant_medication_raw |>
    filter(CMCAT == "CONCOMITANT MEDICATIONS") |>
    # dates that are NA indicate that medication was concomitant at baseline
    mutate(CMSTDY = ifelse(is.na(CMSTDY), -15, CMSTDY)) |>
    mutate(CMENDY = ifelse(is.na(CMENDY), 2, CMENDY)) |>
    filter(CMSTDY < -14 & CMENDY > 1) |>
    select(USUBJID, CMTRT) |>
    mutate(dopamine_agonists = ifelse(CMTRT %in% doag, 1, 0)) |>
    mutate(COMT_inhibitors = ifelse(CMTRT %in% COMT, 1, 0)) |>
    mutate(amantadine = ifelse(CMTRT %in% aman, 1, 0)) |>
    mutate(anticholinergics = ifelse(CMTRT %in% anti, 1, 0)) |>
    mutate(droxidopa = ifelse(CMTRT %in% drox, 1, 0)) |>
    mutate(istradefylline = ifelse(CMTRT %in% istr, 1, 0)) |>
    mutate(zonisamide = ifelse(CMTRT %in% zoni, 1, 0)) |>
    mutate(carbidopa = ifelse(CMTRT %in% carb, 1, 0)) |>
    mutate(levodopa = 1) |> # inclusion criterium for this trial
    group_by(USUBJID) |>
    summarize(across(where(is.numeric), max))
  
  # number of concomitant medications per person
  cm_n <- concomitant_medication_raw |>
    filter(CMSTRTPT == "BEFORE") |>
    filter(CMCAT == "CONCOMITANT MEDICATIONS") |>
    select(USUBJID, CMTRT) |>
    group_by(USUBJID) |>
    summarize(concomitant_medications = n_distinct(CMTRT))
  
  # combine all concomitant medication features
  all_ids <- demographics_raw |> select(USUBJID)
  
  cm <- all_ids |>
    left_join(cm_p, by = "USUBJID") |>
    left_join(cm_n, by = "USUBJID") |>
    # assume that people who were not in cm had zero concomitant medications
    mutate(across(everything(), ~ replace_na(., 0)))
  
  return(cm)
}

#' Extract the modified Hoehn and Yahr stage from questionnaires
#'
#' @param questionnaires_raw A data frame with raw questionnaires for this trial
#' @return A data frame with preprocessed Hoehn and Yahr stage
get_qs_hys_NCT02337738 <- function(questionnaires_raw) {
  
  qs_hys <- questionnaires_raw |>
    select(USUBJID, QSTESTCD, QSSTRESC, QSTPT) |>
    # Modified Hoehn and Yahr stage
    # on-state because UPDRS is measured in on-state
    filter(QSTESTCD == "UPD143" & QSTPT == "ON STATE") |> 
    pivot_wider(names_from = c(QSTESTCD), values_from = QSSTRESC)
  
  return(qs_hys)
}

#' Preprocessing of features from questionnaires that are not MDS-UPDRS scores
#'
#' @param questionnaires_raw A data frame with raw questionnaires for this trial
#' @param updrs_map A data frame that maps UPDRS questions to groups
#' @return A data frame with preprocessed questionnaire features
preprocess_qs_other_NCT02337738 <- function(questionnaires_raw, updrs_map) {
  
  # get modified Hoehn and Yahr stage
  qs_hys <- get_qs_hys_NCT02337738(questionnaires_raw)
  
  # get asymmetry
  qs_asymmetry <- get_qs_asymmetry_Japanese_trial(questionnaires_raw, updrs_map) 
  
  # get TD PIGD phenotypes
  qs_td_pigd <- get_qs_td_pigd_Japanese_trial(questionnaires_raw, updrs_map)
  
  qs_other <- qs_hys |>
    full_join(qs_asymmetry, by = "USUBJID") |>
    full_join(qs_td_pigd, by = "USUBJID")
  
  return(qs_other)
}

#' Add last observation carried forward (LOCF) questionnaire score identifiers 
#'   for all visits of interest
#'
#' @param questionnaires A data frame with questionnaires scores per MDS-UPDRS 
#'   part at each visit
#' @return A data frame with LOCF questionnaires score identifiers per 
#'   MDS-UPDRS part at each visit added
add_qs_locf_identifiers_NCT02337738 <- function(questionnaires) {
  
  qs_locf_id <- questionnaires |>
    group_by(USUBJID) |>
    # for visit 4 at week 6 (402 subjects, 1 subject has only baseline visit, 
    # 1 subject has no data in qs)
    mutate(VISIT_LOCF6 = NA) |>
    # ideal cases (376 subjects)
    mutate(VISIT_LOCF6 = ifelse(
      VISIT == "VISIT 4(Week 6)", 
      "VISIT 4(Week 6)", 
      VISIT_LOCF6
    )) |>
    # discontinued before visit 4 (24 subjects)
    mutate(VISIT_LOCF6 = ifelse(
      VISIT == "VISIT 8(Week 26)/Disc." & !"VISIT 4(Week 6)" %in% VISIT & 
        (sum(!is.na(VISIT)) < 3), 
      "VISIT 4(Week 6)", 
      VISIT_LOCF6
    )) |>
    # no visit 4 but visit 6 (2 subjects)
    mutate(VISIT_LOCF6 = ifelse(
      VISIT == "VISIT 6(Week 14)" & !"VISIT 4(Week 6)" %in% VISIT & 
        (sum(!is.na(VISIT)) == 3), 
      "VISIT 4(Week 6)", 
      VISIT_LOCF6
    )) |>
    
    # for visit 6 at week 14 (402 subjects, 1 subject has only baseline visit, 
    # 1 subject has no data in qs)
    mutate(VISIT_LOCF14 = NA) |>
    # ideal cases (339 subjects)
    mutate(VISIT_LOCF14 = ifelse(
      VISIT == "VISIT 6(Week 14)", 
      "VISIT 6(Week 14)", 
      VISIT_LOCF14
    )) |>
    # discontinued before visit 6 (57 subjects)
    mutate(VISIT_LOCF14 = ifelse(
      (VISIT == "VISIT 8(Week 26)/Disc." & !("VISIT 6(Week 14)" %in% VISIT) & 
         QSDY < 140), 
      "VISIT 6(Week 14)", 
      VISIT_LOCF14
    )) |>
    # two subjects have a discontinued visit, but visit 4 is closer to visit 6 
    # than the discontinued one
    mutate(VISIT_LOCF14 = ifelse(
      (VISIT == "VISIT 4(Week 6)" & !("VISIT 6(Week 14)" %in% VISIT) & 
         ("VISIT 8(Week 26)/Disc." %in% VISIT) & max(QSDY) > 140),
      "VISIT 6(Week 14)", 
      VISIT_LOCF14
    )) |>
    # four subjects have only baseline and visit 4
    mutate(VISIT_LOCF14 = ifelse(
      (VISIT == "VISIT 4(Week 6)" & !("VISIT 6(Week 14)" %in% VISIT) & 
         !("VISIT 8(Week 26)/Disc." %in% VISIT)),
      "VISIT 6(Week 14)", 
      VISIT_LOCF14
    )) |>
    
    # for visit 8 at week 26 (402 subjects, 1 subject has only baseline, 
    # 1 subject has no data in qs)
    mutate(VISIT_LOCF26 = NA) |>
    # most subjects (394) already have this identifier
    mutate(VISIT_LOCF26 = ifelse(
      VISIT == "VISIT 8(Week 26)/Disc.", 
      "VISIT 8(Week 26)/Disc.", 
      VISIT_LOCF26
    )) |>
    # four subjects have only baseline and visit 4
    mutate(VISIT_LOCF26 = ifelse(
      (VISIT == "VISIT 4(Week 6)" & !("VISIT 6(Week 14)" %in% VISIT) & 
         !("VISIT 8(Week 26)/Disc." %in% VISIT)), 
      "VISIT 8(Week 26)/Disc.", 
      VISIT_LOCF26
    )) |>
    # four subjects have only baseline, visit 4 and visit 6
    mutate(VISIT_LOCF26 = ifelse(
      (VISIT == "VISIT 6(Week 14)" & !("VISIT 8(Week 26)/Disc." %in% VISIT)),
      "VISIT 8(Week 26)/Disc.", 
      VISIT_LOCF26
    ))
  return(qs_locf_id)
}

#' Add last observation carried forward (LOCF) questionnaire scores for all 
#'    visits of interest
#'
#' @param questionnaires A data frame with questionnaires scores per MDS-UPDRS 
#'   part at each visit
#' @return A data frame with LOCF questionnaires scores per MDS-UPDRS part at 
#'   each visit added
get_qs_locf_NCT02337738 <- function(questionnaires) {
  
  qs_updrs <- questionnaires
  
  # select needed columns
  qs_updrs <- qs_updrs |>
    select(USUBJID, QSDY, names(qs_updrs)[grepl("VISIT", names(qs_updrs))],
           names(qs_updrs)[grepl("UPDRS", names(qs_updrs))])
  
  # extract baseline
  qs_w0_o <- qs_updrs |>
    filter(VISIT == "VISIT 2(Week 0)") |>
    rename_with(~paste0(., "_w0_o"), c(starts_with("MDS_UPDRS"), QSDY)) |>
    select(USUBJID, QSDY_w0_o, c(starts_with("MDS_UPDRS")))
  
  # extract outcome variations
  qs_w6_locf <- qs_updrs |>
    filter(VISIT_LOCF6 == "VISIT 4(Week 6)") |>
    rename_with(~paste0(., "_w6_locf"), c(starts_with("MDS_UPDRS"), QSDY)) |>
    select(USUBJID, QSDY_w6_locf, c(starts_with("MDS_UPDRS")))
  
  qs_w14_locf <- qs_updrs |>
    filter(VISIT_LOCF14 == "VISIT 6(Week 14)") |>
    rename_with(~paste0(., "_w14_locf"), c(starts_with("MDS_UPDRS"), QSDY)) |>
    select(USUBJID, QSDY_w14_locf, c(starts_with("MDS_UPDRS")))
  
  qs_w26_locf <- qs_updrs |>
    filter(VISIT_LOCF26 == "VISIT 8(Week 26)/Disc.") |>
    rename_with(~paste0(., "_w26_locf"), c(starts_with("MDS_UPDRS"), QSDY)) |>
    select(USUBJID, QSDY_w26_locf, c(starts_with("MDS_UPDRS")))
  
  # combine all the variations
  qs_updrs_locf <- qs_w0_o |>
    left_join(qs_w6_locf, by = "USUBJID") |>
    left_join(qs_w14_locf, by = "USUBJID") |>
    left_join(qs_w26_locf, by = "USUBJID")
  
  # add NA columns for the visits present in other RCTs but missing here
  qs_updrs_locf <- qs_updrs_locf |>
    mutate(MDS_UPDRS23_total_diff_w10_locf = NA) |>
    mutate(MDS_UPDRS23_total_diff_w12_locf = NA) |>
    mutate(MDS_UPDRS23_total_diff_w20_locf = NA) |>
    mutate(MDS_UPDRS23_total_diff_w24_locf = NA) |>
    mutate(MDS_UPDRS23_total_diff_w36_locf = NA)
  
  return(qs_updrs_locf)
}

#' Preprocessing of MDS-UPDRS scores from questionnaires
#'
#' @param questionnaires_raw A data frame with raw questionnaires for this trial
#' @param updrs_map A data frame that maps UPDRS questions to groups
#' @return A data frame with preprocessed questionnaire features
preprocess_qs_updrs_NCT02337738 <- function(questionnaires_raw, updrs_map) {
  
  # visits for MDS-UPDRS
  visits_needed <- c("VISIT 2(Week 0)", "VISIT 4(Week 6)", "VISIT 6(Week 14)", 
                     "VISIT 8(Week 26)/Disc.")
  
  # calculate totals per subject for UPDRS parts and total UPDRS score
  # for the 3 subjects who are missing some questions for part 3, scale to the 
  # approximately expected total
  qs_updrs <- get_qs_updrs_totals_Japanese_trial(
    questionnaires_raw = questionnaires_raw, 
    visits = visits_needed, 
    updrs_map = updrs_map
  )
  
  # calculate differences from baseline for UPDRS parts and total UPDRS score
  qs_updrs <- get_qs_updrs_diff_Japanese_trial(qs_updrs)
  
  # normalize baseline to day 0
  qs_updrs <- qs_updrs |>
    group_by(USUBJID) |>
    mutate(QSDY = QSDY - min(QSDY))
  
  # interpolate UPDRS differences between all known values
  qs_updrs_full <- get_qs_full_int_Japanese_trial(qs_updrs)
  
  # get the interpolated outcome at week 14
  qs_updrs_w14_i <- get_qs_int_14_Japanese_trial(qs_updrs_full)
  
  # add LOCF visit identifiers
  qs_updrs <- add_qs_locf_identifiers_NCT02337738(qs_updrs)
  
  # get MDS-UPDRS differences at LOCF visits
  qs_updrs_locf <- get_qs_locf_NCT02337738(qs_updrs)
  
  # combine questionnaire features
  qs <- full_join(qs_updrs_w14_i, qs_updrs_locf, by = "USUBJID")
  
  return(qs)
}
