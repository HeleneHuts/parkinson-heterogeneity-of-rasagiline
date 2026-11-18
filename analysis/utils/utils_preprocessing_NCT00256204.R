# Functions needed for preprocessing of NCT00256204 ADAGIO, a trial in 
# patients with early Parkinson's disease
#
# The functions in this file are specifically written for NCT00256204, due to 
# differences in the structure of this trail compared to the other trials

source("analysis/utils/utils_preprocessing_Japanese_trials.R")

#' Get the trial ID of NCT00256204
#'
#' @param dm_all A data frame with raw demographics for all trials
#' @return A string ID that identifies NCT00256204
get_trial_id_NCT00256204 <- function(dm_all) {
  
  trial_arms <- c("ACTIVE ARM 1", "ACTIVE ARM 2", "ACTIVE ARM 3", 
                  "ACTIVE ARM 4", "Not Treated")
  
  dm_id <- dm_all |>
    group_by(STUDYID) |>
    filter(all(trial_arms %in% ACTARM)) |>
    filter(sum(!is.na(ACTARM)) == 1176) |>
    filter(sum(ACTARM == "ACTIVE ARM 1") == 298) |>
    filter(sum(ACTARM == "ACTIVE ARM 2") == 295) |>
    filter(sum(ACTARM == "ACTIVE ARM 3") == 293) |>
    filter(sum(ACTARM == "ACTIVE ARM 4") == 288) |>
    filter(sum(ACTARM == "Not Treated") == 2) |>
    pull(STUDYID)
  
  trial_id <- dm_id[1]
  
  return(trial_id)
}

#' Preprocessing of features from demographics
#'
#' @param demograhics_raw A data frame with raw demographics for this trial
#' @return A data frame with preprocessed demographic features
preprocess_dm_NCT00256204 <- function(demographics_raw) {
  
  dm <- demographics_raw |>
    select(USUBJID, AGE, SEX, ACTARM, COUNTRY) |>
    mutate(arm_d = recode_values(
      ACTARM, 
      "ACTIVE ARM 1" ~ "placebo_1mg_delayed",
      "ACTIVE ARM 2" ~ "placebo_2mg_delayed",
      "ACTIVE ARM 3" ~ "rasagiline_2mg", 
      "ACTIVE ARM 4" ~ "rasagiline_1mg",
      "Not Treated"~ "not_treated"
    )) |>
    mutate(arm = recode_values(
      arm_d, 
      "placebo_1mg_delayed" ~ "placebo",
      "placebo_2mg_delayed" ~ "placebo",
      "rasagiline_2mg" ~ "rasagiline", 
      "rasagiline_1mg" ~ "rasagiline",
      "not_treated"~ "rasagiline"
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
    mutate(region = recode_values(
      COUNTRY, 
      "ARG" ~ "South America", 
      "AUT" ~ "Europe", 
      "CAN" ~ "North America", 
      "DEU" ~ "Europe", 
      "ESP" ~ "Europe", 
      "FRA" ~ "Europe", 
      "GBR" ~ "Europe", 
      "HUN" ~ "Europe", 
      "ISR" ~ "Middle East", 
      "ITA" ~ "Europe", 
      "NLD" ~ "Europe", 
      "PRT" ~ "Europe", 
      "ROU" ~ "Europe", 
      "USA" ~ "North America"
    )) |>
    mutate(region_num = recode_values(
      region, 
      "Asia" ~ 1,
      "Europe" ~ 2,
      "North America" ~ 3,
      "Middle East" ~ 4,
      "South America" ~ 5
    )) |>
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
preprocess_lb_NCT00256204 <- function(laboratory_results_raw) {
  
  lb_features <- c("ALB", "ALT", "AST", "BASO", "BILI", "CA", "CHOL", "CL", 
                   "CREAT", "EOS", "GGT", "GLUC", "HGB", "K", "LDH", "LYM", 
                   "MONO", "NEUT", "PHOS", "PLAT", "PROT", "UREAN", "WBC")
  
  # values at screening
  lb_sc <- laboratory_results_raw |>
    filter(VISIT == "Screening") |>
    select(USUBJID, LBTESTCD, LBSTRESC, LBCAT, LBGRPID, LBDY) |>
    filter(LBTESTCD %in% lb_features) |>
    filter(!(LBTESTCD %in% c("PROT",  "WBC") & LBCAT == "URINALYSIS")) |>
    filter(!(LBTESTCD == "GLUC" & LBCAT == "BIOCHEMISTRY")) |>
    group_by(USUBJID, LBTESTCD) |>
    # two subjects miss a date, set some date so they don't get dropped entirely
    mutate(LBDY_nn = ifelse(is.na(LBDY), -1, LBDY)) |>
    # take the lab value that was recorded closest to baseline
    filter(LBDY_nn == ifelse(
      sum(!is.na(LBSTRESC)) > 1, 
      max(LBDY_nn), 
      LBDY_nn
    )) |> 
    # if multiple values on the same day, take LAB RESULT 2
    filter(!(sum(!is.na(LBSTRESC)) > 1 & LBGRPID == "LAB RESULT 1"))  |>
    ungroup() |>
    select(-c(LBCAT, LBGRPID, LBDY, LBDY_nn))|>
    rename(LBSTRESC_sc = LBSTRESC)
  
  # values at baseline
  lb_bl_selection <- laboratory_results_raw |>
    filter(VISIT == "Baseline") |>
    select(USUBJID, LBTESTCD, LBSTRESC, LBCAT, LBGRPID) |>
    # if more than one value per subject per feature drop one of both
    filter(!(LBTESTCD %in% c("PROT", "RBC", "WBC") & LBCAT == "URINALYSIS")) |>
    filter(!(LBTESTCD == "GLUC" & LBCAT == "BIOCHEMISTRY"))
  
  # get IDs for people who have multiple tests and overlapping features
  lb_bl_overlap <- lb_bl_selection |>
    distinct(USUBJID, LBGRPID, LBTESTCD) |>
    group_by(USUBJID, LBTESTCD) |>
    summarise(LBGRPID_n = n_distinct(LBGRPID), .groups = "drop") |>
    filter(LBGRPID_n > 1) |>
    distinct(USUBJID)
  
  # get number of lab features per test
  lb_bl_overlap_counts <- lb_bl_selection |>
    distinct(USUBJID, LBGRPID, LBTESTCD) |>
    count(USUBJID, LBGRPID, name = "LBTESTCD_n")
  
  # determine which overlapping tests should be kept
  # should be the one with the most recorded results
  lb_bl_tests_to_keep <- lb_bl_overlap_counts |>
    group_by(USUBJID) |>
    filter(!USUBJID %in% lb_bl_overlap$USUBJID | LBTESTCD_n == max(LBTESTCD_n)) |>
    group_by(USUBJID) |>
    # if same number of results for both labs, keep lab result 1
    mutate(keep_flag = "keep") |>
    mutate(keep_flag = ifelse(
      length(unique(LBTESTCD_n)) == 1 & length(LBTESTCD_n) > 1 & 
        LBGRPID == "LAB RESULT 2", 
      "drop", 
      keep_flag
    )) |>
    filter(keep_flag == "keep") |>
    ungroup() |>
    select(USUBJID, LBGRPID)
  
  lb_bl <- lb_bl_selection |>
    semi_join(lb_bl_tests_to_keep, by = c("USUBJID", "LBGRPID")) |>
    select(-c(LBCAT, LBGRPID)) |>
    filter(LBTESTCD %in% lb_features) |>
    ungroup() |>
    rename(LBSTRESC_bl = LBSTRESC)
  
  lb <- lb_bl |>
    full_join(lb_sc, by = c("USUBJID", "LBTESTCD")) |>
    mutate(LBSTRESC = ifelse(
      LBSTRESC_bl==""|is.na(LBSTRESC_bl), 
      LBSTRESC_sc, 
      LBSTRESC_bl
    )) |>
    select(-c(LBSTRESC_bl, LBSTRESC_sc)) |>
    pivot_wider(names_from = c(LBTESTCD), values_from = LBSTRESC) |>
    mutate(across(-c(USUBJID, GLUC), as.numeric)) |>
    # create numeric glucose
    mutate(GLUC_num = recode_values(
      GLUC, 
      "NEGATIVE" ~ 0, 
      "TRACE LEVELS" ~ 0, 
      "0" ~ 0, 
      "100" ~ 1, 
      "+" ~ 1, 
      "250" ~ 2, 
      "++" ~ 2,
      "500" ~ 3, 
      "+++" ~ 3, 
      ">=1000" ~ 4, 
      "++++" ~ 4
    )) |>
    mutate(GLUC_text = recode_values(
      GLUC_num, 
      0 ~ "0", 
      1 ~ "1+", 
      2 ~ "2+",
      3 ~ "3+", 
      4 ~ "4+"
    )) |>
    mutate(GLUC_fac = factor(
      GLUC_text, 
      levels = c("0", "1+", "2+", "3+", "4+", "5+")
    )) |>
    mutate(HGB = HGB * 10) # unit g/L
  return(lb)
}


#' Preprocessing of features from concomitant medication
#'
#' @param concomitant_medication_raw A data frame with raw concomitant 
#'   medication for this trial
#' @param demograhics_raw A data frame with raw demographics for this trial
#'   needed for a full list of USUBJIDs
#' @return A data frame with preprocessed concomitant medication features
preprocess_cm_NCT00256204 <- function(concomitant_medication_raw, 
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
    filter(VISIT == "Baseline" | VISIT == "Screening") |>
    filter(CMENDY >= 0 | is.na(CMENDY)) |>
    select(USUBJID, CMDECOD, CMENDY) |>
    mutate(dopamine_agonists = ifelse(CMDECOD %in% doag & !(CMENDY<(-90)), 1, 0)) |>
    mutate(COMT_inhibitors = ifelse(CMDECOD %in% COMT & !(CMENDY<(-90)), 1, 0)) |>
    mutate(amantadine = ifelse(CMDECOD %in% aman & !(CMENDY<(-90)), 1, 0)) |>
    mutate(anticholinergics = ifelse(CMDECOD %in% anti & !(CMENDY<(-90)), 1, 0)) |>
    mutate(droxidopa = ifelse(CMDECOD %in% drox & !(CMENDY<(-90)), 1, 0)) |>
    mutate(istradefylline = ifelse(CMDECOD %in% istr & !(CMENDY<(-90)), 1, 0)) |>
    mutate(zonisamide = ifelse(CMDECOD %in% zoni & !(CMENDY<(-90)), 1, 0)) |>
    mutate(carbidopa = ifelse(CMDECOD %in% carb & !(CMENDY<(-90)), 1, 0)) |>
    mutate(levodopa = ifelse(CMDECOD %in% levo  & !(CMENDY<(-90)), 1, 0)) |>
    select(-CMENDY) |>
    group_by(USUBJID) |>
    summarize(across(where(is.numeric), max))
  
  # number of concomitant medications per person
  cm_n <- concomitant_medication_raw |>
    filter(VISIT == "Baseline" | VISIT == "Screening") |>
    filter(CMENDY >= 0 | is.na(CMENDY)) |>
    select(USUBJID, CMDECOD) |>
    group_by(USUBJID) |>
    summarize(concomitant_medications=n_distinct(CMDECOD))

  # combine all concomitant medication features
  all_ids <- demographics_raw |> select(USUBJID)
  
  cm <- all_ids |>
    left_join(cm_p, by = "USUBJID") |>
    left_join(cm_n, by = "USUBJID") |>
    # assume that people who were not in cm had zero concomitant medications
    mutate(across(everything(), ~ replace_na(., 0)))
  
  return(cm)
}

#' Preprocessing of features from Parkinson duration
#'
#' @param parkinson_duration_raw A data frame with raw Parkinson duration
#'   for this trial
#' @return A data frame with preprocessed Parkinson duration features
preprocess_xd_NCT00256204 <- function(parkinson_duration_raw) {
  
  xd <- parkinson_duration_raw |>
    select(USUBJID, XDORRES, XDORRESU) |>
    mutate(XDORRES = as.numeric(XDORRES)) |>
    mutate(XDORRES = ifelse(
      XDORRESU == "MONTHS", 
      XDORRES*(52/12*7), 
      XDORRES
    )) |>
    mutate(XDORRES = XDORRES/365) |>  # duration in years
    select(c("USUBJID", "XDORRES")) |>
    rename("PDDUR" = "XDORRES")
  
  return(xd)
}

#' Preprocessing of features from vital signs
#'
#' @param vital_signs_raw A data frame with raw vital signs for this trial
#' @return A data frame with preprocessed vital signs features
preprocess_vs_NCT00256204 <- function(vital_signs_raw) {
  
  # weight at baseline (15 missing)
  vs_w <- vital_signs_raw |>
    filter(VISIT == "Baseline" & VSTESTCD == "WEIGHT") |>
    select(USUBJID, VSTESTCD, VSSTRESC) |>
    pivot_wider(names_from = c(VSTESTCD), values_from = VSSTRESC) |>
    mutate(across(-c(USUBJID), as.numeric))
  
  # subjects who have a baseline flag and a non-standing, non-supine measurement
  # diastolic and systolic blood pressure (21 missing)
  vs_sd <- vital_signs_raw |>
    filter(VISIT == "Baseline" & 
             (VSTESTCD == "DIABP" | VSTESTCD == "SYSBP" |VSTESTCD=="PULSE")) |>
    filter(VSBLFL == "Y" & !(VSPOS == "STANDING" | VSPOS == "SUPINE")) |>
    select(USUBJID, VSTESTCD, VSSTRESC) |>
    pivot_wider(names_from = c(VSTESTCD), values_from = VSSTRESC) |>
    mutate(across(-c(USUBJID), as.numeric))
  
  # subjects who have a baseline flag and only a standing and/or supine measurement
  # calculate a mean of the standing and/or supine measurement
  # diastolic and systolic blood pressure (2 missing)
  vs_sd2 <- vital_signs_raw %>%
    filter(VSTESTCD == "DIABP" | VSTESTCD == "SYSBP" | VSTESTCD=="PULSE") |>
    filter(VISIT == "Baseline" & VSBLFL == "Y") |>
    select(USUBJID, VSTESTCD, VSSTRESC, VSPOS) |>
    group_by(USUBJID, VSTESTCD) |>
    filter(length(unique(VSPOS)) < 3 & !"" %in% VSPOS) |>
    mutate(across(c(VSSTRESC), as.numeric)) |>
    group_modify(~ {
      new_value <- mean(.x$VSSTRESC)
      bind_rows(.x, data.frame(VSPOS="", VSSTRESC=new_value))
    }) |>
    filter(VSPOS == "") |>
    select(USUBJID, VSTESTCD, VSSTRESC) |>
    pivot_wider(names_from = c(VSTESTCD), values_from = VSSTRESC) 
  
  # join and keep all non-NA values
  vs_sd_all <- full_join(vs_sd, vs_sd2, by = "USUBJID") |>
    mutate(PULSE = ifelse(is.na(PULSE.x), PULSE.y, PULSE.x)) |>
    mutate(DIABP = ifelse(is.na(DIABP.x), DIABP.y, DIABP.x)) |>
    mutate(SYSBP = ifelse(is.na(SYSBP.x), SYSBP.y, SYSBP.x)) |>
    select(USUBJID, DIABP, SYSBP) |>
    mutate(BP_mean = (DIABP + SYSBP)/2)
  
  vs <- full_join(vs_w, vs_sd_all, by = "USUBJID")
  
  return(vs)
}

#' Preprocessing of features from ECG results
#'
#' @param vital_signs_raw A data frame with raw ECG results for this trial
#' @return A data frame with preprocessed ECG features
preprocess_eg_NCT00256204 <- function(ecg_results_raw) {
  
  eg <- ecg_results_raw |>
    select(USUBJID, EGTESTCD, EGSTRESC, VISIT, EGDTC, EGREPNUM) |>
    filter(VISIT=="Screening") |>
    mutate(EGDTC = substring(EGDTC, 7, 11)) |>
    mutate(EGDTC = as.time(EGDTC)) |>
    mutate(EGREPNUM = as.numeric(EGREPNUM)) |>
    group_by(USUBJID, EGTESTCD) |>
    filter(EGDTC == min(EGDTC)) |>
    filter(EGREPNUM == min(EGREPNUM)) |>
    select(-c(EGDTC, EGREPNUM, VISIT)) |>
    pivot_wider(names_from = c(EGTESTCD), values_from = EGSTRESC) |>
    select(USUBJID, INTP, PRAG, QRSAG, QTAG, QTCBAG) |>
    mutate(across(c(PRAG, QRSAG, QTAG, QTCBAG), as.numeric)) |>
    mutate(INTP_num = recode_values(
      INTP, 
      "NORMAL" ~ 0, 
      "ABNORMAL" ~ 1
    )) |>
    mutate(INTP_text = recode_values(
      INTP, 
      "NORMAL" ~ "Normal", 
      "ABNORMAL" ~ "Abnormal"
    )) |>
    mutate(INTP_fac = factor(INTP_text, levels = c("Normal", "Abnormal")))
  
  return(eg)
}

#' Extract the modified Hoehn and Yahr stage from questionnaires
#'
#' @param questionnaires_raw A data frame with raw questionnaires for this trial
#' @return A data frame with preprocessed Hoehn and Yahr stage
get_qs_hys_NCT00256204 <- function(questionnaires_raw) {
  
  qs_hys <- questionnaires_raw %>%
    select(USUBJID, QSTESTCD, QSSTRESC, VISIT) |>
    # Modified Hoehn and Yahr stage, measured only at screening
    filter(VISIT == "Screening" & QSTESTCD == "UPD143") |>
    select(-c(VISIT)) |>
    pivot_wider(names_from = c(QSTESTCD), values_from = QSSTRESC)
  
  return(qs_hys)
}

#' Calculate a standardize motor asymmetry score from questionnaires
#'
#' @param questionnaires_raw A data frame with raw questionnaires for this trial
#' @param updrs_map A data frame that maps UPDRS questions to groups
#' @return A data frame with preprocessed standardize motor asymmetry score
get_qs_asymmetry_NCT00256204 <- function(questionnaires_raw, updrs_map) {
  
  right_q <- updrs_map |> 
    filter(version == "UPDRS_ADAGIO") |>
    filter(asymmetry == "right") |> 
    pull(question)
  
  left_q <- updrs_map |> 
    filter(version == "UPDRS_ADAGIO") |>
    filter(asymmetry == "left") |> 
    pull(question)
  
  qs_asymmetry <- questionnaires_raw |>
    select(USUBJID, QSTESTCD, QSSTRESC, VISIT) |>
    mutate(QSSTRESC = as.numeric(QSSTRESC)) |>
    filter(VISIT == "Baseline") |>
    select(-VISIT) |>
    pivot_wider(names_from = c(QSTESTCD), values_from = QSSTRESC) |>
    mutate(right_sum = rowSums(across(all_of(right_q)))) |>
    mutate(left_sum = rowSums(across(all_of(left_q)))) |>
    mutate(asymmetry = ifelse(
      left_sum == 0 & right_sum == 0, 
      0, 
      (left_sum - right_sum)/(left_sum + right_sum))) |>
    select(USUBJID, asymmetry)
  
  return(qs_asymmetry)
}

#' Calculate TD and PIGD subtypes from questionnaires
#'
#' @param questionnaires_raw A data frame with raw questionnaires for this trial
#' @param updrs_map A data frame that maps UPDRS questions to groups
#' @return A data frame with preprocessed TD and PIGD subtypes
get_qs_td_pigd_NCT00256204 <- function(questionnaires_raw, updrs_map) {
  
  td_q <- updrs_map |> 
    filter(version == "UPDRS_ADAGIO") |>
    filter(td_pigd == "TD") |> 
    pull(question)
  
  pigd_q <- updrs_map |> 
    filter(version == "UPDRS_ADAGIO") |>
    filter(td_pigd == "PIGD") |> 
    pull(question)
  
  qs_td_pigd <- questionnaires_raw |>
    select(USUBJID, QSTESTCD, QSSTRESC, VISIT) |>
    mutate(QSSTRESC = as.numeric(QSSTRESC)) |>
    filter(VISIT == "Baseline") |>
    select(-VISIT) |>
    pivot_wider(names_from = c(QSTESTCD), values_from = QSSTRESC) |>
    # TD and PIGD phenotypes
    mutate(td_mean = rowSums(across(all_of(td_q)))/length(td_q)) |>
    mutate(pigd_mean = rowSums(across(all_of(pigd_q)))/length(pigd_q)) |>
    mutate(td_pigd = case_when(
      td_mean == 0 & pigd_mean == 0 ~ "Indeterminate",
      td_mean > 0 & pigd_mean == 0 ~ "TD",
      td_mean == 0 & pigd_mean > 0 ~ "PIGD",
      td_mean > 0 & pigd_mean > 0 & td_mean/pigd_mean >= 1.5 ~ "TD",
      td_mean > 0 & pigd_mean > 0 & td_mean/pigd_mean <= 1 ~ "PIGD",
      td_mean > 0 & pigd_mean > 0 & td_mean/pigd_mean >= 1 & 
        td_mean/pigd_mean <= 1.5 ~ "Indeterminate"
    )) |>
    mutate(td_pigd_fac = factor(
      td_pigd, 
      levels = c("TD", "PIGD", "Indeterminate")
    )) |>
    mutate(td_pigd_num = recode_values(
      td_pigd, 
      "TD" ~ 1, 
      "Indeterminate" ~ 2, 
      "PIGD" ~ 3
    )) |>
    select(USUBJID, td_pigd, td_pigd_num, td_pigd_fac)
  
  return(qs_td_pigd)
}

#' Calculate an approximation of MDS-UPDRS part 1 based on UPDRS part 4 
#'
#' @param questionnaires_raw A data frame with raw questionnaires for this trial
#' @param updrs_map A data frame that maps UPDRS questions to groups
#' @return A data frame with preprocessed MDS-UPDRS part 1 scores
get_qs_mds_updrs1_NCT00256204 <- function(questionnaires_raw, updrs_map) {
  
  # part 4 of the ADAGIO updrs fits with part 1 of the MDS-UPDRS
  # but has two questions less so add some scaling to compensate
  UPDRS4_q <- updrs_map |> 
    filter(version == "UPDRS_ADAGIO") |>
    filter(part == "part 4") |> 
    pull(question)
  
  qs_updrs1 <- questionnaires_raw |>
    select(USUBJID, QSTESTCD, QSSTRESC, VISIT) |>
    mutate(QSSTRESC = as.numeric(QSSTRESC)) |>
    filter(VISIT == "Baseline") |>
    pivot_wider(names_from = c(QSTESTCD), values_from = QSSTRESC) |>
    select(USUBJID, all_of(UPDRS4_q)) |>
    mutate(UPDRS4_total_w0_o = rowSums(across(all_of(UPDRS4_q)))) |>
    mutate(MDS_UPDRS1_total_w0_o = UPDRS4_total_w0_o/44*52) |>
    select(USUBJID, UPDRS4_total_w0_o, MDS_UPDRS1_total_w0_o)
  
  return(qs_updrs1)
}

#' Preprocessing of features from questionnaires that are not MDS-UPDRS scores
#'
#' @param questionnaires_raw A data frame with raw questionnaires for this trial
#' @param updrs_map A data frame that maps UPDRS questions to groups
#' @return A data frame with preprocessed questionnaire features
preprocess_qs_other_NCT00256204 <- function(questionnaires_raw, updrs_map) {
  
  # get modified Hoehn and Yahr stage
  qs_hys <- get_qs_hys_NCT00256204(questionnaires_raw)
  
  # get asymmetry
  qs_asymmetry <- get_qs_asymmetry_NCT00256204(questionnaires_raw, updrs_map)
  
  # get TD PIGD phenotypes
  qs_td_pigd <- get_qs_td_pigd_NCT00256204(questionnaires_raw, updrs_map)
  
  # get an approximated baseline MDS-UPDRS 1
  qs_updrs1 <-  get_qs_mds_updrs1_NCT00256204(questionnaires_raw, updrs_map)
  
  qs_other <- qs_hys |>
    full_join(qs_asymmetry, by = "USUBJID") |>
    full_join(qs_td_pigd, by = "USUBJID") |>
    full_join(qs_updrs1, by = "USUBJID")
}

#' Resolve some issues in the raw questionnaires data
#'
#' @param questionnaires_raw A data frame with raw questionnaires for this trial
#' @param visits Visits of interest
#' @return A data frame with resolved issues
resolve_qs_issues_NCT00256204 <- function(questionnaires_raw, visits) {
  
  qs_updrs <- questionnaires_raw |>
    filter(VISIT %in% visits) |>
    select(USUBJID, QSTESTCD, QSSTRESC, VISIT, QSDY, EPOCH) |>
    mutate(QSSTRESC = as.numeric(QSSTRESC)) |>
    mutate(QSDY = as.numeric(QSDY)) |>
    # one column per UPDRS question
    pivot_wider(names_from = c(QSTESTCD), values_from = QSSTRESC) |>
    # remove issues with dates for visit 4 and part 4, which is not needed here
    filter(!(VISIT == "Visit 4/Week 36" & is.na(UPD105))) |>
    # remove issues with negative data for PFS, which is not needed here
    filter(!(VISIT != "Baseline" & QSDY < 0))
  
  # identify issues with double baseline
  qs_updrs_baseline_issue <- qs_updrs |>
    filter(VISIT == "Baseline") |>
    group_by(USUBJID) |>
    filter(length(VISIT) > 1)
  
  qs_updrs <- qs_updrs |>
    # remove issues with double baseline, keep only the "SCREENING" entry
    filter(!(USUBJID %in% qs_updrs_baseline_issue$USUBJID & 
               VISIT == "Baseline" & EPOCH == "PLACEBO-CONTROLLED PHASE")) |>
    # remove issues with double early termination, keep only relevant UPDRS dates
    filter(!(VISIT == "Early Termination PC Phase" & is.na(UPD105))) |>
    select(-c(names(qs_updrs)[startsWith(names(qs_updrs), "PFS")])) |>
    select(-c(EPOCH, UPD143, UPD148))
  
  return(qs_updrs)
}

#' Calculate UPDRS totals per part at visit times
#'
#' @param questionnaires_raw A data frame with raw questionnaires for this trial
#' @param visits Visits of interest
#' @param updrs_map A data frame that maps UPDRS questions to groups
#' @return A data frame with UPDRS part totals at visit times
get_qs_updrs_totals_NCT00256204 <- function(questionnaires, visits, updrs_map) {
  
  UPDRS2_q <- updrs_map |>
    filter(version == "UPDRS_ADAGIO") |>
    filter(part == "part 2") |> 
    pull(question)
  
  UPDRS3_q <- updrs_map |> 
    filter(version == "UPDRS_ADAGIO") |>
    filter(part == "part 3") |> 
    pull(question)
  
  # calculate score totals
  # for 1 subjects who is missing some questions for part 3, scale to the 
  # approximately expected total
  qs_updrs_totals <- questionnaires |> 
    mutate(UPDRS3_na = rowSums(is.na(across(all_of(UPDRS3_q))))) |>
    mutate(
      UPDRS2_total = rowSums(across(all_of(UPDRS2_q))),
      UPDRS3_total = ifelse(
        UPDRS3_na > 0 & (length(UPDRS3_q) - UPDRS3_na) > 0,
        rowSums(across(all_of(UPDRS3_q)), na.rm = TRUE) *
          (length(UPDRS3_q)/(length(UPDRS3_q)- UPDRS3_na)),
        rowSums(across(all_of(UPDRS3_q)))),
      UPDRS23_total = rowSums(across(c(UPDRS2_total, UPDRS3_total)))
    )
  
  # drop entries that are empty because all questions for part 3 are missing
  qs_updrs_totals <- qs_updrs_totals |>
    filter(!(UPDRS3_na == length(UPDRS3_q)))
  
  return(qs_updrs_totals)
}

#' Normalize the baseline visit to day 0
#'
#' @param questionnaires A data frame outcomes at recorded visit times
#' @return A data frame with normalized visit times
normalize_qs_NCT00256204 <- function(questionnaires) {
  
  qs_updrs_normalized <- questionnaires |>
    # fill QSDY for two baseline visits that do not have a day
    mutate(QSDY = ifelse((is.na(QSDY) & VISIT=="Baseline"), -1, QSDY)) |>
    # fix QSDY for one visit 2 that has day -1 which should be day 83
    mutate(QSDY = ifelse((QSDY == -1 & VISIT == "Visit 2/Week 12"), 83, QSDY)) |>
    group_by(USUBJID) |>
    # normalize baseline to day 0
    mutate(QSDY = QSDY - min(QSDY))
  
  return(qs_updrs_normalized)
}

#' Interpolate scores to all possible days
#'
#' @param questionnaires A data frame outcomes at recorded visit times
#' @return A data frame with outcomes at all possible days
get_qs_full_int_NCT00256204 <- function(questionnaires) {
  
  # outcomes to use for interpolation
  qs_updrs_int <- questionnaires |>
    select(USUBJID, QSDY, MDS_UPDRS2_total_diff, MDS_UPDRS3_total_diff, 
           MDS_UPDRS23_total_diff)
  
  # interpolate between all known values
  # one row for each patient and each interpolated day
  qs_updrs_full <- expand.grid(
    USUBJID = unique(questionnaires$USUBJID), 
    QSDY = c(min(questionnaires$QSDY):max(questionnaires$QSDY))
  ) |>
    merge(qs_updrs_int, by=c("USUBJID", "QSDY"), all.x=TRUE) |>
    group_by(USUBJID) |>
    # exclude patients who only have baseline
    filter(sum(!is.na(MDS_UPDRS23_total_diff)) > 1) |>
    # interpolate from the known values
    mutate(MDS_UPDRS2_total_diff_i = approx(
      x = QSDY[!is.na(MDS_UPDRS2_total_diff)], 
      y = MDS_UPDRS2_total_diff[!is.na(MDS_UPDRS2_total_diff)], 
      xout = QSDY, 
      rule = 2
    )$y
    ) |>
    mutate(MDS_UPDRS3_total_diff_i = approx(
      x = QSDY[!is.na(MDS_UPDRS3_total_diff)], 
      y = MDS_UPDRS3_total_diff[!is.na(MDS_UPDRS3_total_diff)], 
      xout = QSDY, 
      rule = 2
    )$y
    ) |>
    mutate(MDS_UPDRS23_total_diff_i = approx(
      x = QSDY[!is.na(MDS_UPDRS23_total_diff)],
      y = MDS_UPDRS23_total_diff[!is.na(MDS_UPDRS23_total_diff)], 
      xout = QSDY, 
      rule = 2
    )$y
    )
  return(qs_updrs_full)
}

#' Get interpolated outcomes at week 14
#'
#' @param questionnaires A data frame with outcomes at all possible days
#' @return A data frame with outcomes at week 14
get_qs_int_14_NCT00256204 <- function(questionnaires) {
  
  # get the interpolated outcome at week 12, then scale to week 14
  qs_updrs_i <- questionnaires |>
    # extract the interpolated value at week 12 (actual measurement date in adagio)
    filter(QSDY == 84) |>
    # scale to week 14
    mutate(
      MDS_UPDRS2_total_diff_w14_i = MDS_UPDRS2_total_diff_i * (14/12),
      MDS_UPDRS3_total_diff_w14_i = MDS_UPDRS3_total_diff_i * (14/12),
      MDS_UPDRS23_total_diff_w14_i = MDS_UPDRS23_total_diff_i * (14/12)
    ) |>
    select(c(USUBJID, QSDY, MDS_UPDRS2_total_diff_w14_i, 
             MDS_UPDRS3_total_diff_w14_i, MDS_UPDRS23_total_diff_w14_i))
  
  return(qs_updrs_i)
}

#' Add last observation carried forward (LOCF) questionnaire score identifiers 
#'   for all visits of interest
#'
#' @param questionnaires A data frame with questionnaires scores per MDS-UPDRS 
#'   part at each visit
#' @return A data frame with LOCF questionnaires score identifiers per 
#'   MDS-UPDRS part at each visit added
add_qs_locf_NCT00256204 <- function(questionnaires) {
  
  # extract one needed row for one subject who only misses visit 2
  qs_updrs_missing2 <- questionnaires |>  
    group_by(USUBJID) |>
    filter(!"Visit 2/Week 12" %in% VISIT) |>
    filter("Visit 3/Week 24" %in% VISIT) |>
    filter("Visit 4/Week 36" %in% VISIT) |>
    filter(VISIT == "Baseline") |>  # use baseline as LOCF
    mutate(VISIT_LOCF12 = "Visit 2/Week 12") |>
    relocate(VISIT_LOCF12, .after = VISIT) |>
    mutate(VISIT = "Visit 2/Week 12")
  
  qs_locf_id <- questionnaires |>
    # LOCF for visit 2 (1164 subjects total, 12 subjects only have baseline)
    mutate(VISIT_LOCF12 = NA) |>
    relocate(VISIT_LOCF12, .after = VISIT) |>
    group_by(USUBJID) |>
    # merge 1 subject who was missing only visit 2
    bind_rows(qs_updrs_missing2) |>
    # ideal cases (1108 subjects)
    mutate(VISIT_LOCF12 = ifelse(
      VISIT == "Visit 2/Week 12", 
      "Visit 2/Week 12", 
      VISIT_LOCF12
    )) |>
    # early termination or early transfer (56 subjects)
    mutate(VISIT_LOCF12 = ifelse(
      (VISIT == "Visit 4/Week 36" | VISIT == "Early Termination PC Phase") & 
        !"Visit 2/Week 12" %in% VISIT & QSDY < 168, 
      "Visit 2/Week 12", 
      VISIT_LOCF12
    )) |>
    
    # LOCF for visit 3 (1164 subjects total, 12 subjects only have baseline)
    mutate(VISIT_LOCF24 = NA) |>
    # ideal cases (1007 subjects)
    mutate(VISIT_LOCF24 = ifelse(
      VISIT=="Visit 3/Week 24", 
      "Visit 3/Week 24", 
      VISIT_LOCF24
    )) |>
    # early termination (57 subjects)
    mutate(VISIT_LOCF24 = ifelse(
      VISIT=="Early Termination PC Phase" & !"Visit 3/Week 24" %in% VISIT, 
      "Visit 3/Week 24", 
      VISIT_LOCF24
    )) |>
    # early transfer (97 subjects)
    mutate(VISIT_LOCF24 = ifelse(
      VISIT=="Visit 4/Week 36" & !"Visit 3/Week 24" %in% VISIT, 
      "Visit 3/Week 24", 
      VISIT_LOCF24
    )) |>
    # 3 subjects have only visit 2
    mutate(VISIT_LOCF24 = ifelse(
      VISIT=="Visit 2/Week 12" & !"Visit 3/Week 24" %in% VISIT & 
        !"Early Termination PC Phase" %in% VISIT & !"Visit 4/Week 36" %in% VISIT, 
      "Visit 3/Week 24", 
      VISIT_LOCF24
    )) |>
    
    # LOCF for visit 4 (for subjects who get past the 2nd visit)
    mutate(VISIT_LOCF36_long = NA) |>
    relocate(VISIT_LOCF36_long, .after = VISIT_LOCF12) |>
    # ideal cases (1070 subjects)
    mutate(VISIT_LOCF36_long = ifelse(
      VISIT == "Visit 4/Week 36" & "Visit 2/Week 12" %in% VISIT, 
      "Visit 4/Week 36", 
      VISIT_LOCF36_long
    )) |>
    # 32 subjects early termination after visit 2
    mutate(VISIT_LOCF36_long = ifelse(
      VISIT == "Early Termination PC Phase" & "Visit 2/Week 12" %in% VISIT & 
        !"Visit 4/Week 36" %in% VISIT, 
      "Visit 4/Week 36", 
      VISIT_LOCF36_long
    )) |>
    # 3 subjects visit 3 after visit 2, but no visit 4 or early termination
    mutate(VISIT_LOCF36_long = ifelse(
      VISIT == "Visit 3/Week 24" & "Visit 2/Week 12" %in% VISIT & 
        !"Visit 4/Week 36" %in% VISIT & !"Early Termination PC Phase" %in% VISIT,
      "Visit 4/Week 36", 
      VISIT_LOCF36_long
    )) |>
    
    # LOCF for visit 4 (1164 subjects total, 12 subjects only have baseline)
    mutate(VISIT_LOCF36 = VISIT_LOCF36_long) |>
    relocate(VISIT_LOCF36, .after = VISIT_LOCF12) |>
    # 56 early termination before visit 2
    mutate(VISIT_LOCF36 = ifelse(
      (VISIT == "Visit 4/Week 36" | VISIT == "Early Termination PC Phase") & 
        !"Visit 2/Week 12" %in% VISIT & QSDY < 168, 
      "Visit 4/Week 36", 
      VISIT_LOCF36
    )) |>
    # nothing after visit 2 (3 subjects)
    mutate(VISIT_LOCF36 = ifelse(
      VISIT == "Visit 2/Week 12" & !"Visit 3/Week 24" %in% VISIT & 
        !"Visit 4/Week 36" %in% VISIT & !"Early Termination PC Phase" %in% VISIT, 
      "Visit 4/Week 36", 
      VISIT_LOCF36
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
get_qs_locf_NCT00256204 <- function(questionnaires) {
  
  qs_updrs <- questionnaires
  
  # select needed columns
  qs_updrs <- qs_updrs |>
    select(USUBJID, QSDY, names(qs_updrs)[grepl("VISIT", names(qs_updrs))], 
           names(qs_updrs)[grepl("UPDRS", names(qs_updrs))])
  
  # extract baseline variations
  qs_updrs_w0_o <- qs_updrs |>
    filter(VISIT == "Baseline") |>
    rename_with(
      ~paste0(., "_w0_o"), 
      c(starts_with("UPDRS"), starts_with("MDS_UPDRS"), QSDY)
    ) |>
    select(USUBJID, QSDY_w0_o, c(starts_with("UPDRS")), 
           c(starts_with("MDS_UPDRS")))
  
  # extract outcome variations
  qs_updrs_w12_locf <- qs_updrs |>
    filter(VISIT_LOCF12 == "Visit 2/Week 12") |>
    rename_with(
      ~paste0(., "_w12_locf"), 
      c(starts_with("UPDRS"), starts_with("MDS_UPDRS"), QSDY)
    ) |>
    select(USUBJID, QSDY_w12_locf, c(starts_with("UPDRS")), 
           c(starts_with("MDS_UPDRS")))
  
  qs_updrs_w24_locf <- qs_updrs |>
    filter(VISIT_LOCF24 == "Visit 3/Week 24") |>
    rename_with(
      ~paste0(., "_w24_locf"), 
      c(starts_with("UPDRS"), starts_with("MDS_UPDRS"), QSDY)
    ) |>
    select(USUBJID, QSDY_w24_locf, c(starts_with("UPDRS")), 
           c(starts_with("MDS_UPDRS")))
  
  qs_updrs_w36_locf <- qs_updrs |>
    filter(VISIT_LOCF36 == "Visit 4/Week 36") |>
    rename_with(
      ~paste0(., "_w36_locf"), 
      c(starts_with("UPDRS"), starts_with("MDS_UPDRS"), QSDY)
    ) |>
    select(USUBJID, QSDY_w36_locf, c(starts_with("UPDRS")), 
           c(starts_with("MDS_UPDRS")))
  
  # combine all the variations
  qs_updrs_locf <- qs_updrs_w0_o |>
    left_join(qs_updrs_w12_locf, by = "USUBJID") |>
    left_join(qs_updrs_w24_locf, by = "USUBJID") |>
    left_join(qs_updrs_w36_locf, by = "USUBJID")
  
  # add NA columns for the visits present in other RCTs but missing here
  qs_updrs_locf <- qs_updrs_locf |>
    mutate(MDS_UPDRS23_total_diff_w6_locf = NA) |>
    mutate(MDS_UPDRS23_total_diff_w10_locf = NA) |>
    mutate(MDS_UPDRS23_total_diff_w14_locf = NA) |>
    mutate(MDS_UPDRS23_total_diff_w20_locf = NA) |>
    mutate(MDS_UPDRS23_total_diff_w26_locf = NA)
  
  return(qs_updrs_locf)
}

#' Convert ADAGIO UPDRS parts to MDS-UPDRS parts
#'
#' @param questionnaires_raw A data frame with raw questionnaires for this trial
#'   this is needed to extract Hoehn and Yahr stages
#' @param questionnaires A data frame with questionnaires scores per UPDRS part 
#'   at each visit
#' @return A data frame with questionnaires scores per UPDRS and MDS-UPDRS part 
#'   at each visit
convert_qs_NCT00256204 <- function(questionnaires_raw, questionnaires) {
  
  # get modified Hoehn and Yahr stage
  qs_hys <- get_qs_hys_NCT00256204(questionnaires_raw)
  
  # convert UPDRS totals to MDS-UPDRS totals
  qs_updrs_convert <- questionnaires |> 
    # get Hoehn and Yahr from screening
    left_join(qs_hys, by="USUBJID") |>
    mutate(
      MDS_UPDRS2_total = ifelse(
        UPD143 < 3, 
        (UPDRS2_total*1.1) + 0.2, 
        (UPDRS2_total*1) + 1.5),
      MDS_UPDRS3_total = ifelse(
        UPD143 < 3, 
        (UPDRS3_total*1.2) + 2.3, 
        (UPDRS3_total*1.2) + 1),
      MDS_UPDRS23_total = rowSums(across(c(MDS_UPDRS2_total, MDS_UPDRS3_total)))
    )
  return(qs_updrs_convert)
}

#' Calculate UPDRS differences from baseline at each visit
#'
#' @param questionnaires A data frame with questionnaires scores per UPDRS and 
#'   MDS-UPDRS part at each visit
#' @return A data frame with differences from baseline for questionnaires scores 
#'   per UPDRS part at each visit
get_qs_updrs_diff_NCT00256204 <- function(questionnaires) {
  
  # calculate score difference from baseline
  qs_updrs_diff <- questionnaires |>
    group_by(USUBJID) |>
    mutate(
      UPDRS2_total_diff = if (any(VISIT == "Baseline")) {
        UPDRS2_total - UPDRS2_total[VISIT == "Baseline"]
        } else NA,
      UPDRS3_total_diff = if (any(VISIT == "Baseline")) {
        UPDRS3_total - UPDRS3_total[VISIT == "Baseline"]
        } else NA,
      UPDRS23_total_diff = if (any(VISIT == "Baseline")) {
        UPDRS23_total - UPDRS23_total[VISIT == "Baseline"]
        } else NA,
      MDS_UPDRS2_total_diff = if (any(VISIT == "Baseline")) {
        MDS_UPDRS2_total - MDS_UPDRS2_total[VISIT == "Baseline"]
        } else NA,
      MDS_UPDRS3_total_diff = if (any(VISIT == "Baseline")) {
        MDS_UPDRS3_total - MDS_UPDRS3_total[VISIT == "Baseline"]
        } else NA,
      MDS_UPDRS23_total_diff = if (any(VISIT == "Baseline")) {
        MDS_UPDRS23_total - MDS_UPDRS23_total[VISIT == "Baseline"]
        } else NA
    ) |>
    ungroup()
  return(qs_updrs_diff)
}

#' Preprocessing of MDS-UPDRS scores from questionnaires
#'
#' @param questionnaires_raw A data frame with raw questionnaires for this trial
#' @param updrs_map A data frame that maps UPDRS questions to groups
#' @return A data frame with preprocessed questionnaire features
preprocess_qs_updrs_NCT00256204 <- function(questionnaires_raw, updrs_map) {
  
  # visits for placebo-controlled phase
  visits_pc <- c("Baseline", "Visit 2/Week 12", "Visit 3/Week 24", 
                 "Visit 4/Week 36", "Early Termination PC Phase")
  
  # fix some issues in the raw data
  qs_updrs <- resolve_qs_issues_NCT00256204(questionnaires_raw, visits_pc)
  
  # calculate totals per subject for UPDRS parts and total UPDRS score
  qs_updrs <- get_qs_updrs_totals_NCT00256204(qs_updrs, visits_pc, updrs_map)
  
  # convert to MDS-UPDRS
  qs_updrs <- convert_qs_NCT00256204(questionnaires_raw, qs_updrs)
  
  # calculate differences from baseline for UPDRS parts and total UPDRS score
  qs_updrs <- get_qs_updrs_diff_NCT00256204(qs_updrs)
  
  # normalize baseline to day 0
  qs_updrs <- normalize_qs_NCT00256204(qs_updrs)
  
  # interpolate UPDRS differences between all known values
  qs_updrs_full <- get_qs_full_int_NCT00256204(qs_updrs)
  
  # get the interpolated outcome at week 14
  qs_updrs_w14_i <- get_qs_int_14_NCT00256204(qs_updrs_full)
  
  # add LOCF visit identifiers
  qs_updrs <- add_qs_locf_NCT00256204(qs_updrs)
  
  # get MDS-UPDRS differences at LOCF visits
  qs_updrs_locf <- get_qs_locf_NCT00256204(qs_updrs)
  
  # combine questionnaire features
  qs <- full_join(qs_updrs_w14_i, qs_updrs_locf, by = "USUBJID")
  
  return(qs)
}

