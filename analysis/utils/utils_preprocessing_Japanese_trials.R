#' Functions needed for preprocessing of NCT02337725 and NCT02337738, two 
#' Japanese trials in patients with early and advanced Parkinson's disease
#'
#' The functions in this file can be used for both Japanese trials

#' Preprocessing of features from vital signs
#'
#' @param vital_signs_raw A data frame with raw vital signs for this trial
#' @return A data frame with preprocessed vital signs features
preprocess_vs_Japanese_trial <- function(vital_signs_raw) {
  
  vs <- vital_signs_raw |>
    select(USUBJID, VSTESTCD, VSSTRESC, VISIT) |>
    mutate(VSSTRESC = as.numeric(VSSTRESC)) |>
    filter(VSTESTCD %in% c("DIABP", "SYSBP", "WEIGHT")) |>
    filter(VISIT == "VISIT 2(Week 0)") |>
    select(-VISIT) |>
    pivot_wider(names_from = c(VSTESTCD), values_from = VSSTRESC) |>
    # create mean blood pressure feature
    mutate(BP_mean = (DIABP + SYSBP)/2) |>
    select(-c(DIABP, SYSBP))
  
  return(vs)
}

#' Preprocessing of features from Parkinson duration
#'
#' @param parkinson_duration_raw A data frame with raw Parkinson duration
#'   for this trial
#' @return A data frame with preprocessed Parkinson duration features
preprocess_xd_Japanese_trial <- function(parkinson_duration_raw) {
  
  xd <- parkinson_duration_raw |>
    select(USUBJID, XDORRES) |>
    mutate(XDORRES = as.numeric(XDORRES)) |>
    rename(PDDUR = XDORRES) # duration in years
  
  return(xd)
}

#' Preprocessing of features from ECG results
#'
#' @param vital_signs_raw A data frame with raw ECG results for this trial
#' @return A data frame with preprocessed ECG features
preprocess_eg_Japanese_trial <- function(ecg_results_raw) {
  
  # values at screening
  eg_sc <- ecg_results_raw |>
    select(USUBJID, EGTESTCD, EGSTRESC, VISIT) |>
    filter(VISIT=="VISIT 1(Week -2)") |>
    select(-VISIT) |>
    rename(EGSTRESC_sc = EGSTRESC)
  
  # values at baseline
  eg_bl <- ecg_results_raw |>
    select(USUBJID, EGTESTCD, EGSTRESC, VISIT) |>
    filter(VISIT=="VISIT 2(Week 0)") |>
    select(-VISIT) |>
    rename(EGSTRESC_bl = EGSTRESC)
  
  # use the values at screening to complete baseline if missing
  eg <- left_join(eg_bl, eg_sc, by=c("USUBJID", "EGTESTCD")) |>
    mutate(EGSTRESC = ifelse(EGSTRESC_bl == "", EGSTRESC_sc, EGSTRESC_bl)) |>
    select(-c(EGSTRESC_bl, EGSTRESC_sc)) |>
    pivot_wider(names_from = c(EGTESTCD), values_from = EGSTRESC) |>
    mutate(across(-c(USUBJID, INTP), as.numeric)) |>
    mutate(INTP_num = recode_values(
      INTP, 
      "WITHIN NORMAL LIMITS" ~ 0, 
      "ABNORMAL" ~ 1
    )) |>
    mutate(INTP_text = recode_values(
      INTP, 
      "WITHIN NORMAL LIMITS" ~ "Normal", 
      "ABNORMAL" ~ "Abnormal"
    )) |>
    mutate(INTP_fac = factor(INTP_text, levels = c("Normal", "Abnormal")))
  
  return(eg)
}

#' Fill gaps in laboratory results from baseline with values from screening
#'
#' @param laboratory_results_screening A data frame with laboratory results at 
#'   screening
#' @param laboratory_results_baseline A data frame with laboratory results at 
#'   baseline
#' @return A data frame with joined laboratory results
join_lb_Japanese_trial <- function(laboratory_results_screening, 
                                   laboratory_results_baseline) {
  
  # use the values at screening to complete baseline if missing
  lb <- laboratory_results_baseline |>
    full_join(laboratory_results_screening, by=c("USUBJID", "LBTESTCD")) |>
    mutate(LBSTRESC = ifelse(
      LBSTRESC_bl==""|is.na(LBSTRESC_bl), 
      LBSTRESC_sc, 
      LBSTRESC_bl
    )) |>
    select(-c(LBSTRESC_bl, LBSTRESC_sc)) |>
    pivot_wider(names_from = c(LBTESTCD), values_from = LBSTRESC) |>
    mutate(across(-c(USUBJID, GLUC), as.numeric)) |>
    # harmonize units
    mutate(BASO = WBC * BASOLE) |>  # for conversion to 10^9/L
    mutate(EOS = WBC * EOSLE) |>
    mutate(LYM = WBC * LYMLE) |>
    mutate(MONO = WBC * MONOLE) |>
    mutate(NEUT = WBC * NEUTLE) |>
    select(-WBC) |>
    # create numeric glucose
    mutate(GLUC_num = recode_values(
      GLUC, 
      "NEGATIVE" ~ 0, 
      "1+" ~ 1, 
      "2+" ~ 2, 
      "3+" ~ 3, 
      "4+" ~ 4, 
      "5+" ~ 5
    )) |>
    mutate(GLUC_text = recode_values(
      GLUC, 
      "NEGATIVE" ~ "0", 
      "1+" ~ "1+", 
      "2+" ~ "2+", 
      "3+" ~ "3+", 
      "4+" ~ "4+", 
      "5+" ~ "5+"
    )) |>
    mutate(GLUC_fac = factor(
      GLUC_text, 
      levels = c("0", "1+", "2+", "3+", "4+", "5+")
    ))
  return(lb)
}

#' Calculate TD and PIGD subtypes from questionnaires
#'
#' @param questionnaires_raw A data frame with raw questionnaires for this trial
#' @param updrs_map A data frame that maps UPDRS questions to groups
#' @return A data frame with preprocessed TD and PIGD subtypes
get_qs_td_pigd_Japanese_trial <- function(questionnaires_raw, updrs_map) {
  
  td_q <- updrs_map |> 
    filter(version == "MDS-UPDRS") |>
    filter(td_pigd == "TD") |> 
    pull(question)
  
  pigd_q <- updrs_map |> 
    filter(version == "MDS-UPDRS") |>
    filter(td_pigd == "PIGD") |> 
    pull(question)
  
  qs_td_pigd <- questionnaires_raw |>
    select(USUBJID, QSTESTCD, QSSTRESC, VISIT) |>
    mutate(QSSTRESC = as.numeric(QSSTRESC)) |>
    filter(VISIT == "VISIT 2(Week 0)") |>
    select(-VISIT) |>
    pivot_wider(names_from = c(QSTESTCD), values_from = QSSTRESC) |>
    # TD and PIGD phenotypes
    mutate(td_mean = rowSums(across(all_of(td_q)))/length(td_q)) |>
    mutate(pigd_mean = rowSums(across(all_of(pigd_q)))/length(pigd_q)) |>
    mutate(td_pigd = case_when(
      td_mean == 0 & pigd_mean == 0 ~ "Indeterminate",
      td_mean > 0 & pigd_mean == 0 ~ "TD",
      td_mean == 0 & pigd_mean > 0 ~ "PIGD",
      td_mean > 0 & pigd_mean > 0 & td_mean/pigd_mean >= 1.15 ~ "TD",
      td_mean > 0 & pigd_mean > 0 & td_mean/pigd_mean <= 0.9 ~ "PIGD",
      td_mean > 0 & pigd_mean > 0 & td_mean/pigd_mean >= 0.9 & 
        td_mean/pigd_mean <= 1.15 ~ "Indeterminate"
    )) |>
    mutate(td_pigd_num = recode_values(
      td_pigd, 
      "TD" ~ 1, 
      "Indeterminate" ~ 2, 
      "PIGD" ~ 3
    )) |>
    mutate(td_pigd_fac = factor(
      td_pigd, 
      levels = c("TD", "PIGD", "Indeterminate")
    )) |>
    select(USUBJID, td_pigd, td_pigd_num, td_pigd_fac)
  
  return(qs_td_pigd)
}

#' Calculate a standardize motor asymmetry score from questionnaires
#'
#' @param questionnaires_raw A data frame with raw questionnaires for this trial
#' @param updrs_map A data frame that maps UPDRS questions to groups
#' @return A data frame with preprocessed standardize motor asymmetry score
get_qs_asymmetry_Japanese_trial <- function(questionnaires_raw, updrs_map) {
  
  right_q <- updrs_map |> 
    filter(version == "MDS-UPDRS") |>
    filter(asymmetry == "right") |> 
    pull(question)
  
  left_q <- updrs_map |> 
    filter(version == "MDS-UPDRS") |>
    filter(asymmetry == "left") |> 
    pull(question)
  
  qs_asymmetry <- questionnaires_raw |>
    select(USUBJID, QSTESTCD, QSSTRESC, VISIT) |>
    mutate(QSSTRESC = as.numeric(QSSTRESC)) |>
    filter(VISIT == "VISIT 2(Week 0)") |>
    select(-VISIT) |>
    pivot_wider(names_from = c(QSTESTCD), values_from = QSSTRESC) |>
    mutate(right_sum = rowSums(across(all_of(right_q)))) |>
    mutate(left_sum = rowSums(across(all_of(left_q)))) |>
    mutate(asymmetry = ifelse(
      left_sum == 0 & right_sum == 0, 
      0, 
      (left_sum - right_sum)/(left_sum + right_sum))
    ) |>
    select(USUBJID, asymmetry)
  
  return(qs_asymmetry)
}

#' Calculate MDS-UPDRS totals per part at visit times
#'
#' @param questionnaires_raw A data frame with raw questionnaires for this trial
#' @param visits Visits of interest
#' @param updrs_map A data frame that maps UPDRS questions to groups
#' @return A data frame with MDS-UPDRS part totals at visit times
get_qs_updrs_totals_Japanese_trial <- function(questionnaires_raw, visits, updrs_map) {
  
  UPDRS1_q <- updrs_map |> 
    filter(version == "MDS-UPDRS") |>
    filter(part == "part 1") |> 
    pull(question)
  
  UPDRS2_q <- updrs_map |>
    filter(version == "MDS-UPDRS") |>
    filter(part == "part 2") |> 
    pull(question)
  
  UPDRS3_q <- updrs_map |> 
    filter(version == "MDS-UPDRS") |>
    filter(part == "part 3") |> 
    pull(question)
  
  # calculate totals per subject for UPDRS parts and total UPDRS score
  # for the 3 subjects who are missing some questions for part 3, scale to the 
  # approximately expected total
  qs_updrs_totals <- questionnaires_raw |>
    filter(VISIT %in% visits) |>
    select(USUBJID, QSTESTCD, QSSTRESC, VISIT, QSDY) |>
    mutate(QSSTRESC = as.numeric(QSSTRESC)) |>
    # one column per UPDRS question
    pivot_wider(names_from = c(QSTESTCD), values_from = QSSTRESC) |>
    mutate(UPDRS3_na = rowSums(is.na(across(all_of(UPDRS3_q))))) |>
    mutate(
      MDS_UPDRS1_total = rowSums(across(all_of(UPDRS1_q))),
      MDS_UPDRS2_total = rowSums(across(all_of(UPDRS2_q))),
      MDS_UPDRS3_total = ifelse(
        UPDRS3_na > 0 & (length(UPDRS3_q) - UPDRS3_na) > 0, 
        rowSums(across(all_of(UPDRS3_q)), na.rm = TRUE) *
          (length(UPDRS3_q)/(length(UPDRS3_q)-UPDRS3_na)), 
        rowSums(across(all_of(UPDRS3_q)))
      ),
      MDS_UPDRS23_total = rowSums(across(c(MDS_UPDRS2_total, MDS_UPDRS3_total))),
      MDS_UPDRS123_total = rowSums(across(c(MDS_UPDRS1_total, MDS_UPDRS2_total, 
                                            MDS_UPDRS3_total))))
  
  # drop entries that are empty because all questions for part 3 are missing
  qs_updrs_totals <- qs_updrs_totals |>
    filter(!(UPDRS3_na == length(UPDRS3_q)))
  
  return(qs_updrs_totals)
}

#' Calculate UPDRS differences from baseline at each visit
#'
#' @param questionnaires A data frame with questionnaires scores per MDS-UPDRS 
#'   part at each visit
#' @return A data frame with differences from baseline for questionnaires scores 
#'   per UPDRS part at each visit
get_qs_updrs_diff_Japanese_trial <- function(questionnaires) {
  
  qs_updrs_diff <- questionnaires |>
    group_by(USUBJID) |>
    mutate(
      MDS_UPDRS1_total_diff = if (any(VISIT == "VISIT 2(Week 0)")) {
        MDS_UPDRS1_total - MDS_UPDRS1_total[VISIT == "VISIT 2(Week 0)"]
        } else NA,
      MDS_UPDRS2_total_diff = if (any(VISIT == "VISIT 2(Week 0)")) {
        MDS_UPDRS2_total - MDS_UPDRS2_total[VISIT == "VISIT 2(Week 0)"]
        } else NA,
      MDS_UPDRS3_total_diff = if (any(VISIT == "VISIT 2(Week 0)")) {
        MDS_UPDRS3_total - MDS_UPDRS3_total[VISIT == "VISIT 2(Week 0)"] 
        }  else NA,
      MDS_UPDRS23_total_diff = if (any(VISIT == "VISIT 2(Week 0)")) {
        MDS_UPDRS23_total - MDS_UPDRS23_total[VISIT == "VISIT 2(Week 0)"]
        }  else NA,
      MDS_UPDRS123_total_diff = if (any(VISIT == "VISIT 2(Week 0)")) {
        MDS_UPDRS123_total - MDS_UPDRS123_total[VISIT == "VISIT 2(Week 0)"]
        }  else NA
      )
  return(qs_updrs_diff)
}

#' Interpolate scores to all possible days
#'
#' @param questionnaires A data frame outcomes at recorded visit times
#' @return A data frame with outcomes at all possible days
get_qs_full_int_Japanese_trial <- function(questionnaires) {
  
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
    group_by(USUBJID)|>
    # exclude patients who only have baseline
    filter(sum(!is.na(MDS_UPDRS23_total_diff)) > 1) |>
    # interpolate from the known values
    # use last known value for days after the last visit
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
get_qs_int_14_Japanese_trial <- function(questionnaires) {
  
  # get the interpolated outcome at the specified week
  qs_updrs_i <- questionnaires |>
    filter(QSDY == 14*7) |>
    rename("MDS_UPDRS2_total_diff_w14_i" = "MDS_UPDRS2_total_diff_i") |>
    rename("MDS_UPDRS3_total_diff_w14_i" = "MDS_UPDRS3_total_diff_i") |>
    rename("MDS_UPDRS23_total_diff_w14_i" = "MDS_UPDRS23_total_diff_i") |>
    select(c(USUBJID,  MDS_UPDRS2_total_diff_w14_i, MDS_UPDRS3_total_diff_w14_i, 
             MDS_UPDRS23_total_diff_w14_i))
  
  return(qs_updrs_i)
}
