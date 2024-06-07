library(CVrisk)

#' Framingham 2008 ASCVD risk score (no lab measurement)
#'
#' Computes 10-year risk for ASCVD event (coronary death, myocardial
#' infarction (MI),coronary insufficiency, angina, ischemic stroke,
#' hemorrhagic stroke, transient ischemic attack, peripheral artery
#' disease, or heart failure).
#'
#' @param gender patient gender (male, female)
#' @param age patient age (years), between 30 and 74
#' @param bmi Body mass index (kg/m2)
#' @param sbp Systolic blood pressure (mm Hg)
#' @param bp_med Patient is on a blood pressure medication (1=Yes, 0=No)
#' @param smoker Current smoker (1=Yes, 0=No)
#' @param diabetes Diabetes (1=Yes, 0=No)
#' @param ... Additional predictors can be passed and will be ignored
#'
#' @return Estimated 10-Y Risk for hard ASCVD (percent)
#'
#' @export
#'
#' @examples
#' library(CVrisk)
#' ascvd_10y_frs_simple(
#'   gender = "male", age = 55,
#'   bmi = 30, sbp = 140,
#'   bp_med = 0, smoker = 0, diabetes = 0
#' )
#'
#' # 16.7
#' @references
#' D’agostino, R.B., Vasan, R.S., Pencina, M.J., Wolf, P.A., Cobain, M.,
#' Massaro, J.M. and Kannel, W.B., 2008. General cardiovascular risk
#' profile for use in primary care: the Framingham Heart Study.
#' Circulation, 117(6), pp.743-753.

adapted_ascvd_10y_frs_simple <- function(gender = c("male", "female"),
                                 age, bmi, sbp,
                                 bp_med, smoker, diabetes, ...) {
  gender <- tolower(gender)
  gender <- ifelse(gender == "m", "male", gender)
  gender <- ifelse(gender == "f", "female", gender)
  
  if (!all(gender %in% c("male", "female")) | missing(gender)) {
    stop("gender must be either 'male' or 'female'")
  }
  
  if (!is.numeric(age) | missing(age)) {
    stop("age must be a valid numeric value'")
  }
  
  # age <- ifelse(age < 30 | age > 74, NA, age)
  
  if (missing(bmi)) {
    bmi <- NA
  }
  
  if (missing(sbp)) {
    sbp <- NA
  }
  
  if (missing(bp_med)) {
    bp_med <- NA
  }
  
  if (missing(smoker)) {
    smoker <- NA
  }
  
  if (missing(diabetes)) {
    diabetes <- NA
  }
  
  # retrieve model coefficients
  # frs_simple_coef <- NULL
  # utils::data(frs_simple_coef, envir = environment())
  frs_coef <- CVrisk::frs_simple_coef
  
  # Generate data.frame of coefficients based on input `gender` vector.
  # We lose the original order after the merge operation, so will need
  # to re-sort the output based on the original order of `sex_df`.
  sex_df <- data.frame(gender)
  sex_df$id <- as.numeric(row.names(sex_df))
  
  model_coef <- merge(sex_df, frs_simple_coef)
  model_coef <- model_coef[order(model_coef$id), ]
  
  sbp_treated <- ifelse(bp_med == 1, sbp, 1)
  sbp_untreated <- ifelse(bp_med == 0, sbp, 1)
  
  indv_sum <- log(age) * model_coef$ln_age +
    log(bmi) * model_coef$ln_bmi +
    log(sbp_treated) * model_coef$ln_treated_sbp +
    log(sbp_untreated) * model_coef$ln_untreated_sbp +
    smoker * model_coef$smoker +
    diabetes * model_coef$diabetes
  
  risk_score <- (1 - (model_coef$baseline_survival^
                              exp(indv_sum - model_coef$group_mean))) * 100.0
  
  # ifelse(risk_score < 1, 1, ifelse(risk_score > 30, 30, risk_score))
}

#' Framingham 2008 ASCVD risk score (with lab measurement)
#'
#' Computes 10-year risk for ASCVD event (coronary death, myocardial
#' infarction (MI), coronary insufficiency, angina, ischemic stroke,
#' hemorrhagic stroke, transient ischemic attack, peripheral artery disease,
#' or heart failure).
#'
#' @param gender patient gender (male, female)
#' @param age patient age (years), between 30 and 74
#' @param hdl HDL cholesterol (mg/dL)
#' @param totchol Total cholesterol (mg/dL)
#' @param sbp Systolic blood pressure (mm Hg)
#' @param bp_med Patient is on a blood pressure medication (1=Yes, 0=No)
#' @param smoker Current smoker (1=Yes, 0=No)
#' @param diabetes Diabetes (1=Yes, 0=No)
#' @param ... Additional predictors can be passed and will be ignored
#'
#' @return Estimated 10-Y Risk for hard ASCVD event (percent)
#'
#' @export
#'
#' @examples
#' library(CVrisk)
#' ascvd_10y_frs(
#'   gender = "male", age = 55,
#'   hdl = 50, totchol = 213, sbp = 140,
#'   bp_med = 0, smoker = 0, diabetes = 0
#' )
#'
#' # 16.7
#' @references
#' D’agostino, R.B., Vasan, R.S., Pencina, M.J., Wolf, P.A., Cobain, M.,
#' Massaro, J.M. and Kannel, W.B., 2008. General cardiovascular risk
#' profile for use in primary care: the Framingham Heart Study.
#' Circulation, 117(6), pp.743-753.

adapted_ascvd_10y_frs <- function(gender = c("male", "female"),
                          age, hdl, totchol, sbp,
                          bp_med, smoker, diabetes, ...) {
  gender <- tolower(gender)
  gender <- ifelse(gender == "m", "male", gender)
  gender <- ifelse(gender == "f", "female", gender)
  
  if (!all(gender %in% c("male", "female")) | missing(gender)) {
    stop("gender must be either 'male' or 'female'")
  }
  
  if (!is.numeric(age) | missing(age)) {
    stop("age must be a valid numeric value'")
  }
  
  # age <- ifelse(age < 30 | age > 74, NA, age)
  
  if (missing(hdl)) {
    hdl <- NA
  }
  
  if (missing(totchol)) {
    totchol <- NA
  }
  
  if (missing(sbp)) {
    sbp <- NA
  }
  
  if (missing(bp_med)) {
    bp_med <- NA
  }
  
  if (missing(smoker)) {
    smoker <- NA
  }
  
  if (missing(diabetes)) {
    diabetes <- NA
  }
  
  # retrieve model coefficients
  # frs_coef <- NULL
  # utils::data(frs_coef, envir = environment())
  frs_coef <- CVrisk::frs_coef
  
  # Generate data.frame of coefficients based on input `gender` vector.
  # We lose the original order after the merge operation, so will need
  # to re-sort the output based on the original order of `sex_df`.
  sex_df <- data.frame(gender)
  sex_df$id <- as.numeric(row.names(sex_df))
  
  model_coef <- merge(sex_df, frs_coef)
  model_coef <- model_coef[order(model_coef$id), ]
  
  sbp_treated <- ifelse(bp_med == 1, sbp, 1)
  sbp_untreated <- ifelse(bp_med == 0, sbp, 1)
  
  indv_sum <- log(age) * model_coef$ln_age +
    log(hdl) * model_coef$ln_hdl +
    log(totchol) * model_coef$ln_totchol +
    log(sbp_treated) * model_coef$ln_treated_sbp +
    log(sbp_untreated) * model_coef$ln_untreated_sbp +
    smoker * model_coef$smoker +
    diabetes * model_coef$diabetes
  
  risk_score <- (1 - (model_coef$baseline_survival^
                              exp(indv_sum - model_coef$group_mean))) * 100.0
  
  # ifelse(risk_score < 1, 1, ifelse(risk_score > 30, 30, risk_score))
  
  # ifelse(risk_score < 1, 1, risk_score)
}

#' ACC/AHA 2013 ASCVD risk score
#'
#' Computes 10-year risk for hard ASCVD event (defined as first occurrence of
#' non-fatal myocardial infarction (MI), congestive heart disease (CHD) death,
#' or fatal or nonfatal stroke).
#'
#' @param race patient race (white, aa, other)
#' @param gender patient gender (male, female)
#' @param age patient age (years)
#' @param totchol Total cholesterol (mg/dL)
#' @param hdl HDL cholesterol (mg/dL)
#' @param sbp Systolic blood pressure (mm Hg)
#' @param bp_med Patient is on a blood pressure medication (1=Yes, 0=No)
#' @param smoker Current smoker (1=Yes, 0=No)
#' @param diabetes Diabetes (1=Yes, 0=No)
#' @param ... Additional predictors can be passed and will be ignored
#'
#'
#' @return Estimated 10-Y Risk for hard ASCVD (percent)
#'
#' @export
#'
#' @examples
#' library(CVrisk)
#' ascvd_10y_accaha(
#'   race = "aa", gender = "male", age = 55,
#'   totchol = 213, hdl = 50, sbp = 140,
#'   bp_med = 0, smoker = 0, diabetes = 0
#' )
#' @references
#' Goff, David C., et al. "2013 ACC/AHA guideline on the assessment of
#' cardiovascular risk: a report of the American College of
#' Cardiology/American Heart Association Task Force on Practice
#' Guidelines." Journal of the American College of Cardiology 63.25
#' Part B (2014): 2935-2959.

adapted_ascvd_10y_accaha <- function(race = "white", gender = c("male", "female"),
                             age, totchol, hdl, sbp,
                             bp_med, smoker, diabetes, ...) {
  
  if (!all(gender %in% c("male", "female")) | missing(gender)) {
    stop("gender must be either 'male' or 'female'")
  }
  
  # age <- ifelse(!is.numeric(age) |
  #                 age < 20 | age > 79 | missing(age), NA, age)
  # 
  # totchol <- ifelse(!is.numeric(totchol) |
  #                     totchol < 130 | totchol > 320 | missing(totchol), NA, totchol)
  # 
  # hdl <- ifelse(!is.numeric(hdl) |
  #                 hdl < 20 | hdl > 100 | missing(hdl), NA, hdl)
  # 
  # sbp <- ifelse(!is.numeric(sbp) |
  #                 sbp < 90 | sbp > 200 | missing(sbp), NA, sbp)
  
  # ascvd_pooled_coef <- NULL
  # utils::data(ascvd_pooled_coef, envir = environment())
  ascvd_pooled_coef <- CVrisk::ascvd_pooled_coef
  
  # Generate data.frame of coefficients based on input `race` and `gender`
  # vectors. We lose the original order after the merge operation, so will
  # need to re-sort the output based on the original order of `race_sex`.
  
  race <- ifelse(race %in% c("white", "aa"), race, "white")
  
  race_sex <- data.frame(race, gender)
  race_sex$id <- as.numeric(row.names(race_sex))
  
  pooled_coef <- merge(race_sex, ascvd_pooled_coef)
  pooled_coef <- pooled_coef[order(pooled_coef$id), ]
  
  sbp_treated <- ifelse(bp_med == 1, sbp, 1)
  sbp_untreated <- ifelse(bp_med == 0, sbp, 1)
  
  indv_sum <- log(age) * pooled_coef$ln_age +
    log(age)^2 * pooled_coef$ln_age_squared +
    log(totchol) * pooled_coef$ln_totchol +
    log(age) * log(totchol) * pooled_coef$ln_age_totchol +
    log(hdl) * pooled_coef$ln_hdl +
    log(age) * log(hdl) * pooled_coef$ln_age_hdl +
    log(sbp_treated) * pooled_coef$ln_treated_sbp +
    log(sbp_treated) * log(age) * pooled_coef$ln_age_treated_sbp +
    log(sbp_untreated) * pooled_coef$ln_untreated_sbp +
    log(sbp_untreated) * log(age) * pooled_coef$ln_age_untreated_sbp +
    smoker * pooled_coef$smoker +
    smoker * log(age) * pooled_coef$ln_age_smoker +
    diabetes * pooled_coef$diabetes
  
  risk_score <- (1 - (pooled_coef$baseline_survival^
                              exp(indv_sum - pooled_coef$group_mean))) * 100.0
  
  # ifelse(risk_score < 1, 1, risk_score)
  
}

#' Compute multiple CV risk scores
#'
#' @param df input dataframe
#' @param scores scores to compute, default is all scores
#' @param age patient age in years (required for all scores)
#' @param gender patient gender (male or female)
#' @param race character string for patient race (white, aa, other) column
#' @param sbp character string of systolic blood pressure (in mm Hg) column
#' @param bmi character string of Body mass index (kg/m2) column
#' @param hdl character string of HDL column
#' @param totchol character string of total cholesterol column
#' @param bp_med character string of blood pressure medication column
#' @param smoker character string of smoking status column
#' @param diabetes character string of diabetes status column
#' @param lipid_med character string of lipid medication column
#' @param fh_heartattack character string of fh of heart attack status column
#' @param cac character string of cac column
#'
#' @return input data frame with risk score results appended as columns
#'
#' @examples
#'
#' library(CVrisk)
#' compute_CVrisk(sample_data,
#'   age = "age", race = "race", gender = "gender", bmi = "BMI", sbp = "sbp",
#'   hdl = "hdl", totchol = "totchol", bp_med = "bp_med", smoker = "smoker",
#'   diabetes = "diabetes", lipid_med = "lipid_med",
#'   fh_heartattack = "fh_heartattack", cac = "cac"
#' )
#' @export
adapted_compute_CVrisk <- function(df, scores = c(
  "adapted_ascvd_10y_accaha",
  "adapted_ascvd_10y_frs", 
  "adapted_ascvd_10y_frs_simple"
),
age, gender, race, sbp = NULL, bmi = NULL,
hdl = NULL, totchol = NULL,
bp_med = NULL, smoker = NULL, diabetes = NULL,
lipid_med = NULL, fh_heartattack = NULL,
cac = NULL) {
  all_args <- as.list(environment())
  valid_pred <- c(
    "age", "gender", "race", "sbp", "bmi", "hdl", "totchol",
    "bp_med", "smoker", "diabetes", "lipid_med", "fh_heartattack", "cac"
  )
  
  pred_args <- list()
  for (var in valid_pred) {
    if (!is.null(all_args[[var]])) pred_args[[var]] <- df[[all_args[[var]]]]
  }
  
  
  results <- sapply(scores, function(x) do.call(x, pred_args))
  
  row.names(results) <- NULL
  return(cbind(df, results))
}