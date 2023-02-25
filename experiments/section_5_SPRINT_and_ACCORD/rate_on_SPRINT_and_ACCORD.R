#!/usr/bin/env Rscript

# Install grf if not already installed
# devtools::install_github("erikcs/grf", subdir = "r-package/grf", ref = "RATE")
library(tidyverse, warn.conflicts=FALSE)

# If calling from the command line, assign train.str and test.str to user specs
args = commandArgs(trailingOnly=TRUE)
train.str <- tolower(args[1])
test.str <- tolower(args[2])

seed <- 42
set.seed(seed)  # Required to get consistent results with the RATE

SPRINT.df <- read_csv("../../data/sprint/sprint_cut.csv", show_col_types=FALSE)
ACCORD.df <- read_csv("../../data/accord/accord_cut.csv", show_col_types=FALSE)

# Add in a diabetes flag
SPRINT.df$Diabetes <- 0
ACCORD.df$Diabetes <- 1

# Set the train/test data according to command-line argument
if (train.str == "sprint" && test.str == "accord") {
  train.df <- SPRINT.df
  test.df <- ACCORD.df
} else if (train.str == "accord" && test.str == "sprint") {
  train.df <- ACCORD.df
  test.df <- SPRINT.df
}

# A function to both select input columns and rename them to something reasonable
select.input.vars <- function(x) {
  x %>% select(
    Age = AGE,
    Female = FEMALE,
    Black = RACE_BLACK,
    Hispanic = hisp,
    SBP = SBP.y,
    DBP = DBP.y,
    N_Medications = N_AGENTS,
    CurrentSmoker = currentsmoker,
    FormerSmoker = formersmoker,
    Aspirin = ASPIRIN,
    Statin = STATIN,
    Creatinine = SCREAT,
    Cholesterol = CHR,
    HDL_Cholesterol = HDL,
    Triglycerides = TRR,
    BMI = BMI,
    Diabetes = Diabetes
  )
}

# Prepare the train data
X.train <- select.input.vars(train.df)
Y.train <- round(train.df$t_cvds)  # Coarsen failure times for efficiency
W.train <- train.df$INTENSIVE
D.train <- train.df$cvd

# Prepare the test data
X.test <- select.input.vars(test.df)
Y.test <- round(test.df$t_cvds)  # Coarsen failure times for efficiency
W.test <- test.df$INTENSIVE
D.test <- test.df$cvd

# Standardize/coarsen the failure time grid across the two datasets
Y.max <- round(365.25 * 3)
failure.times <- 0:Y.max  # ACCORD-BP has some individuals that are censored on day 0.

#########################################################################
# TRAINING PRIORITIZATION RULES ON TRAIN, ESTIMATING PRIORITIES ON TEST #
#########################################################################

# Estimating Restricted Mean Survival Time (RMST)
estimate.rmst <- function(S, unique.times, end.time) {
  S <- S[, unique.times <= end.time]
  unique.times <- unique.times[unique.times <= end.time]
  time.diffs <- c(unique.times[1], diff(unique.times))
  time.diffs <- as.matrix(time.diffs)
  rmst <- S %*% time.diffs
  rmst
}

priorities.on.test <- data.frame(
  CSF = numeric(nrow(X.test)),
  RSF = numeric(nrow(X.test)),
  CoxPHSLearner = numeric(nrow(X.test)),
  ASCVD = numeric(nrow(X.test)),
  Framingham = numeric(nrow(X.test))
)

#------------------------#
# CAUSAL SURVIVAL FOREST #
#------------------------#

# Training a Causal Survival Forest model on train
train.model.CSF <- grf::causal_survival_forest(
  X = X.train,
  Y = Y.train,
  W = W.train,
  D = D.train,
  W.hat = rep(mean(W.train), length(W.train)),  # Use prop of treated individuals as E[W | X_i] for everyone (RCT)
  target = "RMST",
  horizon = Y.max,
  honesty = FALSE,
  failure.times = failure.times,
  alpha = 0.05,  # Default value
  tune.parameters = "all",
  seed = seed
)
# Generating predictions of the CSF train model on test
priorities.on.test$CSF <- predict(train.model.CSF, newdata=X.test)$predictions

#------------------------#
# RANDOM SURVIVAL FOREST #
#------------------------#

D.train.truncated <- D.train
Y.train.truncated <- Y.train
D.train.truncated[Y.train.truncated >= Y.max] <- 1
Y.train.truncated[Y.train.truncated >= Y.max] <- Y.max

# Training a Random Survival Forest model on train
train.model.RSF <- grf::survival_forest(
  X = X.train[W.train == 0,],
  Y = Y.train.truncated[W.train == 0],
  D = D.train.truncated[W.train == 0],
)

# Generating predictions of the RSF train model on test
RSF.pred.results <- predict(train.model.RSF, newdata=X.test)
priorities.on.test$RSF <- estimate.rmst(
  S = RSF.pred.results$predictions, 
  unique.times = RSF.pred.results$failure.times, 
  end.time = Y.max
) * -1  # We use negative of the RMST because we want to prioritize highest those individuals at greatest risk

#----------------------------------------------------#
# FRAMINGHAM HEART RISK SCORE & ASCVD RISK ESTIMATOR #
#----------------------------------------------------#

# Generate estimates of risk from Framingham Risk Score (simple) and 
# ACC/AHA 2013 ASCVD Risk Calculator
X.test.modified <- X.test
X.test.modified$Race <- ifelse(X.test$Black, "aa", "white")
X.test.modified$Gender <- ifelse(X.test$Female, "female", "male")

# For the purposes of the Framingham Heart Risk Score and 
# ASCVD Risk Calculator, we assume that subjects were not 
# on a blood pressure medication at the start of the trial
X.test.modified$BP_Medications <- 0

# Generate predictions of the Framingham 2008 10-year ASCVD risk model (with BMI)
library(CVrisk)
ascvd.frs.results <- CVrisk::compute_CVrisk(
  df = X.test.modified,
  age = "Age",
  gender = "Gender",
  race = "Race",
  sbp = "SBP",
  bmi = "BMI",
  hdl = "HDL_Cholesterol",
  totchol = "Cholesterol",
  bp_med = "BP_Medications",
  smoker = "CurrentSmoker",
  diabetes = "Diabetes"
)
priorities.on.test$ASCVD <- ascvd.frs.results$ascvd_10y_accaha
priorities.on.test$Framingham <- ascvd.frs.results$ascvd_10y_frs_simple

#--------------#
# COX PH MODEL # 
#--------------#

design.matrix.train <- X.train %>% select(-Diabetes)
design.matrix.train$Y <- Y.train
design.matrix.train$D <- D.train
design.matrix.train$W <- W.train
train.model.CPH.train <- survival::coxph(
  survival::Surv(time = Y, event = D) ~ . * W,
  data = design.matrix.train
)

design.matrix.test <- X.test %>% select(-Diabetes)
design.matrix.test.under.treatment <- design.matrix.test
design.matrix.test.under.treatment$W <- 1
design.matrix.test.under.control <- design.matrix.test
design.matrix.test.under.control$W <- 0

coxph.test.preds.treat <- predict(train.model.CPH.train, design.matrix.test.under.treatment)
coxph.test.preds.control <- predict(train.model.CPH.train, design.matrix.test.under.control)
priorities.on.test$CoxPHSLearner <- coxph.test.preds.treat - coxph.test.preds.control

###########################################
# ESTIMATING DOUBLY ROBUST SCORES ON TEST #
###########################################

cat("Training the forest to be used for evaluation...\n")
eval.model.CSF <- grf::causal_survival_forest(
  X = X.test,
  Y = Y.test,
  W = W.test,
  D = D.test,
  W.hat = rep(mean(W.test), length(W.test)),
  target = "RMST",
  horizon = Y.max,
  honesty = TRUE,
  failure.times = failure.times,
  alpha = 0.05,
  seed = seed
)

#############################################
# ESTIMATING THE RATE (AUTOC, QINI) ON TEST #
#############################################

tmp.n.rows <- 2 * ncol(priorities.on.test)
auc.results <- data.frame(
  prioritization_rule = character(tmp.n.rows),
  target = character(tmp.n.rows),
  point_estimate = numeric(tmp.n.rows),
  ci_lb = numeric(tmp.n.rows),
  ci_ub = numeric(tmp.n.rows),
  std_err = numeric(tmp.n.rows),
  p_value = numeric(tmp.n.rows),
  stringsAsFactors = FALSE
)

cat("Estimating the RATE under the learned prioritization rules...\n")
i <- 1
for(col in colnames(priorities.on.test)) {
  for(target in c("QINI", "AUTOC")) {
    tmp.RATE <- grf::rank_average_treatment_effect(
      forest = eval.model.CSF,
      target = target,
      priorities = priorities.on.test[, col],
      R = 1000
    )
    auc.results[i, "prioritization_rule"] <- col
    auc.results[i, "target"] <- target
    auc.results[i, "point_estimate"] <- tmp.RATE$estimate
    auc.results[i, "ci_lb"] <- tmp.RATE$estimate - 1.96 * tmp.RATE$std.err
    auc.results[i, "ci_ub"] <- tmp.RATE$estimate + 1.96 * tmp.RATE$std.err
    auc.results[i, "std_err"] <- tmp.RATE$std.err
    auc.results[i, "p_value"] <- pnorm(-abs(tmp.RATE$estimate) / tmp.RATE$std.err)
    i <- i + 1
  }
}

auc.results <- auc.results %>% arrange(target)  # Sort by AUTOC vs. QINI
auc.results %>% write_csv(paste0("train_on_", train.str, 
                                 "_test_on_", test.str,
                                 ".csv"))
