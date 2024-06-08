#!/usr/bin/env Rscript

# Install grf if not already installed
# devtools::install_github("erikcs/grf", subdir = "r-package/grf", ref = "RATE")
library(tidyverse, warn.conflicts=FALSE)
library(CVrisk)

# If calling from the command line, assign train.str and test.str to user specs
args = commandArgs(trailingOnly=TRUE)
train.str <- tolower(args[1])
test.str <- tolower(args[2])
estimand <- args[3]  # Should be "survival.probability" or "RMST"

# If not calling from the command line, set manually as follows:
# train.str <- "combined"
# test.str <- "combined"
# estimand <- "survival.probability"

seed <- 42
set.seed(seed)  # Required to get consistent results with the RATE
SPRINT.df <- read_csv("data/sprint_preprocessed.csv", show_col_types=FALSE)
ACCORD.df <- read_csv("data/accord_preprocessed.csv", show_col_types=FALSE)

# Set the train/test data according to command-line argument
if (train.str == "sprint" && test.str == "accord") {
  train.df <- SPRINT.df
  test.df <- ACCORD.df
} else if (train.str == "accord" && test.str == "sprint") {
  train.df <- ACCORD.df
  test.df <- SPRINT.df
} else if (train.str == "combined" && test.str == "combined") {
  combined.df <- rbind(SPRINT.df, ACCORD.df)
  n.combined <- nrow(combined.df)
  indxs <- 1:n.combined
  train.indxs <- sample(indxs, size=round(n.combined / 2), replace=FALSE)
  train.df <- combined.df[train.indxs,]
  test.df <- combined.df[-train.indxs,]
}

#' Select and Rename Input Variables
#'
#' This function selects a predefined set of input columns from a data frame 
#' and renames them to more reasonable names if necessary.
#'
#' @param x A data frame containing the variables to be selected.
#' @return A data frame with the selected columns
select.input.vars <- function(x) {
  x %>% select(
    age,
    female,
    black,
    hispanic,
    SBP,
    DBP,
    BP_medications,
    current_smoker,
    former_smoker,
    aspirin,
    statin,
    creatinine,
    cholesterol,
    HDL_cholesterol,
    triglycerides,
    BMI,
    diabetes
  )
}

# Prepare the train data
X.train <- select.input.vars(train.df)
Y.train <- round(train.df$t_primary)  # Coarsen failure times for efficiency
W.train <- train.df$intensive
D.train <- train.df$event_primary

# Prepare the test data
X.test <- select.input.vars(test.df)
Y.test <- round(test.df$t_primary)  # Coarsen failure times for efficiency
W.test <- test.df$intensive
D.test <- test.df$event_primary

# Standardize/coarsen the failure time grid across the two datasets
Y.max <- round(365.25 * 3)
failure.times <- 0:Y.max  # ACCORD-BP has some individuals that are censored on day 0.

#########################################################################
# TRAINING PRIORITIZATION RULES ON TRAIN, ESTIMATING PRIORITIES ON TEST #
#########################################################################

#' Estimating Restricted Mean Survival Time (RMST)
#'
#' Computes the Restricted Mean Survival Time (RMST) up to a specified end time.
#'
#' @param S A matrix of survival probabilities, where each row represents an individual and each column represents a time point.
#' @param unique.times A vector of unique time points corresponding to the columns of the survival matrix S.
#' @param end.time The end time up to which the RMST is calculated.
#' @return A vector of RMST values for each individual in the survival matrix S.
#' @examples
#' S <- matrix(c(1, 0.8, 0.6, 1, 0.9, 0.7), nrow = 2, byrow = TRUE)
#' unique.times <- c(1, 2, 3)
#' end.time <- 2
#' estimate.rmst(S, unique.times, end.time)
#' # [1] 1.8 1.9
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
  CoxPHSLearner = numeric(nrow(X.test)),
  RSF = numeric(nrow(X.test)),
  Framingham.w.Labs = numeric(nrow(X.test)),
  ASCVD = numeric(nrow(X.test))
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
  target = estimand,
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
) * -1  # We use negative RMST because we want to prioritize highest those individuals at greatest risk

#----------------------------------------------------#
# FRAMINGHAM HEART RISK SCORE & ASCVD RISK ESTIMATOR #
#----------------------------------------------------#

# Generate estimates of risk from Framingham Risk Score (simple) and 
# ACC/AHA 2013 ASCVD Risk Calculator
X.test.modified <- X.test
X.test.modified$race <- ifelse(X.test$black, "aa", "white")
X.test.modified$gender <- ifelse(X.test$female, "female", "male")

# See https://github.com/vcastro/CVrisk/blob/master/R/ascvd_10y_frs.R#LL91C45-L91C45
X.test.modified$BP_Medications <- (X.test.modified$BP_medications > 0) * 1

# Generate predictions of the Framingham 2008 10-year ASCVD risk model (with BMI)
ascvd.frs.results <- CVrisk::compute_CVrisk(
  df = X.test.modified,
  age = "age",
  gender = "gender",
  race = "race",
  sbp = "SBP",
  bmi = "BMI",
  hdl = "HDL_cholesterol",
  totchol = "cholesterol",
  bp_med = "BP_medications",
  smoker = "current_smoker",
  diabetes = "diabetes"
)
priorities.on.test$ASCVD <- ascvd.frs.results$ascvd_10y_accaha
priorities.on.test$Framingham.w.Labs <- ascvd.frs.results$ascvd_10y_frs

#--------------#
# COX PH MODEL # 
#--------------#

design.matrix.train <- X.train
if (train.str != "combined") {
  design.matrix.train <- design.matrix.train %>% select(-diabetes)
}
design.matrix.train$Y <- Y.train
design.matrix.train$D <- D.train
design.matrix.train$W <- W.train
train.model.CPH.train <- survival::coxph(
  survival::Surv(time = Y, event = D) ~ . * W,
  data = design.matrix.train
)

design.matrix.test <- X.test
if (train.str != "combined") {
  # Otherwise `diabetes` column is all 1's or all 0's
  design.matrix.test <- design.matrix.test %>% select(-diabetes)
}
design.matrix.test.under.treatment <- design.matrix.test
design.matrix.test.under.treatment$W <- 1
design.matrix.test.under.control <- design.matrix.test
design.matrix.test.under.control$W <- 0

coxph.test.preds.treat <- predict(train.model.CPH.train, design.matrix.test.under.treatment)
coxph.test.preds.control <- predict(train.model.CPH.train, design.matrix.test.under.control)
priorities.on.test$CoxPHSLearner <- coxph.test.preds.treat - coxph.test.preds.control

#-----------------#
# TRIAL INDICATOR #
#-----------------#

if(train.str == "combined") {
  priorities.on.test$TrialIndicator <- X.test$diabetes == 0
}

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
  target = estimand,
  horizon = Y.max,
  honesty = TRUE,
  failure.times = failure.times,
  alpha = 0.05,
  seed = seed
)

#############################################
# ESTIMATING THE RATE (AUTOC, QINI) ON TEST #
#############################################

auc.results <- data.frame(
  matrix(,nrow=2*ncol(priorities.on.test),ncol=0),
  stringsAsFactors = FALSE
)

cat("Estimating the RATE under the learned prioritization rules...\n")
i <- 1
for(col in colnames(priorities.on.test)) {
  for(target in c("QINI", "AUTOC")) {
    tmp.RATE <- grf::rank_average_treatment_effect(
      forest = eval.model.CSF,
      target = target,  # AUTOC vs. QINI
      priorities = priorities.on.test[, col],
      subset = which(!is.na(priorities.on.test[, col])),
      R = 1000
    )
    auc.results[i, "Prioritization Rule"] <- col
    auc.results[i, "Train Set"] <- train.str
    auc.results[i, "Test Set"] <- test.str
    auc.results[i, "Estimand"] <- estimand
    auc.results[i, "RATE Metric"] <- target
    auc.results[i, "RATE Point Estimate"] <- tmp.RATE$estimate
    auc.results[i, "CI Lower Bound"] <- tmp.RATE$estimate - 1.96 * tmp.RATE$std.err
    auc.results[i, "CI Upper Bound"] <- tmp.RATE$estimate + 1.96 * tmp.RATE$std.err
    auc.results[i, "RATE Std Err"] <- tmp.RATE$std.err
    auc.results[i, "P-value"] <- 2 * pnorm(-abs(tmp.RATE$estimate) / tmp.RATE$std.err)
    i <- i + 1
  }
}

hist(priorities.on.test$CSF)
hist(priorities.on.test$RSF)

auc.results <- auc.results %>% arrange("Prioritization Rule")  # Sort by AUTOC vs. QINI
auc.results[auc.results$`RATE Metric` == "AUTOC",]
auc.results %>% write_csv(paste0("train_on_", train.str, 
                                 "_test_on_", test.str,
                                 "_estimating_", estimand,
                                 ".csv"))