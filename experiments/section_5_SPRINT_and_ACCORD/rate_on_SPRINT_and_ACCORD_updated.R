#!/usr/bin/env Rscript

# Install grf if not already installed
# devtools::install_github("erikcs/grf", subdir = "r-package/grf", ref = "RATE")
library(tidyverse, warn.conflicts=FALSE)

# If calling from the command line, assign train.str and test.str to user specs
args = commandArgs(trailingOnly=TRUE)
train.str <- tolower(args[1])
test.str <- tolower(args[2])
estimand <- args[3]  # Should be "survival.probability" or "RMST"

# If not calling from the command line, set manually
# setwd("~/Documents/GitHub/RATE-experiments/experiments/section_5_SPRINT_and_ACCORD/")
# train.str <- "accord"
# test.str <- "sprint"
# estimand <- "survival.probability"
# train.str <- "combined"
# test.str <- "combined"
# estimand <- "survival.probability"

seed <- 42
set.seed(seed)  # Required to get consistent results with the RATE
SPRINT.df <- read_csv("../../data/sprint/sprint_cut_updated.csv", show_col_types=FALSE)
ACCORD.df <- read_csv("../../data/accord/accord_cut_updated.csv", show_col_types=FALSE)

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

# A function to both select input columns and rename them to something reasonable
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
    # eGFR,
    # angina,
    # revascularization,
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

#################################################
# ESTIMATING SURVIVAL CURVE NUISANCE PARAMETERS #
#################################################

# The conditional survival function for the censoring process S_C(t, x, w).
C.hat.train <- NULL
C.hat.test <- NULL
W.hat.train <- rep(mean(W.train), length(W.train))  # Assuming cross-trial prediction
W.hat.test <- rep(mean(W.test), length(W.test))
if (train.str == "combined" && test.str == "combined") {
  
  # # Generate nuisance parameter estimates of the conditional censoring distribution for each trial
  # args.nuisance <- list(failure.times = 0:round(365.25 * 3),
  #                       num.trees = 2000,
  #                       honesty = TRUE,
  #                       honesty.fraction = 0.5,
  #                       honesty.prune.leaves = TRUE,
  #                       prediction.type = "Nelson-Aalen", # to guarantee non-zero estimates.
  #                       compute.oob.predictions = TRUE,
  #                       seed = seed)
  # X.comb <- rbind(X.train, X.test)
  # X.comb$W <- c(W.train, W.test)
  # Y.comb <- c(Y.train, Y.test)
  # D.comb <- c(D.train, D.test)
  # # Note: Diabetes is an exact proxy for whether sample is from SPRINT/ACCORD-BP
  # sf.censor.SPRINT <- do.call(
  #   grf::survival_forest, 
  #   c(
  #     list(X=X.comb[X.comb$Diabetes == 0,], 
  #          Y = Y.comb[X.comb$Diabetes == 0], 
  #          D = 1 - D.comb[X.comb$Diabetes == 0]), 
  #     args.nuisance
  #   )
  # )
  # sf.censor.ACCORD <- do.call(
  #   grf::survival_forest, 
  #   c(
  #     list(X=X.comb[X.comb$Diabetes == 1,], 
  #          Y = Y.comb[X.comb$Diabetes == 1], 
  #          D = 1 - D.comb[X.comb$Diabetes == 1]), 
  #     args.nuisance
  #   )
  # )
  # C.hat.SPRINT <- predict(sf.censor.SPRINT, failure.times = failure.times)$predictions
  # C.hat.ACCORD <- predict(sf.censor.ACCORD, failure.times = failure.times)$predictions
  # C.hat <- matrix(NA, nrow=nrow(X.comb), ncol=ncol(C.hat.SPRINT))
  # C.hat[X.comb$Diabetes == 0,] <- C.hat.SPRINT
  # C.hat[X.comb$Diabetes == 1,] <- C.hat.ACCORD
  # C.hat.train <- C.hat[train.indxs,]
  # C.hat.test <- C.hat[-train.indxs,]
  
  # Generate nuisance parameter estimates of the probability of receiving treatment
  # W.hat.SPRINT <- mean(SPRINT.df$INTENSIVE)
  # W.hat.ACCORD <- mean(ACCORD.df$INTENSIVE)
  # W.hat <- numeric(nrow(combined.df))
  # W.hat[X.comb$Diabetes == 0] <- W.hat.SPRINT
  # W.hat[X.comb$Diabetes == 1] <- W.hat.ACCORD
}

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
  # C.hat = C.hat.train,
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
  # tune.parameters="all"
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
X.test.modified$race <- ifelse(X.test$black, "aa", "white")
X.test.modified$gender <- ifelse(X.test$female, "female", "male")

# For the purposes of the Framingham Heart Risk Score and 
# ASCVD Risk Calculator, we assume that subjects were not 
# on a blood pressure medication at the start of the trial
X.test.modified$BP_Medications <- (X.test.modified$BP_medications > 0) * 1

source("risk_calculators.R")

adapted.ascvd.frs.results <- adapted_compute_CVrisk(
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
priorities.on.test$Framingham <- ascvd.frs.results$ascvd_10y_frs_simple
priorities.on.test$Framingham.w.Labs <- ascvd.frs.results$ascvd_10y_frs

priorities.on.test$Extrapolated.ASCVD <- adapted.ascvd.frs.results$adapted_ascvd_10y_accaha
priorities.on.test$Extrapolated.Framingham <- adapted.ascvd.frs.results$adapted_ascvd_10y_frs_simple
priorities.on.test$Extrapolated.Framingham.w.Labs <- adapted.ascvd.frs.results$adapted_ascvd_10y_frs
# View(cbind(test.df[, c("black", "female", "age", "cholesterol", "HDL_cholesterol", "SBP", "DBP", "BP_medications", "current_smoker", "diabetes")], priorities.on.test))

# Check against https://tools.acc.org/ldl/ascvd_risk_estimator/index.html
# library(CVrisk)
# CVrisk::ascvd_10y_accaha(
#   race = "aa",
#   gender = "male",
#   age = 70,
#   totchol = 220,
#   hdl = 50,
#   sbp = 150,
#   bp_med = 0,
#   smoker = 1,
#   diabetes = 1
# )

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
  # C.hat = C.hat.test,
  target = estimand,
  horizon = Y.max,
  honesty = TRUE,
  failure.times = failure.times,
  alpha = 0.05,
  seed = seed
)

# ate <- grf::average_treatment_effect(eval.model.CSF)
# grf::average_treatment_effect(eval.model.CSF, subset=X.test$diabetes == 0)
# grf::average_treatment_effect(eval.model.CSF, subset=X.test$diabetes == 1)

#############################################
# ESTIMATING THE RATE (AUTOC, QINI) ON TEST #
#############################################

auc.results <- data.frame(
  matrix(,nrow=2*ncol(priorities.on.test),ncol=0),
  stringsAsFactors = FALSE
)

# TODO: Pull request to add censorship weights explicitly
# TODO: Pull request to allow for NaN-valued priorities
cat("Estimating the RATE under the learned prioritization rules...\n")
i <- 1
for(col in colnames(priorities.on.test)) {
  for(target in c("QINI", "AUTOC")) {
    tmp.RATE <- grf::rank_average_treatment_effect(
      forest = eval.model.CSF,
      target = target,
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

auc.results <- auc.results %>% arrange("Prioritization Rule")  # Sort by AUTOC vs. QINI
auc.results %>% write_csv(paste0("train_on_", train.str, 
                                 "_test_on_", test.str,
                                 "_estimating_", estimand,
                                 "_updated.csv"))

# tmp.RATE <- grf::rank_average_treatment_effect(
#   forest = eval.model.CSF,
#   target = "AUTOC",
#   q = seq(0.01, 1, 0.005),
#   priorities = priorities.on.test[, "Framingham"],
#   # subset = which(priorities.on.test[, col] > -999999999999999999999999),
#   subset = which(!is.na(priorities.on.test[, "Framingham"])),
#   R = 1000
# )
# plot(
#   tmp.RATE, 
#   main = sprintf(
#     "TOC - Framingham\nRATE = %.3f +/- %.3f, P = %.3f", 
#     tmp.RATE$estimate, 
#     1.96 * tmp.RATE$std.err, 
#     pnorm(-abs(tmp.RATE$estimate)/tmp.RATE$std.err))
# )
# 
# tmp.RATE <- grf::rank_average_treatment_effect(
#   forest = eval.model.CSF,
#   target = "AUTOC",
#   q = seq(0.005, 1, 0.005),
#   priorities = priorities.on.test[, "Framingham.w.Labs"],
#   # subset = which(priorities.on.test[, col] > -999999999999999999999999),
#   subset = which(!is.na(priorities.on.test[, "Framingham.w.Labs"])),
#   R = 1000
# )
# plot(
#   tmp.RATE, 
#   main = sprintf(
#     "TOC - Framingham with Labs\nRATE = %.3f +/- %.3f, P = %.3f", 
#     tmp.RATE$estimate, 
#     1.96 * tmp.RATE$std.err, 
#     2*pnorm(-abs(tmp.RATE$estimate)/tmp.RATE$std.err))
# )
# 
# tmp.RATE <- grf::rank_average_treatment_effect(
#   forest = eval.model.CSF,
#   target = "AUTOC",
#   priorities = priorities.on.test[, "CSF"],
#   # subset = which(priorities.on.test[, col] > -999999999999999999999999),
#   subset = which(!is.na(priorities.on.test[, "CSF"])),
#   R = 1000
# )
# plot(
#   tmp.RATE, 
#   main = sprintf(
#     "TOC - Causal Survival Forest\nRATE = %.3f +/- %.3f, P = %.3f", 
#     tmp.RATE$estimate, 
#     1.96 * tmp.RATE$std.err, 
#     2*pnorm(-abs(tmp.RATE$estimate)/tmp.RATE$std.err))
# )
# 
# tmp.RATE <- grf::rank_average_treatment_effect(
#   forest = eval.model.CSF,
#   target = "AUTOC",
#   q = seq(0.005, 1, 0.005),
#   priorities = priorities.on.test[, "CSF"],
#   # subset = which(priorities.on.test[, col] > -999999999999999999999999),
#   subset = which(!is.na(priorities.on.test[, "CSF"])),
#   R = 1000
# )
# plot(
#   tmp.RATE, 
#   main = sprintf(
#     "TOC - Causal Survival Forest\nTrain on %s, Test on %s\nRATE = %.3f +/- %.3f, P = %.3f", 
#     toupper(train.str),
#     toupper(test.str),
#     tmp.RATE$estimate, 
#     1.96 * tmp.RATE$std.err, 
#     2*pnorm(-abs(tmp.RATE$estimate)/tmp.RATE$std.err)
#   )
# )
# 
# length(priorities.on.test$CSF)
# sum(tmp.RATE$TOC > 0.2)
# 
# reso <- 1200
# length <- 3.25*reso/72
# 
# png(
#   paste0("hist_of_Framingham_simple_train_on_",train.str,"_test_on_",test.str,"_updated.png"),
#   res=reso,
#   height=length,
#   width=2*length
# )
# hist(
#   # priorities.on.test$Framingham[priorities.on.test$Framingham > -999999999999999999999999], 
#   priorities.on.test$Extrapolated.Framingham[!is.na(priorities.on.test$Extrapolated.Framingham)],
#   main=paste0("Framingham Scores (Simple) on ", test.str, " (test set)"),
#   xlab="Framingham Score (Simple)"
# )
# dev.off()
# 
# png(paste0("hist_of_Framingham_with_labs_train_on_",train.str,"_test_on_",test.str,"_updated.png"))
# hist(
#   # priorities.on.test$Framingham.w.Labs[priorities.on.test$Framingham.w.Labs > -999999999999999999999999], 
#   main=paste0("Framingham Scores w/ Labs on ", test.str, " (test set)"),
#   xlab="Framingham Score w/ Labs"
# )
# dev.off()
# 
# png(paste0("hist_of_ASCVD_train_on_",train.str,"_test_on_",test.str,".png"))
# hist(
#   # priorities.on.test$ASCVD[priorities.on.test$ASCVD > -999999999999999999999999], 
#   main=paste0("ASCVD Scores w/ Labs on ", test.str, " (test set)"),
#   xlab="ASCVD Score w/ Labs"
# )
# dev.off()
# 
# # priorities.on.test[priorities.on.test == -999999999999999999999999] <- NA
# 
# # blp <- grf::best_linear_projection(train.model.CSF, X.train)
# # blp
# 
# var.imp.df <- data.frame(
#   feature = colnames(X.train),
#   importance = grf::variable_importance(train.model.CSF)
# )
# var.imp.df <- var.imp.df %>% arrange(desc(importance))
# var.imp.df %>% write_csv(paste0("feat_importances_train_on_", train.str, 
#                       "_test_on_", test.str,
#                       "_estimating_", estimand,
#                       ".csv"))
# 
# # dr.scores <- grf::get_scores(eval.model.CSF)
# # dr.scores.by.ntiles <- data.frame(
# #   dr_scores = dr.scores,
# #   n_tile_CSF = factor(ntile(priorities.on.test$CSF, n=n.tiles)),
# #   n_tile_Framingham = factor(ntile(priorities.on.test$Framingham, n=n.tiles))
# # )
# # ggplot(dr.scores.by.ntiles) + 
# #   geom_point(aes(x=n_tile_CSF, y=dr_scores))
# 
# n.tiles <- 5
# eval.by.n.tile <- function(priorities, n.tiles = 4) {
#   estimated.aipw.ate <- lapply(
#     seq(n.tiles),
#     function(w) {
#       ate <- grf::average_treatment_effect(
#         eval.model.CSF,
#         subset = factor(ntile(priorities, n=n.tiles)) == w
#       )
#     }
#   )
#   data.frame(do.call(rbind, estimated.aipw.ate))
# }
# ate.by.ntiles <- rbind(
#   data.frame(eval.by.n.tile(priorities.on.test$CSF, n.tiles = n.tiles), Model="CSF", n_tile = seq(1, n.tiles)),
#   data.frame(eval.by.n.tile(priorities.on.test$Framingham, n.tiles = n.tiles), Model="Framingham", n_tile = seq(1, n.tiles)),
#   data.frame(eval.by.n.tile(priorities.on.test$Framingham.w.Labs, n.tiles = n.tiles), Model="Framingham.w.Labs", n_tile = seq(1, n.tiles)),
#   data.frame(eval.by.n.tile(priorities.on.test$ASCVD, n.tiles = n.tiles), Model="ASCVD", n_tile = seq(1, n.tiles))
# )
# ate.by.ntiles$ci.lb <- ate.by.ntiles$estimate + qnorm(0.025) * ate.by.ntiles$std.err
# ate.by.ntiles$ci.ub <- ate.by.ntiles$estimate + qnorm(0.975) * ate.by.ntiles$std.err
# 
# # ate.by.ntiles <- rbind(
# #   data.frame(eval.by.n.tile(X.test$HDL_Cholesterol), Model="HDL Cholesterol", n_tile = seq(1, 4)),
# #   data.frame(eval.by.n.tile(X.test$Triglycerides), Model="Triglycerides", n_tile = seq(1, 4))
# # )
# # for(m in unique(ate.by.ntiles$Model)) {
# g <- ggplot(ate.by.ntiles, aes(x=n_tile, y=estimate, group=Model, color=Model)) +
#   geom_errorbar(
#     aes(ymin=estimate - 1.96 * std.err, ymax=estimate + 1.96 * std.err), 
#     width=0.1, position=position_dodge(0.3)
#   ) + 
#   geom_hline(
#     aes(yintercept=ate["estimate"], color="ATE"), 
#     linetype="dashed",
#     # color = "grey"
#   ) + 
#   # geom_line(position = position_dodge(0.3)) + 
#   geom_point(position = position_dodge(0.3)) + 
#   geom_ribbon(
#     aes(
#       # x = seq_len(15) * (5/16),
#       ymin = ate["estimate"] - 1.96 * ate["std.err"],
#       ymax = ate["estimate"] + 1.96 * ate["std.err"],
#       color = "ATE",
#       # fill = "ATE",
#     ),
#     # fill = "grey",
#     # color = "grey",
#     linetype = 0,
#     alpha = 0.01
#   ) +
#   xlab("Quintile by Model Priority") + 
#   ylab(paste("ATE in Quintile\n(", estimand, "at 3 Years)")) +
#   theme_bw()
# ggsave(filename = paste0("ate_by_ntiles_train_on_", train.str, "_test_on_", test.str, "_estimating_", estimand, ".png"), 
#        plot = g, 
#        device = "png",
#        dpi=320)
# 



