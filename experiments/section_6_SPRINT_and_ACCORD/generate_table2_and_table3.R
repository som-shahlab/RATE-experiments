#!/usr/bin/env Rscript

library(tidyverse)

#' Extract data for a specific prioritization rule
#'
#' @param prio.rule The abbreviation of the prioritization rule
#' @param df The data frame containing the results
#' @return A subset of the data frame filtered by the specified prioritization rule and target "AUTOC"
get.row.data <- function(prio.rule, df) {
  df[df$`Prioritization Rule` == prio.rule & df$`RATE Metric` == "AUTOC",]
}

#' Generate a LaTeX table row for a given prioritization rule
#'
#' @param prio.rule.long The full name of the prioritization rule
#' @param prio.rule.abbrev The abbreviation of the prioritization rule
#' @param res.df The data frame containing the results
#' @return A string representing a LaTeX table row with the point estimate, confidence interval, and p-value
gen.row.tex <- function(prio.rule.long, prio.rule.abbrev, res.df) {
  tmp.res <- get.row.data(prio.rule.abbrev, res.df)
  tmp.str <- ""
  tmp.str <- paste0(tmp.str, "\t", prio.rule.long, " & ")
  tmp.str <- paste0(tmp.str, round(tmp.res$`RATE Point Estimate`, 2), " ")
  tmp.str <- paste0(tmp.str, "(", round(tmp.res$`CI Lower Bound`, 2), ", ")
  tmp.str <- paste0(tmp.str, round(tmp.res$`CI Upper Bound`, 2), ")")
  tmp.str <- paste0(tmp.str, " & ", round(tmp.res$`P-value`, 2), " \\\\\n")
  tmp.str
}

#' Generate the LaTeX code for a table
#'
#' @param train.str The name of the training dataset (either "sprint" or "accord")
#' @param test.str The name of the testing dataset (either "sprint" or "accord")
#' @return A string representing the LaTeX code for the table
#' @examples
#' gen.table.tex("accord", "sprint")
#' gen.table.tex("sprint", "accord")
gen.table.tex <- function(train.str, test.str, estimand = "RMST") {
  fname <- paste0(
    "train_on_", train.str, 
    "_test_on_", test.str,
    "_estimating_", estimand,
    ".csv"
  )
  res.df <- read_csv(fname)
  if (train.str == "sprint" && test.str == "accord") {
    n.train <- 9069
    n.test <- 4535
  } else if (train.str == "accord" && test.str == "sprint") {
    n.train <- 4535
    n.test <- 9069
  } else {
    stop("Train and/or test strings: ", train.str, " / ", test.str,
        "not defined. Only 'sprint' and 'accord' supported.")
  }

  tab <- ""
  tab <- paste0(tab, "\\begin{table}\n")
  tab <- paste0(tab, "\\centering\n")
  tab <- paste0(tab, "\\begin{tabular}{|r|cc|}\n")
  tab <- paste0(tab, "\t\\hline\n")
  tab <- paste0(tab, "\tPrioritization Rule & AUTOC (95\\% CI) & $p$-value \\\\\n")
  tab <- paste0(tab, "\t\\hline\n")
  tab <- paste0(tab, gen.row.tex("Causal Survival Forest (grf)", "CSF", res.df))
  tab <- paste0(tab, gen.row.tex("Cox PH S-learner", "CoxPHSLearner", res.df))
  tab <- paste0(tab, "\t\\hline\n")
  tab <- paste0(tab, gen.row.tex("Random Survival Forest Risk (grf)", "RSF", res.df))
  tab <- paste0(tab, gen.row.tex("Framingham Risk Score", "Framingham", res.df))
  tab <- paste0(tab, gen.row.tex("ACC/AHA Pooled Cohort Equations", "ASCVD", res.df))
  tab <- paste0(tab, "\t\\hline\n")
  tab <- paste0(tab, "\\end{tabular}\n")
  tab <- paste0(tab, "\\caption{AUTOC estimates obtained using data from ")
  tab <- paste0(tab, toupper(test.str), " ($n = ", n.test, "$),")
  tab <- paste0(tab, " with prioritization rules trained on ")
  tab <- paste0(tab, toupper(train.str), " ($n = ", n.train, "$),")
  tab <- paste0(tab, " if necessary. We also show 95\\% confidence intervals")
  tab <- paste0(tab, " obtained using the half-sample bootstrap, ")
  tab <- paste0(tab, " along with associated $p$-values.}\n")
  tab <- paste0(tab, "\\label{train_", train.str, "_test_", test.str, "}\n")
  tab <- paste0(tab, "\\end{table}")
}

# Generate LaTeX tables and write to files
writeLines(gen.table.tex("accord", "sprint"), "table2.tex")
writeLines(gen.table.tex("sprint", "accord"), "table3.tex")