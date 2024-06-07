#!/usr/bin/env Rscript

library(tidyverse)

get.row.data <- function(prio.rule, tmp.df) {
  tmp.df[tmp.df$prioritization_rule == prio.rule & tmp.df$target == "AUTOC",]
}

gen.row.tex <- function(prio.rule.long, prio.rule.abbrev, res.df) {
  tmp.res <- get.row.data(prio.rule.abbrev, res.df)
  tmp.str <- ""
  tmp.str <- paste0(tmp.str, "\t", prio.rule.long, " & ")
  tmp.str <- paste0(tmp.str, round(tmp.res$point_estimate, 2), " ")
  tmp.str <- paste0(tmp.str, "(", round(tmp.res$ci_lb, 2), ", ")
  tmp.str <- paste0(tmp.str, round(tmp.res$ci_ub, 2), ")")
  tmp.str <- paste0(tmp.str, " & ", round(tmp.res$p_value, 2), " \\\\\n")
  tmp.str
}

gen.table.tex <- function(train.str, test.str) {
  fname <- paste0("train_on_", train.str, "_test_on_", test.str, ".csv")
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

writeLines(gen.table.tex("accord", "sprint"), "table1.tex")
writeLines(gen.table.tex("sprint", "accord"), "table2.tex")