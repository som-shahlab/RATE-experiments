#!/bin/bash

echo "============================================================\n\n"

# Preprocess the SPRINT data
echo "Preprocessing the SPRINT data...\n"
Rscript preprocess_sprint.R

echo "============================================================\n\n"

# Preprocess the ACCORD-BP data
echo "Preprocessing the ACCORD-BP data...\n"
Rscript preprocess_accord.R

echo "============================================================\n\n"

# Train on ACCORD-BP, test on SPRINT (Table 2 in the paper)
echo "Running with train=ACCORD-BP, test=SPRINT, estimand=RMST"
Rscript rate_on_SPRINT_and_ACCORD.R accord sprint RMST

echo "============================================================\n\n"

# Train on SPRINT, test on ACCORD-BP (Table 3 in the paper)
echo "Running with train=SPRINT, test=ACCORD-BP, estimand=RMST"
Rscript rate_on_SPRINT_and_ACCORD.R sprint accord RMST

echo "============================================================\n\n"

# Train on a 50/50 split of combined SPRINT/ACCORD-BP (Section D.3 in Supplement)
echo "Running with train=Combined, test=Combined, estimand=RMST"
Rscript rate_on_SPRINT_and_ACCORD.R combined combined RMST

# Generate LaTeX tables from the output of the above
Rscript generate_table2_and_table3.R
