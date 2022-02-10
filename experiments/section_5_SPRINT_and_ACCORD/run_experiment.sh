#!/bin/bash

# Preprocess the SPRINT data
Rscript ../../data/accord/preprocess_accord.R

# Preprocess the ACCORD-BP data
Rscript ../../data/sprint/preprocess_sprint.R

# Train on SPRINT, test on ACCORD-BP
Rscript rate_on_SPRINT_and_ACCORD.R sprint accord

# Train on ACCORD-BP, test on SPRINT
Rscript rate_on_SPRINT_and_ACCORD.R accord sprint

# Generate tables from the output of the above
Rscript generate_table1_and_table2.R
