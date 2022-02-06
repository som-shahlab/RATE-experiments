#!/bin/bash

# Train on SPRINT, test on ACCORD-BP
Rscript rate_on_SPRINT_and_ACCORD.R sprint accord

# Train on ACCORD-BP, test on SPRINT
Rscript rate_on_SPRINT_and_ACCORD.R accord sprint
