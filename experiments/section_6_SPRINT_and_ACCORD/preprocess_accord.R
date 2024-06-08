# TODO: If not calling from the command line, set the working directory 
# according to your local environment. You should place in this same directory 
# the unzipped ACCORD_2017b_2 folder.
# setwd("/path/to/source_dir/")
library(tidyverse)

# Read in all the relevant datasets
tmp.prefix <- "data/ACCORD_2017b_2/Main_Study/3-Data_Sets-Analysis/3a-Analysis_Data_Sets/csv/"
accord.key <- read.csv(paste0(tmp.prefix, "accord_key.csv"))
blood.pressure <- read.csv(paste0(tmp.prefix, "bloodpressure.csv"))
concomitant.meds <- read.csv(paste0(tmp.prefix, "concomitantmeds.csv"))
lipids <- read.csv(paste0(tmp.prefix, "lipids.csv"))
other.labs <- read.csv(paste0(tmp.prefix, "otherlabs.csv"))
cvd.outcomes <- read.csv(paste0(tmp.prefix, "cvdoutcomes.csv"))

tmp.prefix <-  "data/ACCORD_2017b_2/Main_Study/4-Data_Sets-CRFs/4a-CRF_Data_Sets/csv/"
history <- read.csv(paste0(tmp.prefix, "f07_baselinehistoryphysicalexam.csv"))
incl.excl <- read.csv(paste0(tmp.prefix, "f01_inclusionexclusionsummary.csv"))

# Filter tables down to just baseline covariates (rather than follow-up visits)
blood.pressure <- blood.pressure %>% filter(Visit == "BLR")
concomitant.meds <- concomitant.meds %>% filter(Visit == "BLR")
lipids <- lipids %>% filter(Visit == "BLR")
other.labs <- other.labs %>% 
  filter(Visit == "BLR") %>%
  rename(potassium_lab = potassium)

# Filter subjects down to just those subjects who were assigned to the
# Intensive BP vs. Standard BP arms
accord.key <- accord.key %>% filter(
  treatment %in% c("Standard Gylcemia/Intensive BP",  # Arm 1
                   "Standard Gylcemia/Standard BP",  # Arm 2
                   "Intensive Gylcemia/Intensive BP",  # Arm 3
                   "Intensive Gylcemia/Standard BP")  # Arm 4
)

# Join all the datasets into a single dataframe
df <- accord.key %>% 
  left_join(blood.pressure, by="MaskID") %>% 
  left_join(concomitant.meds, by="MaskID") %>% 
  left_join(lipids, by="MaskID") %>%
  left_join(other.labs, by="MaskID") %>%
  left_join(cvd.outcomes, by="MaskID") %>%
  left_join(history, by="MaskID") %>%
  left_join(incl.excl, by="MaskID")

# Collect all the variables that are easy to derive/compute
df <- df %>% mutate(
  intensive = treatment %in% c(  # treatment = Treatment Assignment: description
    "Standard Gylcemia/Intensive BP",  # Arm 1
    "Intensive Gylcemia/Intensive BP"  # Arm 3
  ) * 1,
  age = baseline_age,
  female = female,
  black = (raceclass == "Black") * 1,
  hispanic = (raceclass == "Hispanic") * 1,
  SBP = sbp,
  DBP = dbp,
  aspirin = aspirin,
  statin = statin,
  creatinine = screat,
  cholesterol = chol,
  HDL_cholesterol = hdl,
  triglycerides = trig,
  BMI = wt_kg / ((ht_cm / 100)^2),
  eGFR = gfr,
  angina = (x2angina == 1) * 1,  # angina or ischemic changes on Graded Exercise Tolerance Test or positive imaging
  revascularization = (
    (cabg == 1) |  # CABG
    (ptci == 1) |  # PTCI/PTCA/Atherectomy
    (orevasc == 1)  # other revascularization procedure
  ) * 1
)

# Derive the number of anti-hypertensive medications prescribed
df <- df %>% mutate(
  BP_medications = loop + 
    thiazide + 
    ksparing + 
    potassium +
    a2rb + 
    acei + 
    dhp_ccb +
    nondhp_ccb + 
    alpha_blocker + 
    central_agent + 
    beta_blocker + 
    vasodilator + 
    reserpine + 
    other_bpmed
)

df <- df %>% mutate(
  current_smoker = (cigarett == 1) * 1,  # cigarett: smoked cigarettes in last 30 days
  former_smoker = ((smokelif == 1) & (cigarett != 1)) * 1  # smokelif: smoked more than 100 cigarettes during lifetime
)

# Derive the primary outcome for each subject
# The original primary outcome was a composite of (1) nonfatal MI, 
# (2) nonfatal stroke, and (3) death from cardiovascular causes. We expand
# this outcome to include "fatal or hospitalized congestive heart failure (CHF)"
# and "major coronary events (fatal coronary heart disease (CHD), non‐fatal 
# MI or unstable angina)" which is listed as "Major CHD" in the codebook.
# See Notes.CVDOutcomes.pdf for reference.
df <- df %>% mutate(
  event_primary = (
    (censor_nmi == 0) |  # Nonfatal MI censoring flag (0=event, 1=censored)
    (censor_nst == 0) |  # Nonfatal stroke censoring flag (0=event, 1=censored)
    (censor_cm == 0) |   # CVD mortality censoring flag (0=event, 1=censored)
    (censor_chf == 0) |  # CHF censoring flag (0=event, 1=censored)
    (censor_maj == 0)    # Major CHD censoring flag (0=event, 1=censored)
  ) * 1  
)

# Derive the censoring time for each subject
# See ACCORD Data Dictionary pg. 6 on cvdoutcomes.sas7bdat for definitions
df <- df %>% 
  rowwise() %>%
  mutate(
    t_primary = min(  
      fuyrs_nmi,  # fuyrs_nmi: Nonfatal MI Follow-up time (years)
      fuyrs_nst,  # fuyrs_nst: Nonfatal stroke Follow-up time (years)
      fuyrs_cm,   # fuyrs_cm: CVD mortality Follow-up time (years)
      fuyrs_chf,  # fuyrs_chf: CHF Follow-up time (years)
      fuyrs_maj   # fuyrs_maj: Major CHD Follow-up time (years)
    ) * 365.25
  )

df$diabetes <- 1

df <- df %>% 
  select(
    intensive,
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
    eGFR,
    angina,
    revascularization,
    diabetes,
    event_primary,
    t_primary
  )

df %>% drop_na() %>% write_csv("data/accord_cut.csv")