# TODO: If not calling from the command line, set this according to your 
# local environment. You should place in this directory the 
# unzipped SPRINT_2020b folder.
# setwd("/path/to/source_dir/")
library(tidyverse)

# See SPRINT_2020b/SPRINT-POP/documentation/codebook.xlsx and 
# SPRINT_2020b/SPRINT-POP/SPRINT-POP Data Dictionary.pdf 
# for details onwhat the variable values actually mean
baseline.df <- read.csv("data/SPRINT_2020b/SPRINT-POP/data/baseline.csv")
outcomes.df <- read.csv("data/SPRINT_2020b/SPRINT-POP/data/outcomes.csv")
history.df <- read.csv("data/SPRINT_2020b/SPRINT/data/CSV/bl_history.csv")
incl.excl.df <- read.csv("data/SPRINT_2020b/SPRINT/data/CSV/incl_excl.csv")

# We use fields from SPRINT-POP except when available only in SPRINT
history.df <- history.df %>% select(-ASPIRIN)  # Use ASPIRIN from SPRINT-POP
incl.excl.df <- incl.excl.df %>% select(-HDL)  # Use HDL from SPRINT-POP

# Merge all the dataframes into a single dataframe 
df <- baseline.df %>% 
  left_join(outcomes.df, by="MASKID") %>%
  left_join(history.df, by="MASKID") %>%
  left_join(incl.excl.df, by="MASKID")

# To compare SPRINT and ACCORD-BP we need to slightly modify
# the primary outcome for sprint, specifically we need to 
# exclude "other acute coronary syndromes", listed as EVENT_NONMIACS
# in the codebook (see abstract of the original manuscript)
df <- df %>% 
  rowwise() %>%
  mutate(  # Compare to `EVENT_PRIMARY` (ours excludes EVENT_NONMIACS)
    event_primary = (
      EVENT_MI |  # Myocardial infarction outcome (1=event, 0=censored)
        EVENT_STROKE |  # Stroke outcome (1=event, 0=censored)
        EVENT_HF |  # Heart failure outcome (1=event, 0=censored)
        EVENT_CVDDEATH  # CVD death outcome (1=event, 0=censored)
    ) * 1,  
    t_primary = min(  # T_MI (etc.) is the censoring time if EVENT_MI (etc.) == 0
      T_MI,  # Elapsed time for myocardial infarction (days)
      T_STROKE,  # Elapsed time for stroke (days)
      T_CVDDEATH,  # Elapsed time for CVD death (days)
      T_HF  # Elapsed time for heart failure (days)
    )
  )

df <- df %>% 
  transmute(
    intensive = INTENSIVE,  # Assigned to intensive BP arm
    age = AGE,  # Derived: Age at randomization top-coded at 90 years
    female = FEMALE,  # Derived: Female Gender (0/1 indicator of female gender)
    black = (RACE4 == "BLACK") * 1,  # Incl/Excl: self-reported race/ethnicity
    hispanic = (RACE4 == "HISPANIC") * 1,  # if Hispanic ethnicity then value is "HISPANIC" all other values are non-Hispanic
    SBP = SBP,  # Derived: Seated Systolic Blood Pressure (mm Hg)
    DBP = DBP,  # Derived: Seated Diastolic Blood Pressure (mm Hg)
    BP_medications = N_AGENTS,  # Derived: Number of anti-hypertensive medications prescribed (if 0 then participants on no anti-hypertensive agents)
    # SMOKE_3CAT: Baseline smoking status = 1 - Never;  2 - Former; 3 - Current; 4 - Missing
    # (We use this formulation, but ACCORD-BP uses an alternative derived variable)
    current_smoker = (SMOKE_3CAT == 3) * 1,
    former_smoker = (SMOKE_3CAT == 2) * 1,
    # current_smoker = !is.na(NOWSMOKE) & (NOWSMOKE == 1) * 1,
    # former_smoker = (SMOKED100 == 1) & (NOWSMOKE == 0) * 1,
    aspirin = ASPIRIN,  # BSL Hist: Daily Aspirin Use
    statin = STATIN,  # Derived: on any statin
    creatinine = SCREAT,  # Lab: Serum creatinine (mg/dL)
    cholesterol = CHR,  # Lab: Cholesterol, mg/dL
    HDL_cholesterol = HDL,  # Lab: HDL-cholesterol direct, mg/dL
    triglycerides = TRR,  # Lab: Triglycerides, mg/dL
    BMI = BMI,  # Derived: Body Mass Index (kg/m^2)
    eGFR = case_when(!is.na(EGFR) ~ EGFR,  # EGFR is from SPRINT-POP (also see RESULT_GFR)
                     is.na(EGFR) ~ GFRESTIMATE),  # GFRESTIMATE is estimated gfr within past 6 months
    angina = ANGINA,  # Have you ever been told by a physician that you have: Blockage in the blood flow to your heart, also called coronary artery disease?
    revascularization = CORONARYREVAS,  # CVD History: Coronary Revascularization (CABG, PCI)
    diabetes = 0,
    event_primary,
    t_primary
  )

df %>% drop_na() %>% write_csv("data/sprint_cut.csv")