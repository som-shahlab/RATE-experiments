# ===========================================================================
# Code adapted from https://github.com/sanjaybasu/sprint-challenge.
# ===========================================================================

# ==== load relevant libraries

library(Hmisc)
library(matrixStats)
library(sas7bdat)

# ==== load relevant data

# TODO: If not calling from the command line, set this according to your 
# local environment. You should place in this folder the 
# unzipped SPRINT_2019a folder.
setwd("~/Documents/GitHub/RATE-experiments/data/sprint/")

##########################
# SANJAY'S ORIGINAL CODE #
##########################

# Copied from https://github.com/sanjaybasu/sprint-challenge/blob/master/sprintorg.R
# Copyright (c) 2020 Sanjay Basu, under MIT License

baseline = read.sas7bdat("SPRINT_2020b/SPRINT-POP/data/baseline.sas7bdat")
bp = read.sas7bdat("SPRINT_2020b/SPRINT-POP/data/bp.sas7bdat")
outcomes = read.sas7bdat("SPRINT_2020b/SPRINT-POP/data/outcomes.sas7bdat")
retention = read.sas7bdat("SPRINT_2020b/SPRINT-POP/data/retention.sas7bdat")
safety = read.sas7bdat("SPRINT_2020b/SPRINT-POP/data/safety.sas7bdat")
bp_cut = bp[which(bp$VISITCODE=="RZ"),]
sprint_set = merge(baseline,bp_cut,by="MASKID")
sprint_set = merge(sprint_set,outcomes,by="MASKID")
sprint_set = merge(sprint_set,retention,by="MASKID")
sprint_set = merge(sprint_set,safety,by="MASKID")

####################################
# TONY'S ADDITION TO SANJAY'S CODE #
####################################

# TODO: Remove the attach (for better style and less ambiguity)
attach(sprint_set)

hisp = (RACE4 == "HISPANIC")
currentsmoker = (SMOKE_3CAT == 3)
formersmoker = (SMOKE_3CAT == 2)
cvd = (EVENT_MI == 1) | (EVENT_STROKE == 1) | (EVENT_CVDDEATH == 1) | (EVENT_HF == 1)
t_censor = rowMaxs(cbind(T_MI,T_STROKE,T_CVDDEATH,T_HF))
t_cvds = rowMaxs(cbind(T_MI*EVENT_MI,T_STROKE*EVENT_STROKE,T_CVDDEATH*EVENT_CVDDEATH,T_HF*EVENT_HF))
t_cvds[t_cvds==0] = t_censor[t_cvds==0]
t_cvds[t_cvds==0] = 'NA'
t_cvds = as.numeric(t_cvds)
cOutcome = Surv(time=t_cvds, event = cvd)
sae = (HYP_SAE_EVNT==1)|(SYN_SAE_EVNT==1)|(ELE_SAE_EVNT==1)|(AKI_SAE_EVNT==1)|(BRA_SAE_EVNT==1)
t_censor = rowMaxs(cbind(HYP_SAE_DAYS,SYN_SAE_DAYS,ELE_SAE_DAYS,AKI_SAE_DAYS,BRA_SAE_DAYS))
t_saes = rowMaxs(cbind(HYP_SAE_DAYS*HYP_SAE_EVNT,SYN_SAE_DAYS*SYN_SAE_EVNT,ELE_SAE_DAYS*ELE_SAE_EVNT,AKI_SAE_DAYS*AKI_SAE_EVNT,BRA_SAE_DAYS*BRA_SAE_EVNT))
t_saes[t_saes==0] = t_censor[t_saes==0]
t_saes[t_saes==0] = 'NA'
t_saes = as.numeric(t_saes)
dOutcome = Surv(time=t_saes, event = sae)
testsubset = data.frame(cOutcome,
                        INTENSIVE,AGE,FEMALE,RACE_BLACK,hisp,
                        SBP.y,DBP.y,N_AGENTS,currentsmoker,formersmoker,
                        ASPIRIN,STATIN,
                        SCREAT,CHR,HDL,TRR,BMI,
                        INTENSIVE*AGE,INTENSIVE*FEMALE,INTENSIVE*RACE_BLACK,INTENSIVE*hisp,
                        INTENSIVE*SBP.y,INTENSIVE*DBP.y,INTENSIVE*N_AGENTS,INTENSIVE*currentsmoker,INTENSIVE*formersmoker,
                        INTENSIVE*ASPIRIN,INTENSIVE*STATIN,
                        INTENSIVE*SCREAT,INTENSIVE*CHR,INTENSIVE*HDL,INTENSIVE*TRR,INTENSIVE*BMI)
testsubset=testsubset[complete.cases(testsubset),]

# For descriptions of variables see SPRINT_2019a/SPRINT-POP/codebook.xlsx
c<-data.frame(
    cvd,
    t_cvds,
    INTENSIVE,  # Assigned to intensive BP arm
    AGE,  #
    FEMALE,  # Derived: Female Gender (0/1 indicator of female gender)
    RACE_BLACK,  # Incl/Excl: 0/1 Indicator of African American race, by self-report
    hisp,  # Derived: Patient is hispanic (by self-reported race/ethnicity, px may also be white/black/other if hisp)
    SBP.y,  # Derived: Seated Systolic Blood Pressure (mm Hg)
    DBP.y,  # Derived: Seated Diastolic Blood Pressure (mm Hg)
    N_AGENTS,  # Number of medications prescribed
    currentsmoker,  # Derived: Patient is a current smoker
    formersmoker,  # Derived: Patient is a former smoker
    ASPIRIN,  # BSL Hist: Daily Aspirin Use
    STATIN,  # Derived: on any statin
    SCREAT,  # Lab: Serum creatinine (mg/dL)
    CHR,  # Lab: Cholesterol, mg/dL
    HDL,  # Lab: HDL-cholesterol direct, mg/dL
    TRR,  # Lab: Triglycerides, mg/dL
    BMI   # Derived: Body Mass Index (kg/m^2)
)
c=c[complete.cases(c),]

write.csv(c, "sprint_cut.csv", row.names=FALSE)
