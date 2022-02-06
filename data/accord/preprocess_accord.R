library(sas7bdat)
library(doBy)  # Needed for "summaryBy"
library(matrixStats)  # Needed for "rowMaxs"
library(survival)  # NEeded for "Surv"
# library(Hmisc)

# TODO: Set this according to your local environment. You should
# place in this folder the unzipped ACCORD_2017b_2 folder.
path_to_accord_data <- "~/Documents/GitHub/risk-vs-hte/data/accord/"
setwd(
  paste0(
    path_to_accord_data, 
    "ACCORD_2017b_2/Main_Study/3-Data_Sets-Analysis/3a-Analysis_Data_Sets/"
  )
)

##########################
# SANJAY'S ORIGINAL CODE ##
##########################

# Copied from https://github.com/sanjaybasu/sprint-challenge/blob/master/accordorg.R

accord_key = read.sas7bdat("accord_key.sas7bdat")
accord_key = accord_key[(accord_key$arm>=1)&(accord_key$arm<=4),]
bloodpressure = read.sas7bdat("bloodpressure.sas7bdat")
concomitantmeds = read.sas7bdat("concomitantmeds.sas7bdat")
cvdoutcomes = read.sas7bdat("cvdoutcomes.sas7bdat")
lipids = read.sas7bdat("lipids.sas7bdat")
otherlabs = read.sas7bdat("otherlabs.sas7bdat")

setwd(paste0(path_to_accord_data, "ACCORD_2017b_2/Main_Study/4-Data_Sets-CRFs/4a-CRF_Data_Sets"))
f13_intensivebpmanagement = read.sas7bdat("f13_intensivebpmanagement.sas7bdat")
f14_standardbpmanagement = read.sas7bdat("f14_standardbpmanagement.sas7bdat")
f07_baselinehistoryphysicalexam = read.sas7bdat("f07_baselinehistoryphysicalexam.sas7bdat")
accord_key_cut = accord_key[which((accord_key$arm==1)|(accord_key$arm==2)|(accord_key$arm==3)|(accord_key$arm==4)),]
accord_key_cut$INTENSIVE = (accord_key_cut$arm==1)|(accord_key_cut$arm==3)
bloodpressure_cut = bloodpressure[which(bloodpressure$Visit=="BLR"),]
concomitantmeds_cut = concomitantmeds[which(concomitantmeds$Visit=="BLR"),]
f13_intensivebpmanagement_cut = f13_intensivebpmanagement[which(f13_intensivebpmanagement$Visit=="BLR"|
                                                                  f13_intensivebpmanagement$Visit=="F01"|
                                                                  f13_intensivebpmanagement$Visit=="F02"|
                                                                  f13_intensivebpmanagement$Visit=="F03"|
                                                                  f13_intensivebpmanagement$Visit=="F04"|
                                                                  f13_intensivebpmanagement$Visit=="F05"|
                                                                  f13_intensivebpmanagement$Visit=="F06"|
                                                                  f13_intensivebpmanagement$Visit=="F07"|
                                                                  f13_intensivebpmanagement$Visit=="F08"|
                                                                  f13_intensivebpmanagement$Visit=="F09"|
                                                                  f13_intensivebpmanagement$Visit=="F10"|
                                                                  f13_intensivebpmanagement$Visit=="F11"|
                                                                  f13_intensivebpmanagement$Visit=="F12"|
                                                                  f13_intensivebpmanagement$Visit=="F13"|
                                                                  f13_intensivebpmanagement$Visit=="F14"|
                                                                  f13_intensivebpmanagement$Visit=="F15"|
                                                                  f13_intensivebpmanagement$Visit=="F16"|
                                                                  f13_intensivebpmanagement$Visit=="F17"|
                                                                  f13_intensivebpmanagement$Visit=="F18"|
                                                                  f13_intensivebpmanagement$Visit=="F19"|
                                                                  f13_intensivebpmanagement$Visit=="F20"|
                                                                  f13_intensivebpmanagement$Visit=="F21"|
                                                                  f13_intensivebpmanagement$Visit=="F22"|
                                                                  f13_intensivebpmanagement$Visit=="F23"|
                                                                  f13_intensivebpmanagement$Visit=="F24"|
                                                                  f13_intensivebpmanagement$Visit=="F25"|
                                                                  f13_intensivebpmanagement$Visit=="F26"|
                                                                  f13_intensivebpmanagement$Visit=="F27"|
                                                                  f13_intensivebpmanagement$Visit=="F28"|
                                                                  f13_intensivebpmanagement$Visit=="F29"|
                                                                  f13_intensivebpmanagement$Visit=="F30"|
                                                                  f13_intensivebpmanagement$Visit=="F31"|
                                                                  f13_intensivebpmanagement$Visit=="F32"|
                                                                  f13_intensivebpmanagement$Visit=="F33"|
                                                                  f13_intensivebpmanagement$Visit=="F34"|
                                                                  f13_intensivebpmanagement$Visit=="F35"|
                                                                  f13_intensivebpmanagement$Visit=="F36"|
                                                                  f13_intensivebpmanagement$Visit=="F37"|
                                                                  f13_intensivebpmanagement$Visit=="F38"|
                                                                  f13_intensivebpmanagement$Visit=="F39"|
                                                                  f13_intensivebpmanagement$Visit=="F40"|
                                                                  f13_intensivebpmanagement$Visit=="F41"|
                                                                  f13_intensivebpmanagement$Visit=="F42"|
                                                                  f13_intensivebpmanagement$Visit=="F43"|
                                                                  f13_intensivebpmanagement$Visit=="F44"|
                                                                  f13_intensivebpmanagement$Visit=="F45"|
                                                                  f13_intensivebpmanagement$Visit=="F46"|
                                                                  f13_intensivebpmanagement$Visit=="F47"|
                                                                  f13_intensivebpmanagement$Visit=="F48"|
                                                                  f13_intensivebpmanagement$Visit=="F49"|
                                                                  f13_intensivebpmanagement$Visit=="F50"|
                                                                  f13_intensivebpmanagement$Visit=="F51"|
                                                                  f13_intensivebpmanagement$Visit=="F52"|
                                                                  f13_intensivebpmanagement$Visit=="F53"|
                                                                  f13_intensivebpmanagement$Visit=="F54"|
                                                                  f13_intensivebpmanagement$Visit=="F55"|
                                                                  f13_intensivebpmanagement$Visit=="F56"|
                                                                  f13_intensivebpmanagement$Visit=="F57"|
                                                                  f13_intensivebpmanagement$Visit=="F58"|
                                                                  f13_intensivebpmanagement$Visit=="F59"|
                                                                  f13_intensivebpmanagement$Visit=="F60"),]
f13_intensivebpmanagement_cut$visvec = as.numeric(f13_intensivebpmanagement_cut$Visit)
f13_intensivebpmanagement_cut$visnum = 0*(f13_intensivebpmanagement_cut$visvec==1)+
  1*(f13_intensivebpmanagement_cut$visvec==3)+
  2*(f13_intensivebpmanagement_cut$visvec==4)+
  3*(f13_intensivebpmanagement_cut$visvec==5)+
  4*(f13_intensivebpmanagement_cut$visvec==6)+
  6*(f13_intensivebpmanagement_cut$visvec==7)+
  8*(f13_intensivebpmanagement_cut$visvec==8)+
  12*(f13_intensivebpmanagement_cut$visvec==10)+
  16*(f13_intensivebpmanagement_cut$visvec==12)+
  20*(f13_intensivebpmanagement_cut$visvec==14)+
  24*(f13_intensivebpmanagement_cut$visvec==16)+
  28*(f13_intensivebpmanagement_cut$visvec==18)+
  30*(f13_intensivebpmanagement_cut$visvec==19)+
  32*(f13_intensivebpmanagement_cut$visvec==20)+
  34*(f13_intensivebpmanagement_cut$visvec==21)+
  36*(f13_intensivebpmanagement_cut$visvec==22)+
  38*(f13_intensivebpmanagement_cut$visvec==23)+
  40*(f13_intensivebpmanagement_cut$visvec==24)+
  42*(f13_intensivebpmanagement_cut$visvec==25)+
  44*(f13_intensivebpmanagement_cut$visvec==26)+
  46*(f13_intensivebpmanagement_cut$visvec==27)+
  48*(f13_intensivebpmanagement_cut$visvec==28)+
  50*(f13_intensivebpmanagement_cut$visvec==29)+
  52*(f13_intensivebpmanagement_cut$visvec==30)+
  54*(f13_intensivebpmanagement_cut$visvec==31)+
  56*(f13_intensivebpmanagement_cut$visvec==32)+
  58*(f13_intensivebpmanagement_cut$visvec==33)+
  60*(f13_intensivebpmanagement_cut$visvec==34)+
  10*(f13_intensivebpmanagement_cut$visvec==9)+
  14*(f13_intensivebpmanagement_cut$visvec==11)+
  18*(f13_intensivebpmanagement_cut$visvec==13)+
  22*(f13_intensivebpmanagement_cut$visvec==15)+
  26*(f13_intensivebpmanagement_cut$visvec==17)
f13_intensivebpmanagement_cut$visadv = f13_intensivebpmanagement_cut$visnum*f13_intensivebpmanagement_cut$advexp
f13_intensivebpmanagement_cut1 = summaryBy(advexp~MaskID,data = f13_intensivebpmanagement_cut, na.rm=TRUE,FUN=max)
f13_intensivebpmanagement_cut1$advexp.max[f13_intensivebpmanagement_cut1$advexp.max=="-Inf"]=""
f13_intensivebpmanagement_cut2 = summaryBy(visadv~MaskID,data = f13_intensivebpmanagement_cut, na.rm=TRUE,FUN=min)
f13_intensivebpmanagement_cut2$visadv.min[f13_intensivebpmanagement_cut2$visadv.min=="Inf"]=""
f13_intensivebpmanagement_cut3 = summaryBy(visnum~MaskID,data = f13_intensivebpmanagement_cut, na.rm=TRUE,FUN=max)
f13_intensivebpmanagement_cut = merge(f13_intensivebpmanagement_cut1,f13_intensivebpmanagement_cut2,by="MaskID")
f13_intensivebpmanagement_cut = merge(f13_intensivebpmanagement_cut,f13_intensivebpmanagement_cut3,by="MaskID")
f14_standardbpmanagement_cut = f14_standardbpmanagement[which(f14_standardbpmanagement$Visit=="BLR"|
                                                                f14_standardbpmanagement$Visit=="F01"|
                                                                f14_standardbpmanagement$Visit=="F02"|
                                                                f14_standardbpmanagement$Visit=="F03"|
                                                                f14_standardbpmanagement$Visit=="F04"|
                                                                f14_standardbpmanagement$Visit=="F05"|
                                                                f14_standardbpmanagement$Visit=="F06"|
                                                                f14_standardbpmanagement$Visit=="F07"|
                                                                f14_standardbpmanagement$Visit=="F08"|
                                                                f14_standardbpmanagement$Visit=="F09"|
                                                                f14_standardbpmanagement$Visit=="F10"|
                                                                f14_standardbpmanagement$Visit=="F11"|
                                                                f14_standardbpmanagement$Visit=="F12"|
                                                                f14_standardbpmanagement$Visit=="F13"|
                                                                f14_standardbpmanagement$Visit=="F14"|
                                                                f14_standardbpmanagement$Visit=="F15"|
                                                                f14_standardbpmanagement$Visit=="F16"|
                                                                f14_standardbpmanagement$Visit=="F17"|
                                                                f14_standardbpmanagement$Visit=="F18"|
                                                                f14_standardbpmanagement$Visit=="F19"|
                                                                f14_standardbpmanagement$Visit=="F20"|
                                                                f14_standardbpmanagement$Visit=="F21"|
                                                                f14_standardbpmanagement$Visit=="F22"|
                                                                f14_standardbpmanagement$Visit=="F23"|
                                                                f14_standardbpmanagement$Visit=="F24"|
                                                                f14_standardbpmanagement$Visit=="F25"|
                                                                f14_standardbpmanagement$Visit=="F26"|
                                                                f14_standardbpmanagement$Visit=="F27"|
                                                                f14_standardbpmanagement$Visit=="F28"|
                                                                f14_standardbpmanagement$Visit=="F29"|
                                                                f14_standardbpmanagement$Visit=="F30"|
                                                                f14_standardbpmanagement$Visit=="F31"|
                                                                f14_standardbpmanagement$Visit=="F32"|
                                                                f14_standardbpmanagement$Visit=="F33"|
                                                                f14_standardbpmanagement$Visit=="F34"|
                                                                f14_standardbpmanagement$Visit=="F35"|
                                                                f14_standardbpmanagement$Visit=="F36"|
                                                                f14_standardbpmanagement$Visit=="F31"|
                                                                f14_standardbpmanagement$Visit=="F32"|
                                                                f14_standardbpmanagement$Visit=="F33"|
                                                                f14_standardbpmanagement$Visit=="F34"|
                                                                f14_standardbpmanagement$Visit=="F35"|
                                                                f14_standardbpmanagement$Visit=="F36"|
                                                                f14_standardbpmanagement$Visit=="F37"|
                                                                f14_standardbpmanagement$Visit=="F38"|
                                                                f14_standardbpmanagement$Visit=="F39"|
                                                                f14_standardbpmanagement$Visit=="F40"|
                                                                f14_standardbpmanagement$Visit=="F41"|
                                                                f14_standardbpmanagement$Visit=="F42"|
                                                                f14_standardbpmanagement$Visit=="F43"|
                                                                f14_standardbpmanagement$Visit=="F44"|
                                                                f14_standardbpmanagement$Visit=="F45"|
                                                                f14_standardbpmanagement$Visit=="F46"|
                                                                f14_standardbpmanagement$Visit=="F47"|
                                                                f14_standardbpmanagement$Visit=="F48"|
                                                                f14_standardbpmanagement$Visit=="F49"|
                                                                f14_standardbpmanagement$Visit=="F50"|
                                                                f14_standardbpmanagement$Visit=="F51"|
                                                                f14_standardbpmanagement$Visit=="F52"|
                                                                f14_standardbpmanagement$Visit=="F53"|
                                                                f14_standardbpmanagement$Visit=="F54"|
                                                                f14_standardbpmanagement$Visit=="F55"|
                                                                f14_standardbpmanagement$Visit=="F56"|
                                                                f14_standardbpmanagement$Visit=="F57"|
                                                                f14_standardbpmanagement$Visit=="F58"|
                                                                f14_standardbpmanagement$Visit=="F59"|
                                                                f14_standardbpmanagement$Visit=="F60"),]
f14_standardbpmanagement_cut$visvec = as.numeric(f14_standardbpmanagement_cut$Visit)
f14_standardbpmanagement_cut$visnum = 0*(f14_standardbpmanagement_cut$visvec==1)+
  1*(f14_standardbpmanagement_cut$visvec==3)+
  2*(f14_standardbpmanagement_cut$visvec==4)+
  3*(f14_standardbpmanagement_cut$visvec==5)+
  4*(f14_standardbpmanagement_cut$visvec==6)+
  6*(f14_standardbpmanagement_cut$visvec==7)+
  8*(f14_standardbpmanagement_cut$visvec==8)+
  12*(f14_standardbpmanagement_cut$visvec==9)+
  16*(f14_standardbpmanagement_cut$visvec==10)+
  20*(f14_standardbpmanagement_cut$visvec==11)+
  24*(f14_standardbpmanagement_cut$visvec==12)+
  28*(f14_standardbpmanagement_cut$visvec==13)+
  32*(f14_standardbpmanagement_cut$visvec==14)+
  36*(f14_standardbpmanagement_cut$visvec==15)+
  40*(f14_standardbpmanagement_cut$visvec==16)+
  44*(f14_standardbpmanagement_cut$visvec==17)+
  48*(f14_standardbpmanagement_cut$visvec==18)+
  52*(f14_standardbpmanagement_cut$visvec==19)+
  56*(f14_standardbpmanagement_cut$visvec==20)+
  60*(f14_standardbpmanagement_cut$visvec==21)
f14_standardbpmanagement_cut$visadv = f14_standardbpmanagement_cut$visnum*f14_standardbpmanagement_cut$advexp
f14_standardbpmanagement_cut1 = summaryBy(advexp~MaskID,data = f14_standardbpmanagement_cut, na.rm=TRUE,FUN=max)
f14_standardbpmanagement_cut1$advexp.max[f14_standardbpmanagement_cut1$advexp.max=="-Inf"]=""
f14_standardbpmanagement_cut2 = summaryBy(visadv~MaskID,data = f14_standardbpmanagement_cut, na.rm=TRUE,FUN=min)
f14_standardbpmanagement_cut2$visadv.min[f14_standardbpmanagement_cut2$visadv.min=="Inf"]=""
f14_standardbpmanagement_cut3 = summaryBy(visnum~MaskID,data = f14_standardbpmanagement_cut, na.rm=TRUE,FUN=max)
f14_standardbpmanagement_cut = merge(f14_standardbpmanagement_cut1,f14_standardbpmanagement_cut2,by="MaskID")
f14_standardbpmanagement_cut = merge(f14_standardbpmanagement_cut,f14_standardbpmanagement_cut3,by="MaskID")
bpmg = rbind(f13_intensivebpmanagement_cut,f14_standardbpmanagement_cut)
lipids_cut = lipids[which(lipids$Visit=="BLR"),]
otherlabs_cut = otherlabs[which(otherlabs$Visit=="BLR"),]
accord_set = merge(accord_key_cut,bloodpressure_cut,by="MaskID")
accord_set = merge(accord_set,concomitantmeds_cut,by=c("MaskID","Visit"))
accord_set = merge(accord_set,bpmg,by="MaskID")
accord_set = merge(accord_set,f07_baselinehistoryphysicalexam,by="MaskID")
accord_set = merge(accord_set,cvdoutcomes,by="MaskID")
accord_set = merge(accord_set,lipids_cut,by="MaskID")
accord_set = merge(accord_set,otherlabs_cut,by="MaskID")


###################################
# TONY'S VERSION OF SANJAY'S CODE #
###################################

# Copied from https://github.com/tonyduan/hte-prediction-rcts/blob/master/data/accord/preprocess_accord.R

cvd = (accord_set$censor_nmi==0)|(accord_set$censor_nst==0)|(accord_set$censor_cm==0)|(accord_set$censor_chf==0)|(accord_set$censor_maj==0)
t_censor = rowMaxs(cbind(accord_set$fuyrs_nmi*365.25,accord_set$fuyrs_nst*365.25,accord_set$fuyrs_cm*365.25,accord_set$fuyrs_chf*365.25,accord_set$fuyrs_maj*365.25))
t_cvds = rowMaxs(cbind(accord_set$fuyrs_nmi*365.25*(1-accord_set$censor_nmi),accord_set$fuyrs_nst*365.25*(1-accord_set$censor_nst),accord_set$fuyrs_cm*365.25*(1-accord_set$censor_cm),accord_set$fuyrs_chf*365.25*(1-accord_set$censor_chf),accord_set$fuyrs_maj*365.25*(1-accord_set$censor_maj)))
t_cvds[t_cvds==0] = t_censor[t_cvds==0]
t_cvds[t_cvds==0] = 'NA'
t_cvds = as.numeric(t_cvds)
cOutcome = Surv(time=t_cvds, event = cvd)
sae = (accord_set$advexp.max==1)
accord_set$visadv.min=as.numeric(accord_set$visadv.min)
t_censor = rowMaxs(cbind(accord_set$visnum.max*30.42))
t_saes = rowMaxs(cbind(accord_set$visadv.min*30.42))
t_saes[is.na(t_saes)] = t_censor[is.na(t_saes)]
t_saes[t_saes==0] = 'NA'
t_saes = as.numeric(t_saes)
dOutcome = Surv(time=t_saes, event = sae)
INTENSIVE = as.numeric(accord_set$INTENSIVE)
AGE = accord_set$baseline_age
FEMALE = accord_set$female
RACE_BLACK = as.numeric(accord_set$raceclass=="Black")
hisp = (accord_set$raceclass=="Hispanic")
SBP.y = accord_set$sbp
DBP.y = accord_set$dbp
N_AGENTS = (accord_set$loop+accord_set$thiazide+accord_set$ksparing+accord_set$a2rb+accord_set$acei+accord_set$dhp_ccb+accord_set$nondhp_ccb+accord_set$alpha_blocker+accord_set$central_agent+accord_set$beta_blocker+accord_set$vasodilator+accord_set$reserpine+accord_set$other_bpmed)
currentsmoker = (accord_set$cigarett==1)
formersmoker = (accord_set$smokelif==1)
ASPIRIN= accord_set$aspirin
STATIN = accord_set$statin
SUB_SENIOR = as.numeric(AGE>=75)
SUB_CKD = as.numeric(accord_set$gfr<60)
CHR = accord_set$chol
GLUR = accord_set$fpg
HDL = accord_set$hdl
TRR = accord_set$trig
UMALCR = accord_set$uacr
EGFR = accord_set$gfr
SCREAT = accord_set$screat
BMI = accord_set$wt_kg/((accord_set$ht_cm/1000)^2)/100
c2<-data.frame(cvd,t_cvds,
               INTENSIVE,AGE,FEMALE,RACE_BLACK,hisp,
               SBP.y,DBP.y,N_AGENTS,currentsmoker,formersmoker,
               ASPIRIN,STATIN,
               SCREAT,CHR,HDL,TRR,BMI)
c2=c2[complete.cases(c2),]

setwd(path_to_accord_data)
write.csv(c2, "./accord_cut.csv", row.names=FALSE)


##########################
# SANJAY'S ORIGINAL CODE #
##########################

# cvd = (accord_set$censor_nmi==0)|(accord_set$censor_nst==0)|(accord_set$censor_cm==0)|(accord_set$censor_chf==0)|(accord_set$censor_maj==0)
# t_censor = rowMaxs(cbind(accord_set$fuyrs_nmi*365.25,accord_set$fuyrs_nst*365.25,accord_set$fuyrs_cm*365.25,accord_set$fuyrs_chf*365.25,accord_set$fuyrs_maj*365.25))
# t_cvds = rowMaxs(cbind(accord_set$fuyrs_nmi*365.25*(1-accord_set$censor_nmi),accord_set$fuyrs_nst*365.25*(1-accord_set$censor_nst),accord_set$fuyrs_cm*365.25*(1-accord_set$censor_cm),accord_set$fuyrs_chf*365.25*(1-accord_set$censor_chf),accord_set$fuyrs_maj*365.25*(1-accord_set$censor_maj)))
# t_cvds[t_cvds==0] = t_censor[t_cvds==0]
# t_cvds[t_cvds==0] = 'NA'
# t_cvds = as.numeric(t_cvds)
# cOutcome = Surv(time=t_cvds, event = cvd)
# sae = (accord_set$advexp.max==1)
# accord_set$visadv.min=as.numeric(accord_set$visadv.min)
# t_censor = rowMaxs(cbind(accord_set$visnum.max*30.42))
# t_saes = rowMaxs(cbind(accord_set$visadv.min*30.42))
# t_saes[is.na(t_saes)] = t_censor[is.na(t_saes)]
# t_saes[t_saes==0] = 'NA'
# t_saes = as.numeric(t_saes)
# dOutcome = Surv(time=t_saes, event = sae)
# INTENSIVE = as.numeric(accord_set$INTENSIVE)
# AGE = accord_set$baseline_age
# FEMALE = accord_set$female
# RACE_BLACK = as.numeric(accord_set$raceclass=="Black")
# hisp = (accord_set$raceclass=="Hispanic")
# SBP.y = accord_set$sbp
# DBP.y = accord_set$dbp
# N_AGENTS = (accord_set$loop+accord_set$thiazide+accord_set$ksparing+accord_set$a2rb+accord_set$acei+accord_set$dhp_ccb+accord_set$nondhp_ccb+accord_set$alpha_blocker+accord_set$central_agent+accord_set$beta_blocker+accord_set$vasodilator+accord_set$reserpine+accord_set$other_bpmed)
# currentsmoker = (accord_set$cigarett==1)
# formersmoker = (accord_set$smokelif==1)
# ASPIRIN= accord_set$aspirin
# STATIN = accord_set$statin
# SUB_SENIOR = as.numeric(AGE>=75)
# SUB_CKD = as.numeric(accord_set$gfr<60)
# CHR = accord_set$chol
# GLUR = accord_set$fpg
# HDL = accord_set$hdl
# TRR = accord_set$trig
# UMALCR = accord_set$uacr
# EGFR = accord_set$gfr
# SCREAT = accord_set$screat
# BMI = accord_set$wt_kg/((accord_set$ht_cm/1000)^2)/100