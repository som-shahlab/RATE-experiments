# Rank-Weighted Average Treatment Effects
This repository contains the code that was used to generate figures and tables for the paper, "[Evaluating Treatment Prioritization Rules via Rank-Weighted Average Treatment Effects](https://arxiv.org/abs/2111.07966)" by Steve Yadlowsky, Scott Fleming, Nigam Shah, Emma Brunskill, and Stefan Wager (2021).

## Code
### Abstract
Our [paper](https://arxiv.org/abs/2111.07966) proposes an approach for estimating Rank-weighted Average Treatment Effect (RATE) metrics on data. The empirical properties of this estimator are explored in simulations in Section 4 of the paper and its utility is highlighted in two applications in Section 5 and 6 of the paper. We provide an [implementation](https://grf-labs.github.io/grf/articles/rate.html) of the core RATE estimation machinery via an R package on CRAN ([grf](https://github.com/grf-labs/grf)) and provide a separate repository (this one) with R scripts that can be used to reproduce the results highlighted in the main manuscript.

### Description
The code for estimating the Rank-weighted Average Treatment Effect (RATE) on data is provided as part of the “generalized random forests” (grf) R package (https://github.com/grf-labs/grf) under a GPL-3 license.

Code relying on “grf” in order to generate the empirical results highlighted in sections 4, 5, and 6 of the main manuscript are accessible through a separate GitHub repository available to the public: https://github.com/som-shahlab/RATE-experiments. Code used to extract and preprocess specific features from the SPRINT and ACCORD dataset are given in the “data” folder, and scripts utilizing the “grf” package to analyze these data are provided in the “experiments” folder of the repository. All analysis code is available via a GPL-3 license. 

Experiments were performed on a single CPU, operating system macOS Big Sur. The code used to perform the analyses have a light dependence on the “tidyverse” R package v1.3.1.

## Instructions for Use
### Reproducibility
All tables and figures from the main manuscript can be reproduced from the analysis code repository (https://github.com/som-shahlab/RATE-experiments) by calling the wrapper script “run_experiment.sh” associated with each experiment subfolder. 

### Replication
Additional information regarding use of the “grf” package and RATE estimation specifically in the context of other package functions is provided via package documentation: https://grf-labs.github.io/grf/.

## Data
### Abstract
The experiments in Section 5 of our paper focus on applying our Rank-weighted Average Treatment Effect (RATE) approach to two randomized control clinical trials, SPRINT and ACCORD-BP. Conducted in 2001-2005 and 2010-2013, respectively, the ACCORD-BP and SPRINT trials both sought to identify whether intensive blood-pressure control could reduce the risk of adverse cardiovascular events (e.g., myocardial infarction, stroke, heart failure, death from cardiovascular causes) in individuals with elevated risk for cardiovascular disease. We used a small subset of variables collected in both trials to learn treatment prioritization rules and estimate their respective RATEs on trial data.

We additionally analyze data from a marketing application using the RATE. Criteo released a large benchmark dataset for studying uplift modeling in online digital advertising, based on anonymized results from a number of incrementality trials (Diemert et al., 2018). In combining trials, the interpretation of the treatment is subtle, corresponding to an intent to treat with one of a handful of arbitrarily chosen ads. The purpose of the data is to provide a benchmark for uplift modeling, and therefore, the results are not meant to be used in a particular application. See Diemert et al. (2018) for more information about the dataset construction and validation.


### Availability
Data for the SPRINT and ACCORD-BP clinical trials are available to researchers free of charge via the NHLBI Biologic Specimen and Data Repository Information Coordinating Center (BioLINCC, https://biolincc.nhlbi.nih.gov/home/). Researchers must submit an online data request form, including a study plan/protocol and Institutional Review Board (IRB) approval, to obtain access.

Data for the Criteo Uplift dataset can be found at https://ailab.criteo.com/criteo-uplift-prediction-dataset/.

Data for the International Stroke trial can be accessed via the following link:  
https://static-content.springer.com/esm/art%3A10.1186%2F1745-6215-12-101/MediaObjects/13063_2010_637_MOESM1_ESM.CSV
Data are provided as part of a public release via "The International Stroke Trial database" in _Trials_ by Sandercock et al, 2011.

### Description
SPRINT and ACCORD-BP data were accessed under an NHLBI Research Materials Distribution Agreement (RMDA), signed August 20th, 2019. The document is available upon request. Further information regarding data availability, access, and metadata for the SPRINT and ACCORD-BP clinical trials can be found at the trials’ NHLBI BioLINCC portal pages (https://biolincc.nhlbi.nih.gov/studies/sprint/ and https://biolincc.nhlbi.nih.gov/studies/accord/, respectively). The trial data used in our analysis were up to date as of January 1, 2020. 

According to the Criteo website (https://ailab.criteo.com/criteo-uplift-prediction-dataset/):

_Data description_  
This dataset is constructed by assembling data resulting from several incrementality tests, a particular randomized trial procedure where a random part of the population is prevented from being targeted by advertising. it consists of 25M rows, each one representing a user with 11 features, a treatment indicator and 2 labels (visits and conversions).  

_Privacy_  
For privacy reasons the data has been sub-sampled non-uniformly so that the original incrementality level cannot be deduced from the dataset while preserving a realistic, challenging benchmark. Feature names have been anonymized and their values randomly projected so as to keep predictive power while making it practically impossible to recover the original features or user context.  

_Fields_  
Here is a detailed description of the fields (they are comma-separated in the file):  
f0, f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11: feature values (dense, float)  
treatment: treatment group (1 = treated, 0 = control)  
conversion: whether a conversion occured for this user (binary, label)  
visit: whether a visit occured for this user (binary, label)  
exposure: treatment effect, whether the user has been effectively exposed (binary)  

## Notes
Contact the corresponding author for any additional information regarding data access and reproducibility.