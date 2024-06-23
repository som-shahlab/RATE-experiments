# Rank-Weighted Average Treatment Effects

This repository contains the code that was used to generate figures and tables for the paper, "[Evaluating Treatment Prioritization Rules via Rank-Weighted Average Treatment Effects](https://arxiv.org/abs/2111.07966)" by Steve Yadlowsky, Scott Fleming, Nigam Shah, Emma Brunskill, and Stefan Wager (2021).

## Table of Contents
- [Summary](#summary)
- [Core Rank-weighted Average Treatment Effect (`RATE`) implementation](#core-rank-weighted-average-treatment-effect-rate-implementation)
- [Applications using the `RATE`](#applications-using-the-rate)
- [Computational Environment and Packages](#computational-environment-and-packages)
- [Reproducibility](#reproducibility)
- [Repository Files and Structure](#repository-files-and-structure)
  - [`experiments/`](#experiments)
    - [`experiments/section_4_comparing_weighting_functions`](#experimentssection_4_comparing_weighting_functions)
    - [`experiments/section_5_international_stroke_trial`](#experimentssection_5_international_stroke_trial)
    - [`experiments/section_6_SPRINT_and_ACCORD`](#experimentssection_6_sprint_and_accord)
    - [`experiments/section_7_digital_marketing/`](#experimentssection_7_digital_marketing)
  - [Additional Repository Files](#additional-repository-files)
- [Data](#data)
  - [Section 4: Comparing RATE vs. Qini Coefficient using Synthetic Data](#section-4-comparing-rate-vs-qini-coefficient-using-synthetic-data)
  - [Section 5: The International Stroke Trial (Clinical, no censoring)](#section-5-the-international-stroke-trial-clinical-no-censoring)
  - [Section 6: The SPRINT and ACCORD-BP Trials (Clinical, right-censoring)](#section-6-the-sprint-and-accord-bp-trials-clinical-right-censoring)
  - [Section 7: The Criteo Uplift Dataset (Marketing)](#section-7-the-criteo-uplift-dataset-marketing)
- [Notes](#notes)


## Summary

Our [paper](https://arxiv.org/abs/2111.07966) proposes an approach for estimating Rank-weighted Average Treatment Effect (RATE) metrics on data. The empirical properties of this estimator are explored in simulations in Section 4 of the paper and its utility is highlighted in three applications in Sections 5, 6, and 7 of the paper. We provide an [implementation](https://grf-labs.github.io/grf/articles/rate.html) of the core RATE estimation machinery via an R package on CRAN ([grf](https://github.com/grf-labs/grf)) and provide a separate repository (this one) with R scripts that can be used to reproduce the results highlighted in the main manuscript.

## Core Rank-weighted Average Treatment Effect (`RATE`) implementation

The code for estimating the Rank-weighted Average Treatment Effect (RATE) on data is provided as part of the `generalized random forests` (`grf`) R package (https://github.com/grf-labs/grf) under a GPL-3 license.

Additional information regarding use of the `grf` package and RATE estimation specifically in the context of other package functions is provided via package documentation: https://grf-labs.github.io/grf/.

## Applications using the `RATE`

<!-- Code relying on `grf` in order to generate the empirical results highlighted in sections 4, 5, and 6 of the main manuscript are accessible through a separate GitHub repository available to the public: https://github.com/som-shahlab/RATE-experiments.  -->
Scripts utilizing the `grf` package to analyze data with RATE for each experiment are provided in the `experiments` folder of the repository. All analysis code is available via a GPL-3 license.

## Computational Environment and Packages

Experiments were performed on a single CPU, operating system macOS Big Sur. The code used to perform the analyses have a light dependence on the “tidyverse” R package v1.3.1. Simulations utilizing Python 3.9.7 have lightweight requirements listed in the `requirements.txt` file in the main directory of the repository.

## Reproducibility

All tables and figures from the main manuscript can be reproduced from this analysis code repository (https://github.com/som-shahlab/RATE-experiments) by running the R Markdown notebook (`.Rmd`) or calling the wrapper script `run_experiment.sh` associated with each experiment subfolder. More details below.

## Repository Files and Structure

### `experiments/`

The `experiments/` folder contains all the code necessary to reproduce the results from Sections 4, 5, 6, and 7 in our manuscript. Each subdirectory corresponds to a specific section of the manuscript and includes experiment-specific scripts, data processing, and analysis code.

#### `experiments/section_4_comparing_weighting_functions`

Contains code for generating simulated data and comparing various weighting functions used in our analysis. Invoking `./run_experiment.sh` calls `qini_vs_autoc.py` and generates figures within the `figures/` subdirectory which comprise Figure 2 in Section 4 of the main manuscript.

#### `experiments/section_5_international_stroke_trial`

The file `analysis.Rmd` walks through, in a more didactic and narrative format, our analysis of the International Stroke Trial. The markdown notebook generates Figure 1 in the Introduction and data populating Table 1 in Section 5 of the main manuscript.

#### `experiments/section_6_SPRINT_and_ACCORD`

Scripts for preprocessing and analyzing the SPRINT and ACCORD trials. The code generates Table 2 and Table 3 in Section 6 of the manuscript.

##### `experiments/section_6_SPRINT_and_ACCORD/preprocess_accord.R`

Script for reading in and preprocessing data from the ACCORD-BP clinical trial, specifically the `ACCORD_2017b_2` folder (available upon request from the NHLBI [website](https://biolincc.nhlbi.nih.gov/studies/accord/)). The script expects `ACCORD_2017b_2` to be stored as a subdirectory in `experiments/section_6_SPRINT_and_ACCORD/data/`. Comments are provided throughout the script to explain the preprocessing performed. Preprocessed data are stored at `experiments/section_6_SPRINT_and_ACCORD/data/accord_preprocessed.csv`.

##### `experiments/section_6_SPRINT_and_ACCORD/preprocess_sprint.R`

Script for reading in and preprocessing data from the SPRINT trial, specifically the `SPRINT_2020b` folder (available upon request from the NHLBI [website](https://biolincc.nhlbi.nih.gov/studies/sprint/)). The script expects `SPRINT_2020b` to be stored as a subdirectory in `experiments/section_6_SPRINT_and_ACCORD/data/`. Comments are provided throughout the script to explain the preprocessing performed. Preprocessed data are stored at `experiments/section_6_SPRINT_and_ACCORD/data/sprint_preprocessed.csv`.

##### `experiments/section_6_SPRINT_and_ACCORD/rate_on_SPRINT_and_ACCORD.R`

The core analysis code used to generate Table 2 and Table 3 in Section 6 of the manuscript. At a high level, the script reads in the preprocessed SPRINT and ACCORD-BP data (see above) and performs one of the following, depending on the user-provided arguments: (1) trains prioritization rules (e.g., risk estimators, CATE estimators) on SPRINT and evaluates on ACCORD-BP; (2) trains prioritization rules on ACCORD-BP and evaluates on SPRINT; (3) trains prioritization rules on a subset of the combined ACCORD-BP/SPRINT patient populations, and evaluates on a heldout test set drawn from the same combined trial cohort. The user can also specify the specific estimand to quantify causal effect (restricted mean survival time `RMST` vs. `survival.probability`). The output of this script is a file `train_on_{train_str}_test_on_{test_str}_estimating_{estimand}.csv` where `train_str` and `test_str` are both one of {`accord`, `sprint`, `combined`} as per the above (though `train_str` and `test_str` should not be the same), and `estimand` is one of {`RMST`, `survival.probability`}.

##### `experiments/section_6_SPRINT_and_ACCORD/generate_table2_and_table3.R`

This script takes in `train_on_{train_str}_test_on_{test_str}_estimating_{estimand}.csv` and generates the LaTeX used to produce Table 2 and Table 3 in the main manuscript, stored as `table2.tex` and `table3.tex`, respectively.

##### `experiments/section_6_SPRINT_and_ACCORD/run_experiment.sh`

This script takes in the `train_on_sprint_test_on_accord_estimating_RMST.csv` and `train_on_accord_test_on_sprint_estimating_RMST.csv` files generated from `rate_on_SPRINT_and_ACCORD.R` above and generates LaTeX files `table2.tex` and `table3.tex` which are used to produce tables 2 and 3 in the main manuscript, respectively.

#### `experiments/section_7_digital_marketing/`

Scripts for preprocessing and analyzing the Criteo uplift dataset.

##### `experiments/section_7_digital_marketing/criteo_analysis.R`

This script randomly selects a subset of the Criteo Uplift dataset, trains regression forest and a causal forest prioritization rules estimating potential impact of targeted advertising on website visits and conversion into a sale, then evaluates the Qini coefficient and the AUTOC on a held-out test split. The output of this script provides the results populating Table 4 in the main manuscript. RATE estimates and TOC curves are stored as `rates_and_tocs.Rdata`.

##### `experiments/section_7_digital_marketing/criteo_plot_results.R`

This script takes in the `rates_and_tocs.Rdata` output by `criteo_analysis.R` and produces the TOC curves in Figure 4 of the main manuscript. These figures are saved as `toc_comparison_visit.pdf` (for the website visit outcome) and `toc_comparison_conversion.pdf` (for the conversion to sale outcome). These plots directly compare the Qini coefficient and AUTOC formulations of the RATE in terms of their point estimates, confidence intervals, and TOC curves.

### Additional Repository Files

- `.gitignore:` Specifies files and directories to be ignored by Git.
- `LICENSE:` The repository is licensed under GPL-3.
- `requirements.txt:` Lists all the Python dependencies needed to run the code.
- `README.md:` Provides an overview of the repository structure and usage instructions.

## Data

### Section 4: Comparing RATE vs. Qini Coefficient using Synthetic Data

#### Data Availability

All of the data required for the simulation experiments are generated as part
of running the associated scripts (i.e., `qini_vs_autoc.py` via `run_experiment.sh`).

#### Data Description

We construct simulations representing:

1. A scenario in which almost all subjects exhibit a linearly varying CATE,
2. A scenario in which only a small portion of the population experiences a varying CATE, and
3. A scenario in between (1) and (2).

We draw $(n=400)$ samples from a standard uniform distribution $(X_i \sim \text{Unif}(0, 1))$, and generate potential outcomes according to the following model:

$$
\begin{cases}
Y_i(w) = \mu_w(X_i) + \varepsilon_{i}(w),~\text{where} \\
\mu_0(x) = 0 ~\text{and}~
\mu_1(x) = \max\left(-\frac{2}{p^2} x + \frac{2}{p}, 0\right)
\end{cases}
$$

where $(\varepsilon_{i}(w) \sim \mathcal{N}(0, 0.2))$ represents i.i.d. random noise and $p$ is a simulation-specific parameter representing the proportion of individuals for whom the CATE is non-zero. We draw the treatment assignment randomly with probability 0.5, so that $e = P(W_i = 1) = P(W_i = 0) = 0.5$.

Then, we consider the prioritization rule $S(X_i) = 1 - X_i$, which is "perfect" in the sense that, for all $x_i$ and $x_j$ in $[0, 1]$, $S(x_i) \geq S(x_j)$ implies $\tau(x_i) \geq \tau(x_j)$. When $p=1.0$, the CATE is nonzero for all subjects and varies linearly over quantiles of the prioritization rule. When $p=0.1$, only a small subset of the population has a nonzero CATE, but the treatment effect is large and changes quickly with the quantile.

For each $p$ in $\{1.0, 0.5, 0.1\}$, we simulate a dataset and calculate both the AUTOC and Qini metrics using oracle AIPW scores $\Gamma^{\ast}_i$.

### Section 5: The International Stroke Trial (Clinical, no censoring)

#### Data Availability

Data for the International Stroke trial can be accessed via the following link:  
https://static-content.springer.com/esm/art%3A10.1186%2F1745-6215-12-101/MediaObjects/13063_2010_637_MOESM1_ESM.CSV.
Data are provided as part of a public release via [The International Stroke Trial database](https://trialsjournal.biomedcentral.com/articles/10.1186/1745-6215-12-101) in _Trials_ by Sandercock et al, 2011.

#### Data Description

The International Stroke Trial (IST), originally [published](https://doi.org/10.1016/S0140-6736(97)04011-7) in *The Lancet* in 1997, is one of the largest RCTs ever conducted in acute stroke. It assessed the treatment effects of Aspirin, Heparin (an anticoagulant), both, or neither for patients with presumed acute ischaemic stroke in a factorial design on 19,435 patients across 36 countries. Primary outcomes included (1) death within 14 days of stroke onset, and (2) death or dependency at 6 months. Follow-up for the primary outcomes was 99% complete. We focus our analysis on the outcome of death or dependency at 6 months and restrict our analysis to those patients for whom this outcome was recorded. We also only analyze the effect of Aspirin, irrespective of Heparin assignment status. The exact variables used as treatment indicator, outcome, and predictor variables/model inputs (for risk-based and CATE-based treatment prioritization rules) are documented and explained in `section_5_international_stroke_trial/analysis.Rmd`.

### Section 6: The SPRINT and ACCORD-BP Trials (Clinical, right-censoring)

#### Data Description

Conducted in 2001-2005 and 2010-2013, respectively, the ACCORD-BP and SPRINT randomized controlled clinical trials both sought to identify whether intensive blood-pressure control could reduce the risk of adverse cardiovascular events (e.g., myocardial infarction, stroke, heart failure, death from cardiovascular causes) in individuals with elevated risk for cardiovascular disease. We used a small subset of variables collected in both trials to learn treatment prioritization rules and estimate their respective RATEs on trial data (see `experiments/section_6_SPRINT_and_ACCORD/preprocess_sprint.R` and `experiments/section_6_SPRINT_and_ACCORD/preprocess_accord.R` for details). Given the significant censoring rates in both trials, we accounted for right-censoring in both outcome and treatment effect estimation.

#### Data Availability

Data for the SPRINT and ACCORD-BP clinical trials are available to researchers free of charge via the NHLBI Biologic Specimen and Data Repository Information Coordinating Center (BioLINCC, https://biolincc.nhlbi.nih.gov/home/). Researchers must submit an online data request form, including a study plan/protocol and Institutional Review Board (IRB) approval, to obtain access.

We accessed the SPRINT and ACCORD-BP data under an NHLBI Research Materials Distribution Agreement (RMDA), signed August 20th, 2019. The document is available upon request. Further information regarding data availability, access, and metadata for the SPRINT and ACCORD-BP clinical trials can be found at the trials’ NHLBI BioLINCC portal pages (https://biolincc.nhlbi.nih.gov/studies/sprint/ and https://biolincc.nhlbi.nih.gov/studies/accord/, respectively). The trial data used in our analysis were up to date as of January 1, 2020.

### Section 7: The Criteo Uplift Dataset (Marketing)

#### Data Availability

Data for the Criteo Uplift dataset can be found at https://ailab.criteo.com/criteo-uplift-prediction-dataset/.

#### Data Description

According to the Criteo website (https://ailab.criteo.com/criteo-uplift-prediction-dataset/):
```
This dataset is constructed by assembling data resulting from several incrementality tests, a particular randomized trial procedure where a random part of the population is prevented from being targeted by advertising. It consists of 25M rows, each one representing a user with 11 features, a treatment indicator and 2 labels (visits and conversions).
```

_Privacy_  
For privacy reasons the data has been sub-sampled non-uniformly so that the original incrementality level cannot be deduced from the dataset while preserving a realistic, challenging benchmark. Feature names have been anonymized and their values randomly projected so as to keep predictive power while making it practically impossible to recover the original features or user context.

_Fields_  
Here is a detailed description of the fields (they are comma-separated in the file):  
- `f0`, `f1`, `f2`, `f3`, `f4`, `f5`, `f6`, `f7`, `f8`, `f9`, `f10`, `f11`: feature values (dense, float)  
- `treatment`: treatment group (1 = treated, 0 = control)  
- `conversion`: whether a conversion occured for this user (binary, label)  
- `visit`: whether a visit occured for this user (binary, label)  
- `exposure`: treatment effect, whether the user has been effectively exposed (binary)

## Notes

Contact the corresponding author for any additional information regarding data access and reproducibility.
