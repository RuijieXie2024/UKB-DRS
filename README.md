# Novel Type 2 Diabetes Prediction Score Based on Traditional Risk Factors and Circulating Metabolites

## Overview

This repository contains the R code used for the analysis in the manuscript titled:

**"Novel type 2 diabetes prediction score based on traditional risk factors and circulating metabolites: Model derivation and validation in two large cohort studies."**

The analysis involves developing and validating a predictive model for type 2 diabetes using traditional risk factors and circulating metabolites. The study utilizes data from two large cohort studies.

## Data Access

Due to data sharing agreements and privacy concerns, the datasets used in this analysis are not publicly available. Researchers interested in accessing similar data should apply through the appropriate channels:

- **UK Biobank (Derivation Cohort):** [Apply for Access](https://www.ukbiobank.ac.uk/enable-your-research/apply-for-access)
- **ESTHER Study (Validation Cohort):** Data access requests can be directed to the study coordinators.

## Requirements

- **R version:** 4.0 or higher
- **R Packages:**

  - `survival`
  - `glmnet`
  - `survIDINRI`
  - `nricens`
  - `pec`
  - `dplyr`
  - `tidyr`
  - `doParallel`
  - `missForest`
  - `ggplot2`
  - `gridExtra`
  - `ggcorrplot`
  - `pROC`
  - `rms`

## Installation

Install the required R packages using the following command:

```R
install.packages(c("survival", "glmnet", "survIDINRI", "nricens", "pec", "dplyr", "tidyr",
                   "doParallel", "missForest", "ggplot2", "gridExtra", "ggcorrplot",
                   "pROC", "rms"))
