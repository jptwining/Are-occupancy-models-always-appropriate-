# Are-occupancy-models-always-appropriate-
Simulations examining fitting occupancy models compared to other analytical options including logistic regression, Poisson and NB GLMs, as well as variants of N-mixture models to detection/non-detection and count data across different scenarios with heterogeneity in detection probability and varying sample sizes.

# Summary
Contained within this repository are:
1. Scripts for simulation study 1.1. Simulations exploring modelling choices for count data where naïve occupancy ~ 1
1. Scripts for simulation study 1.2.	Simulations exploring use of logistic regression vs. occupancy model when p <1 and there is heterogeneity in detection probability.

## The working directory

Below you will find descriptions of each folder in this repository and files contained within them.

## The folder directory (./main/...)

This folder has two subfolders

## 1. Simulation Study 1.1 (./main/Simulation1.1)

### 1.1. Simulation_study_modelling_count_data_examining_model_choices_final.R

This script provides functions and code to run a simulation study that explores impact of fitting an occupancy model vs. a range of different models (Poisson GLM, NB GLM, Poisson-Poisson N-mixture model, binomial-Poisson N-mixture) to simulated count data

## 2.  Simulation Study 1.2 (./main/Simulation1.2)

### 2.1. Simulation_study_occupancy_vs_logistic_I=200_final.R

This script provides functions and code to run a simulation study that explores impact of fitting an occupancy model vs. a Bernoulli GLM to detection/non-detection data form sampling at 200 sites

### 2.2. Simulation_study_occupancy_vs_logistic_I=50_final.R
This script provides functions and code to run a simulation study that explores impact of fitting an occupancy model vs. a Bernoulli GLM to detection/non-detection data form sampling at 50 sites
