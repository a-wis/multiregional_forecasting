# Probabilistic Multiregional Population Forecasting

This file describes the use of the code for preparing probabilistic subnational population forecasts presented in the manuscript

**A unified framework for probabilistic forecasts of subnational populations**

Copyright: (anonymised), 2023

This code will be available on GitHub

All calculations were carried out in R 4.2.1 and rstan Version 2.26.18

# Contents of the folder

The repository contains four folders: 

- code (R code for the analysis), 

- data (all data for Australia used in the illustration), 

- models (all models presented in the manuscript written in Stan), 

- plots (plots created using the code). 

# How to use the code

1. File '1_data_processing.R' reads in and prepares all data required for computations. It can be sourced as is.

2. Once the data are read in, one can use frequentist approach to select parameters of the bi-log-linear model by using RMSE of the log-linear model fit or information criteria; this is available in '2_model_selection.R'.

3. File '3_component_estimation.R' contains code to call Stan models in folder 'models' to make forecasts of population components. Each of the models requires preparing data in Stan format, initial values and running the model. The last few lines contain sample code for checking convergence.

4. File '4_goodness_of_fit.R' contains several functions to visualise and check how well the models fit the data. 

5. Function for transforming the results into data frames that can be used for producing outputs are in file '5.0_functions_forecasts_transformations.R', whereas functions for multiregional projection are in file '5.1_functions_multiregional_projection.R'. The results from 3. (outputs of log-bilinear component models) are taken into '5.2_create_projection.R' and produces population forecasts.

6. File '6.1_master_create_outputs.R' reads in the results of population components models (2.) and the population forecasts (3.), transforms the results and joins with input data for each component separately, and reproduces all plots and tables in the manuscript (note: the outputs from the Bayesian bilinear models (3.) need to be in a folder 'outputs'.

7. File '6.2_create_outputs_ageprofiles.R' transforms results and data to reproduce plots of age profiles in the Online Supplement. 

# Other

- running the log-bilinear models for population components can take several hours. By using a five-year-old PC with 16 GB RAM, it takes around 11 hours to obtain the forecasts for the internal migration model. All other models take less time (anything between 2 to 10 hours).

- the raw results of running a stan model take up lots of RAM 

- we recommend using RStudio to use R files as chunks of code can be navigated by using Contents menu in the bottom of the 'Source' panel

- the code contains samples of 1000 (two chains 500 each). These yield maximum Rhat ~1.1. Longer chains reduce it but take up resources (RAM and disc space). 



