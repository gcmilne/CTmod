CTmod
===================

# Description of model #

A mathematical model for estimating the time-varying seroprevalence of anti-Toxoplasma gondii IgG and the incidence of congenital toxoplasmosis (CT). The model is run using age- and country-specific data on demographic factors (fertility distributions, mortality distributions, population size), and immunoassay-specific information on diagnostic sensitivity and specificity. To generate past, current, and future predictions of temporal changes in seroprevalence and CT incidence, the model can be fit (using an Approximate Bayesian framework) to representative country-specific IgG seroprevalence data in pregnancy.

# Description of R scripts #

### (1) seroprev_dat.R ###

A script to load in & clean the .csv file seroprevalence data (remove NAs, apply inclusion & exclusion criteria). The minimum number of timepoints per country for inclusion in further analysis can be specified in this script. The output of the script is an RDS file containing all of the included seroprevalence data as a dataframe. This script isn't required once it has been run once (unless the user wants to change the specified data).

### (2) diagnostics.R ###

A script to load in data on the sensitivity and specificity of various anti-Toxoplasma gondii IgG antibody immunoassays. The output of this script is a dataframe containing all of the information on each diagnostic, which is loaded into subsequent scripts at the relevant stages through source().

### (3) setparms.R ###

For setting the model parameters, loading in demographic data (see demogdat.R), setting initial values & time. Also includes useful toggle parameters to control the type of age- and time-related force of infection (FoI), and other useful parameters to control posterior predictions. Note that the model can work in 3-month or 9-month age intervals (the former was used in all published analyses and has hence the associated code has been more thoroughly tested; note the latter may not work off-the-shelf).

### (4) demogdat.R ###

For getting demographic data for the given country (specified in setparms.R as pars$country), for a given year (automated in demogdat.R to be the first time point of the data (fitting_data, set in setparms.R), rounded to the nearest 5 years).

### (5) linear_mods.R ###

Linear regressions of the relationship between seroprevalence and sampling year, including a plot. Also provides an estimate of the average annual change in seroprevalence.

### (6) model.R ###

The mathematical model, which can be configured with various options by parameters in the setparms.R script. Notably, the age groupings can be changed, as can the form of the FoI over age and over time.

### (7) funcs.R ###

Contains all the functions needed for fitting the model to the data, and for posterior predictions. These functions are source()'d in the relevant scripts.

### (8) fitting.R ###

This script is used to fit the model (model.R) to the data (fitting_data, loaded in setparms.R), using the parameters and demographic data. Options are available to perform the fitting on a local machine (cluster \<- "none"), or on a cluster directory (which can be set at the beginning of this script (current options include cluster \<- "RVC" and cluster \<- "UCL"). Outputs of the fitting procedure are saved as RDS files in the set working directory, labelled by the country being fit to, and the chosen age- and time-specific FoI model.

### (9) gen_postdist.R ###

After fitting the model to the data, this script can be used to load in the .RDS files generated during the fitting procedure, summarise and plot the posterior distributions (including plots of priors vs. posteriors for all countries in the dataset, and dot and whisker plots of the medians + uncertainties of the posterior distributions), and save the posterior distributions for use by subsequent scripts.

### (10) post_pred.R ###

This script uses the posterior distribution generated in the gen_postdist.R script to perform posterior predictions on the seroprevalence of T. gondii exposure and incidence of congenital toxoplasmosis. The cluster object can be used to set the working directory depending on whether the script is being run on a local machine (cluster \<- "none") or cluster (current options: cluster \<- "RVC", or cluster \<- "UCL"). The outputs of this script are multiple .RDS files containing the posterior predictions and labelled by the chosen country and age- and time-related FoI model.

### (11) save_postpred.R ###

This script reads in the saved .RDS files generated following posterior predictions (performed in post_pred.R), transforms lists into easier-to-use dataframes, calculates uncertainty in model posterior predictions, adjusts seroprevalence data to a "true" seroprevalence, using information on diagnostic performance (diagnostics.R), that be compared to modelled seroprevalence estimates. The ouput of the script is 3 datasets: (1) A dataset containing seroprevalence data & model estimates for the years with data; (2) A dataset containing model seroprevalence estimates for a larger range of years, not limited to those years with data; (3) A dataset containing model CT incidence estimates (transformed to be per 10,000 the country-specific number of live births). This final dataset also includes information on the proportion of CT cases leading to particular sequelae (obtained from the literature).

### (12) plot_postpred.R ###

Reads in posterior predictions saved in save_postpred.R and plots modelled estimates of seroprevalence and CT incidence (as well as data where relevant). The outputs of this script are PDF figures, or high-res PNG figures.

### (13) temporal_ct.R ###

Uses saved output of posterior predictions (from save_postpred.R) to identify which countries have had increasing CT incidence relative to some user-defined baseline year. This script is also used to compare the model's estimates to those from a previous study, by extracting the model's country-specific posterior predictions from specific years corresponding to the country-specific years of data collection in the other study.

### (14) suppl_ggplots.R ###

This script is for plotting figures that are supplementary to the main posterior predictions analyses (e.g. plots of demographic data for each country, force of infection predictions, etc.).
