##############
### ToxMod ###
##############
A mathematical model for estimating the incidence of congenital toxoplasmosis (CT) 
and seroprevalence of T. gondii exposure

## Authors ## 
Gregory Milne, Martin Walker

## Affiliation: ##
Royal Veterinary College, Hawkshead campus, U.K.

A deterministic model that tracks the change in the number of susceptible, 
infected individuals and those with maternal antibodies to Toxoplasma gondii (S-I-M). 

## Description of R scripts: ##

# (1) seroprev_dat.R 
A script to load in & clean the .csv file seroprevalence data (remove NAs, apply 
inclusion & exclusion criteria). The minimum number of timepoints per country for 
inclusion in further analysis can be specified in this script. The output of the 
script is an RDS file containing all of the included seroprevalence data as a 
data.frame. This script isn't required once it has been run once (unless the user 
wants to change the specified data).

# (2) diagnostics.R
A script to load in data on the sensitivity and specificity of various anti-
Toxoplasma gondii IgG antibody immunoassays. The output of this script is a 
dataframe containing all of the information on each diagnostic, which is loaded 
into subsequent scripts at the relevant stages.

# (3) setparms.R 
For setting the model parameters. Also includes useful toggle parameters for the 
posterior prediction stage.

# (4) demogdat.R 
For getting demographic data for the given country (specified in setparms.R as 
pars$country), for a given year (automated in demogdat.R to be the first time 
point of the data (fitting_data, set in setparms.R), rounded to the nearest 5 years.

# (5) linear_mods.R 
A short script for testing for a linear relationship between seroprevalence and 
year of data collection.

# (6) model.R 
The mathematical model, which can be configured with various options by parameters 
in the setparms.R script.

# (7) funcs.R
Contains all the functions needed for fitting the model to the data, and for 
posterior predictions. These functions are source()'d in the relevant scripts.

# (8) fitting.R 
This script is used to fit the model (model.R) to the data (fitting_data, loaded 
in setparms.R), using the parameters and demographic data. Options are available
to perform the fitting on a local machine (cluster <- "none"), on the RVC cluster
(cluster <- "RVC") or on the UCL cluster (cluster <- "UCL"). Outputs of the fitting 
procedure are saved as RDS files in the set working directory, labelled by the 
country being fit to, and the chosen age- and time-specific force of infection model

# (9) gen_postdist.R 
After fitting the model to the data, this script can be used to load in the RDS 
files generated during the fitting procedure, summarise and plot the posterior 
distribution, and save the posterior distribution for use by subsequent scripts.

# (10) post_pred.R 
This script uses the posterior distribution generated in the gen_postdist.R script 
to perform posterior predictions on the seroprevalence of Toxoplasma gondii 
exposure and incidence of congenital toxoplasmosis. The cluster object can be 
used to set the working directory depending on whether the script is being run 
on a local machine (cluster <- "none") or cluster (current options: 
cluster <- "RVC", or cluster <- "UCL"). The output of this script is multiple RDS 
files containing the posterior predictions.

# (11) plot_postpred.R 
After performing the posterior predictions using the post.pred.R script, this 
script reads in the generated RDS file output, calculates median, upper and 
lower intervals for the predictions and plots these against the data & 95% CIs 
(data adjusted to true prevalence to match the modelled output). The plots can 
be arranged such that predictions of one kind are all shown together, for all of 
the countries (e.g. seroprevalence over the years for which data were available 
for each given country), or such that all predictions are shown for an individual
country. The output of this script are saved PDF files containing ggplots of the 
modelled posterior predictions.

# (12) suppl_ggplots.R 
This script is for plotting supplementary figures (e.g. plots of demographic data 
for each country, model burn-in period, etc.).
