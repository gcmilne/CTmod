# stan-model
Authors: Gregory Milne, Martin Walker
Affiliation: Royal Veterinary College, Hawkshead campus, U.K.

A deterministic model that tracks the change in the number of susceptible, infected individuals and those with maternal antibodies to Toxoplasma gondii (S-I-M). 

The model aims to predict the number of cases of congenital toxoplasmosis (CT) given:
- (a) Country-level demographic data on population structure over age, age-specific fertility and mortality data
- (b) Country-level representative IgG seroprevalence data (which is used to predict force of infection and no. seroconversions in pregnancy)

The model is written in Stan (stan-mod.stan) and compiled in R (stan-mod.R). The model uses parameters (setparms.R) and demographic data (demogdat.R). 

The model's initial version in R can be viewed in the model.R file. For this version there is a corresponding runmod.R file for model running.
