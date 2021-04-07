#####
# Script to calculate age-specific k and n for each seroprevalence dataset 
# (given only total N & age-specific seroprevalence)
# Assumes sampling is proportional to age distribution of country in given sampling year
#####

# Load scripts #
source("R files/demogdat.R")

# Load packages
library(dplyr)
library(wpp2019)

#read in temporal data from the Netherlands
temporal <- read.csv("data/netherlands_temporal.csv")

#load required datasets
data(popF) # population distribution females
data(popM) # population distribution males

### set country & year ###
country <- "Netherlands"
year1 <- "1995"

#calculate total pop size for given year
tmpdf <- subset(popF, name==country)
pop_f <- tmpdf[,year1]*1000 # total population by age (females)
tmpdf <- subset(popM, name==country)
pop_m <- tmpdf[,year1]*1000 # total population by age (males)
tot_pop <- pop_m + pop_f #total population all ages
age_pop <- age_mort[-1] #set age categories
age_pop[1] <- age_pop[1]-1
tot_pop <- tot_pop[-length(tot_pop)]

#subset Netherlands data by year of data collection
neth_95 <- temporal %>% subset(year==1996)  #1995/96 (1995 has demographic data available)

#distribute n according to demographic data from specific year
f <- spline(age_pop, tot_pop, xout=neth_95$age_mid)
neth_95$n <- round(f$y/sum(f$y)*neth_95$n, 0)
#calculate k
neth_95$k <- round(neth_95$prevalence * neth_95$n, 0)
#recalculate prevalence
neth_95$prevalence <- neth_95$k/neth_95$n

### set country & year ###
country <- "Netherlands"
year1 <- "2005"

#calculate total pop size for given year
tmpdf <- subset(popF, name==country)
pop_f <- tmpdf[,year1]*1000 # total population by age (females)
tmpdf <- subset(popM, name==country)
pop_m <- tmpdf[,year1]*1000 # total population by age (males)
tot_pop <- pop_m + pop_f #total population all ages
age_pop <- age_mort[-1] #set age categories
age_pop[1] <- age_pop[1]-1
tot_pop <- tot_pop[-length(tot_pop)]

#subset Netherlands data by year of data collection
neth_06 <- temporal %>% subset(year==2005)  #2006/07 (2005 has demographic data available)

#distribute n according to demographic data from specific year
f <- spline(age_pop, tot_pop, xout=neth_06$age_mid)
neth_06$n <- round((f$y/sum(f$y))*neth_06$n, 0)
#calculate k
neth_06$k <- round(neth_06$prevalence * neth_06$n, 0)
#recalculate prevalence
neth_06$prevalence <- neth_06$k/neth_06$n

# plot(neth_95$age_mid, neth_95$prevalence, ylim=c(0,1))
# points(neth_06$age_mid, neth_06$prevalence, col="red")
