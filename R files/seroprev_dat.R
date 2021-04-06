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
neth_96 <- temporal %>% subset(year==1996)  #1996 (1995 has demographic data available)

#distribute n according to demographic data from specific year
f <- spline(age_pop, tot_pop, xout=neth_96$age_mid)
neth_96$n <- round(f$y/sum(f$y)*neth_96$n, 0)
#calculate k
neth_96$k <- round(neth_96$prevalence * neth_96$n, 0)
#recalculate prevalence
neth_96$prevalence <- neth_96$k/neth_96$n

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
neth_05 <- temporal %>% subset(year==2005)  #2005

#distribute n according to demographic data from specific year
f <- spline(age_pop, tot_pop, xout=neth_05$age_mid)
neth_05$n <- round((f$y/sum(f$y))*neth_05$n, 0)
#calculate k
neth_05$k <- round(neth_05$prevalence * neth_05$n, 0)
#recalculate prevalence
neth_05$prevalence <- neth_05$k/neth_05$n

# plot(neth_96$age_mid, neth_96$prevalence, ylim=c(0,1))
# points(neth_05$age_mid, neth_05$prevalence, col="red")
