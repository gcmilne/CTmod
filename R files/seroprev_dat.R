#####
# Script to calculate age-specific k and n for each seroprevalence dataset 
# (given only total N & age-specific seroprevalence)
# Assumes sampling is proportional to age distribution of country in given sampling year
#####

## NB make sure demogdat.R is set to get data from the correct country & year

# Load scripts #
source("R files/demogdat.R")  #local
# source("demogdat.R")  #cluster

# Load packages
library(dplyr)
library(wpp2019)

#################################################
################## Netherlands ##################
#################################################
# temporal <- read.csv("data/netherlands_temporal.csv")  #local
# # temporal <- read.csv("netherlands_temporal.csv")  #cluster
# 
# #load required datasets
# data(popF) # population distribution females
# data(popM) # population distribution males
# 
# ### set country & year ###
# country <- "Netherlands"
# year1 <- "1995"
# 
# #calculate total pop size for given year
# tmpdf <- subset(popF, name==country)
# pop_f <- tmpdf[,year1]*1000 # total population by age (females)
# tmpdf <- subset(popM, name==country)
# pop_m <- tmpdf[,year1]*1000 # total population by age (males)
# tot_pop <- pop_m + pop_f #total population all ages
# age_pop <- age_mort[-1] #set age categories
# age_pop[1] <- age_pop[1]-1
# tot_pop <- tot_pop[-length(tot_pop)]
# 
# #subset Netherlands data by year of data collection
# neth_95 <- temporal %>% subset(year==1996)  #1995/96 (1995 has demographic data available)
# 
# #distribute n according to demographic data from specific year
# f <- spline(age_pop, tot_pop, xout=neth_95$age_mid)
# neth_95$n <- round(f$y/sum(f$y)*neth_95$n, 0)
# #calculate k
# neth_95$k <- round(neth_95$prevalence * neth_95$n, 0)
# #recalculate prevalence
# neth_95$prevalence <- neth_95$k/neth_95$n
# 
# ### set country & year ###
# country <- "Netherlands"
# year1 <- "2005"
# 
# #calculate total pop size for given year
# tmpdf <- subset(popF, name==country)
# pop_f <- tmpdf[,year1]*1000 # total population by age (females)
# tmpdf <- subset(popM, name==country)
# pop_m <- tmpdf[,year1]*1000 # total population by age (males)
# tot_pop <- pop_m + pop_f #total population all ages
# age_pop <- age_mort[-1] #set age categories
# age_pop[1] <- age_pop[1]-1
# tot_pop <- tot_pop[-length(tot_pop)]
# 
# #subset Netherlands data by year of data collection
# neth_06 <- temporal %>% subset(year==2005)  #2006/07 (2005 has demographic data available)
# 
# #distribute n according to demographic data from specific year
# f <- spline(age_pop, tot_pop, xout=neth_06$age_mid)
# neth_06$n <- round((f$y/sum(f$y))*neth_06$n, 0)
# #calculate k
# neth_06$k <- round(neth_06$prevalence * neth_06$n, 0)
# #recalculate prevalence
# neth_06$prevalence <- neth_06$k/neth_06$n
# 
# # plot(neth_95$age_mid, neth_95$prevalence, ylim=c(0,1))
# # points(neth_06$age_mid, neth_06$prevalence, col="red")
# 
# ##### set up to run model in "fitting.R"
# ##### Read in data #
# data <- neth_95  # from "R files/seroprev_dat.R" script (NB length of datasets from both years are the same)
# number_of_data_points = length(data$n)
# 
# ### select age groups from model output that match data age groups
# clean_dat <- data.frame("age_mid"=data$age_mid, "k"=data$k, "n"=data$n, "prev"=data$prevalence)
# 
# #create new dataset
# matched_dat <- clean_dat
# # x[,2] increases by one each time data age midpoint is closest match to modelled age midpoint
# x <- cbind(pars$age, findInterval(pars$age, matched_dat$age_mid))
# #returns FALSE if there's change between element i and element i+1
# y1 <- diff(x[,2]) <= 0
# 
# # each time x[,2] increases by 1, save value of x[,1][i+1] to matched_dat$age_mid[i]
# matched_ages <- vector("numeric", pars$agrps)
# for(i in 1:length(y1)){
#   if(y1[i]==T){
#     matched_ages[i] <- NA
#   }else if(y1[i]==F){
#     matched_ages[i] <- x[,1][i+1]
#   }
# }
# 
# #removes last element (which is 0 because of indexing)
# matched_ages <- head(matched_ages, -1)
# 
# # find the indices which match the age group most closely
# matched_indices <- which(!is.na(matched_ages))

#################################################
################## New Zealand ##################
#################################################
## NB make sure demogdat.R is set to get data from NZ

#read in temporal data from New Zealand
nz1 <- data.frame(n=566, k=round(0.6*566), prev=0.60, year=1982)    # Cursons et al 1982
nz2 <- data.frame(n=500, k=round(0.354*500), prev=0.354, year=2004) # Morris et al 2004

# assume representative sampling of pregnant woman w.r.t. age
preg <- subset(pars$propfert, pars$propfert>0)        #keep only pregnant women
min.age <- which(pars$propfert == preg[1])            #save 1st age index of pregnant women
max.age <- which(pars$propfert == preg[length(preg)]) #save last age index of pregnant women
nz_matched_indices <- seq(min.age, max.age, by=1)     #create sequence of age indices for extracting correct modelled seroprevalence

