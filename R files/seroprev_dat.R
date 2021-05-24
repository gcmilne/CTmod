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

# #################################################
# ################## Netherlands ##################
# #################################################
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
# neth_95 <- temporal %>% subset(year==1995)  #1995/96 (1995 has demographic data available)
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
# neth_06 <- temporal %>% subset(year==2006)  #2006/07 (2005 has demographic data available)
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
# ##### calculate age indices that match the model ouput for "fitting.R"
# data <- neth_95
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
# neth_matched_indices <- which(!is.na(matched_ages))
# 
# #################################################
# ################## New Zealand ##################
# #################################################
# ## NB make sure demogdat.R is set to get data from NZ
# 
# #read in temporal data from New Zealand
# nz1 <- data.frame(n=566, k=round(0.6*566), prev=0.60, year=1982)    # Cursons et al 1982
# nz2 <- data.frame(n=500, k=round(0.354*500), prev=0.354, year=2004) # Morris et al 2004
# 
# # assume representative sampling of pregnant woman w.r.t. age
# preg <- subset(pars$propfert, pars$propfert>0)        #keep only pregnant women
# min.age <- which(pars$propfert == preg[1])            #save 1st age index of pregnant women
# max.age <- which(pars$propfert == preg[length(preg)]) #save last age index of pregnant women
# nz_matched_indices <- seq(min.age, max.age, by=1)     #create sequence of age indices for extracting correct modelled seroprevalence
# 

################################################################
################## Global seroprevalence data ##################
################################################################
library(dplyr)
library(ggplot2)
library(patchwork)

data <- read.csv("data/seroprev_global.csv")
min_samples <- 4 #minimum number of longitudinal samples

# clean data
data <- data %>%
  filter(exclude=="n") %>%                            # exclude rows with epidemiological biases or other errors
  filter(!is.na(first_sample)) #%>%                    # remove rows with no info on sampling year
  # filter(method=="ELISA" & method2=="" & method3=="") # only include studies using ELISA (& no other methods)
  # filter(method=="ELISA" & method2=="" & method3=="") # only include studies using ELISA (in combo w other methods)

# calculate prevalence
data$prev <- data$k/data$n

# calculate median sampling year
for(i in 1:nrow(data)){
  data$year[i] <- floor(median(seq(data$first_sample[i],data$last_sample[i], 1)))
}

# include countries with minimum of x observations
data <- data %>%
  group_by(country) %>%
  filter(n() >= min_samples)

# create separate dfs for countries, group by year & calculate weighted prevalence
countries <- unique(data$country)
data_list <- vector(mode = "list", length = length(countries))

for(i in 1:length(countries)){
  data_list[[i]] <- 
    data %>%
    filter(country==countries[i]) %>% 
    filter(!is.na(year)) %>%
    arrange(year) %>%
    group_by(year) %>%
    summarise(w.prev = weighted.mean(prev, n), k=sum(k), n = sum(n), country=countries[i], state=state) %>%
    distinct(k, n, .keep_all = T) %>% #remove duplicate rows (where median sampling years are equal)
    as.data.frame()
}

# make ggplots for each of the countries
p <- vector(mode = "list", length = length(countries))

#for(j in 1:length(who_regions)){  ## if you want to order plots by region too
for(i in 1:length(countries)){
  p[[i]] <- ggplot(data_list[[i]], aes(x=year, y=w.prev)) + 
    geom_point() + 
    geom_smooth(method = "lm", formula = 'y~x') +
    xlab("Year") + ylab("Seroprevalence") +
    scale_y_continuous(breaks = seq(0,1, 0.2)) +
    labs(title = paste(countries[i])) + 
    theme_light()
}

# wrap_plots(p)

####################
## Regional plots ##
####################

# omit rows with no regional data
for(i in 1:length(countries)){
  data_list[[i]] <- data_list[[i]] %>%
    filter(!is.na(year)) %>%
    filter(!is.na(state))
}

# store names of unique states
states <- vector(mode = "list", length = length(countries))
for(i in 1:length(countries)){
  states[[i]] <- unique(data_list[[i]]$state)
}

# find n times each unique state is present in dataset
n_states <- table(data$state)

# store names of states if >n number of samples
n <- 3

freq_states <- vector("numeric", length=length(n_states))

for(i in 1:length(n_states)){
  if(n_states[i] >= n){
    print(n_states[i])
    freq_states[i] <- names(n_states[i])
  }else{
    freq_states[i] <- NA
  }
}

freq_states <- freq_states[!is.na(freq_states)] #omit na

## store data & calculate weighted prevalence for these states
regional_data <- vector(mode = "list", length = length(freq_states))

for(i in 1:length(freq_states)){
  regional_data[[i]] <- data %>%
    filter(!is.na(year)) %>%
    filter(!is.na(state)) %>%
    filter(state==freq_states[[i]]) %>%
    arrange(year) %>%
    group_by(year) %>%
    summarise(w.prev = weighted.mean(prev, n), k=sum(k), n = sum(n), country=country, state=freq_states[[i]]) %>%
    as.data.frame()
}

# omit those with fewer than n observations
for(i in 1:length(freq_states)){
  if(length(unique(regional_data[[i]]$year)) < n){
    regional_data[[i]] <- NA
  }
}

regional_data <- regional_data[!is.na(regional_data)] #omit na

# ggplot for each state
p <- vector(mode = "list", length = length(regional_data))
for(i in 1:length(regional_data)){
  p[[i]] <- vector(mode = "list", length = length(regional_data[[i]]))
}

for(i in 1:length(regional_data)){
  p[[i]] <- ggplot(regional_data[[i]], aes(x=year, y=w.prev)) + 
    geom_point() + 
    geom_smooth(method = "lm", formula = 'y~x') +
    scale_y_continuous(breaks = seq(0,1, 0.2)) +
    labs(title = paste(c(freq_states[i]), regional_data[[i]]$country, sep=", ")) + 
         xlab("Year") +
         ylab("Seroprevalence")
}

# wrap_plots(p)


############
## Linear models
############

## Option (1): Fit linear models for each country separately
lm.mod   <- vector("list", length=length(countries))  #unweighted model
lm.mod_w <- vector("list", length=length(countries))  #weighted model


# save lm model output
for(i in 1:length(lm.mod)){
  #unweighted model
  lm.mod[[i]]           <- lm(w.prev ~ year, data=data_list[[i]])
  lm.mod[[i]]$summary   <- summary(lm.mod[[i]])
  #weighted model
  lm.mod_w[[i]]         <- lm(w.prev ~ year, data=data_list[[i]], weights = n)
  lm.mod_w[[i]]$summary <- summary(lm.mod_w[[i]])
}

# save p values for year term
p_value   <- vector("list", length=length(countries))  #unweighted model
p_value_w <- vector("list", length=length(countries))  #weighted model

for(i in 1:length(p_value)){
  p_value[[i]]   <- lm.mod[[i]]$summary$coefficients["year","Pr(>|t|)"]
  p_value_w[[i]] <- lm.mod_w[[i]]$summary$coefficients["year","Pr(>|t|)"]
  
}

# Probe the significance of temporally changing seroprevalence, for each country separately
## unweighted model
for(i in 1:length(p_value)){
  if(p_value[[i]] < 0.05){#find the P value of the year coefficient
    print(paste("Significant: ", countries[i]))
  } else if (p_value[[i]] >=0.05 & p_value[[i]] < 0.10){
    print(paste("Marginally significant: ", countries[i]))
  } else if (p_value[[i]] >=0.10){
    print(paste("NS: ", countries[i]))
  }
}

## weighted model
for(i in 1:length(p_value_w)){
  if(p_value_w[[i]] < 0.05){#find the P value of the year coefficient
    print(paste("Significant: ", countries[i]))
  } else if (p_value_w[[i]] >=0.05 & p_value_w[[i]] < 0.10){
    print(paste("Marginally significant: ", countries[i]))
  } else if (p_value_w[[i]] >=0.10){
    print(paste("NS: ", countries[i]))
  }
}

## Option (2): Fit a hierarchical model for all countries together

## Fixed effects model ##

# put data list into one dataframe
df <- do.call(rbind.data.frame, data_list)

#unweighted model
lm.mod <- lm(w.prev ~ year, data=df)
summary(lm.mod)

#weighted model
lm.mod_w <- lm(w.prev ~ year, data=df, weights=n)
# summary(lm.mod_w)

# plot fitted points
# plot(df$year, df$w.prev)
# lines(df$year, lm.mod$fitted.values)              #unweighted model
# lines(df$year, lm.mod_w$fitted.values, col="red") #weighted model

## Mixed effects model (with country as random effect) ##
library(nlme)

df$country <- as.factor(df$country)

# weighted mixed effects model (with country as random intercept and random gradient)
me_mod <- lme(fixed=w.prev ~ year, random=~year|country, data=df, weights=~n)
# summary(me_mod)

# plot the output of the me_mod vs the observed data
par(mfrow=c(1,1))
# plot(df$year, df$w.prev, col=df$country, ylim=c(0, 1)) #observed data
# for(i in 1:nrow(country_fits)){ #fitted line for each country 
#   abline(me_mod$coefficients$fixed[1] + me_mod$coefficients$random$country[i,1], 
#          me_mod$coefficients$fixed[2] + me_mod$coefficients$random$country[i,2], col=df$country, 
#          lty=i)
# }
# abline(me_mod$coefficients$fixed[1], me_mod$coefficients$fixed[2], col="red") #best-fit line (all countries)

## Examine residual plots ##
par(mfrow=c(2,1))
# Plot residuals vs. fitted values
# plot(me_mod$fitted[,1], me_mod$residuals[,1], xlab="Fitted values", ylab="Residuals")
# lines(seq(0.2,0.7,by=0.1), rep(0,6), lty=3)
# Plot residuals vs. predictor (year)
# plot(df$year, me_mod$residuals[,1], xlab="Year", ylab="Residuals")
# lines(seq(1980, 2020, by=10), rep(0,5), lty=3)
