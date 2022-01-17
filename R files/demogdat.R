######################################
## Country-specific demographic data ##
#######################################

### set country & year ###
country <- pars$country
# 5-year band
year <- paste(as.character(round(min(fitting_data$year)/5)*5), 
      as.character((round(min(fitting_data$year)/5)*5)+5), sep= "-")  
# single year (first timepoint)
year1 <- as.character(round(min(fitting_data$year)/5)*5)

library(wpp2019)

### load required datasets
data(tfr) # total fertility rate (births per woman over lifetime)
data(percentASFR) # age distribution of fertility rates (percentage)
data(popF) # population distribution females
data(popM) # population distribution males
data(mxF) # age-spcific mortality rate for females
data(mxM) # age-specific mortality rate for males

## 1. mortality 
tmpdf <- subset(mxF, name==country)
mort_f <- tmpdf[,year]

tmpdf <- subset(mxM, name==country)
mort_m <- tmpdf[,year]

## average mortality rate (male and female)
mort <- apply(cbind(mort_m, mort_f), 1, mean)
age_mort <- tmpdf[,"age"]
age_mort <- age_mort[-length(age_mort)]+diff(age_mort)/2
mort <- mort[-length(mort)]
rm(mort_f, mort_m)

## 2. total population size
tmpdf <- subset(popF, name==country)
pop_f <- tmpdf[,year1]*1000 # total population by age (females)

tmpdf <- subset(popM, name==country)
pop_m <- tmpdf[,year1]*1000 # total population by age (males)

tot_pop <- pop_m + pop_f #total population all ages

age_pop <- age_mort[-1]
age_pop[1] <- age_pop[1]-1
tot_pop <- tot_pop[-length(tot_pop)]

## 3. age-specific births rates
age_pop_cat <- tmpdf[,"age"]
tmpdf <- subset(percentASFR, name==country)
age_fert_cat <- tmpdf[,"age"]
index <- !is.na(match(age_pop_cat,age_fert_cat))
index <- index[-length(index)]

prop_fert_age <- rep(0,length=length(index))
prop_fert_age[index] <- tmpdf[,year]/100
age_fert <- age_pop
rm(index, age_pop_cat, age_fert_cat)
