###################################################################
## Set model parameters, demographic data, initial values & time ##
###################################################################

## Parameters to be fit to data ##
pars  <- list(log.lambda0 = log(0.02),  #baseline force of infetion (FoI) 
              log.beta = log(0.60),     #proportional change in FoI
              log.tau = log(20))        #no. years before 1st datapoint that FoI decline begins

## Toggle parameters ##
pars$post_pred     <- 1           #posterior predictions (1) or not? (0)
pars$grps_per_year <- 4           #no. age groups per year of life (options: 4 or 12/9)
pars$temporal_foi  <- "linear"    #TEMPORAL FoI decline (options: "none", "stepwise" or "linear" )
pars$age_foi       <- "constant"  #AGE-related FoI change (options: "constant", "half" or "double")
pars$troubleshoot  <- 0           #model trouble shooting parameter (options: 1 (plots graph of age-foi profile); 0 (model runs without age-foi graph plot))

# For forecasting: if pars$forecast<-1, time interval set to be longer by pars$years_forecast
if (pars$post_pred == 0) {
  pars$forecast <- 0
  
} else if (pars$post_pred == 1) {
  pars$forecast <- 1
}

# No. of years to forecast after the final data sampling year
pars$years_forecast <- 30

## Set age-group-related & other parameters ##
if (pars$grps_per_year == 4) {        #[for 3-month age groups]
  amax        <- 55                   #maximum age (in years)
  pars$mctr   <- c(0.15, 0.44, 0.71)  #mother-child transmission rate. SYROCOT 2007. The Lancet 369
  pars$burnin <- 2000                 #model burn-in period
  
} else if (pars$grps_per_year == 12/9) {    #[for 9-month age groups]
  amax        <- 60                         #maximum age (in years)
  pars$mctr   <- mean(c(0.15, 0.44, 0.71))  #mother-child transmission rate. SYROCOT 2007. The Lancet 369
  # pars$burnin <- XX                       #model burn-in period (needs to be set by user)
}

pars$amax     <- amax                          #set maximum age of population
pars$agrps    <- pars$amax*pars$grps_per_year  #no. age groups
pars$la       <- (pars$amax/pars$agrps)        #no. years in each age group
pars$da       <-  1/pars$la                    #ageing rate
pars$age      <- seq(0+pars$la/2, pars$amax-pars$la/2, length.out=pars$agrps)  #age at midpoints of age groups
pars$r        <- 1/(21/365)  #maternally-derived IgG half life of 3 weeks (Villard et al., 2016. Diagnostic microbiology and infectious disease;84(1):22-33)

if (pars$post_pred == 0) {
  pars$burnin   <- 750
  
} else if (pars$post_pred == 1) {
  pars$burnin   <- 2000
}

## Seroprevalence data ##
if (cluster == "none") { #local
  df <- readRDS("data/global_data.rds")
  
} else if (cluster == "RVC" | cluster == "UCL") { #cluster
  df <- readRDS("global_data.rds")
}

countries    <- sort(unique(df$country))                #countries in the dataset
pars$country <- countries[1]                            #set country
fitting_data <- subset(df, df$country == pars$country)  #subset data

## Demographic data ##
if (cluster == "none"){ #local
  source("R files/demogdat.R")
  
} else if (cluster == "RVC" | cluster == "UCL"){ #cluster
  source("demogdat.R")
}

# No. of years between 1st & last data time points
pars$tdiff <- max(fitting_data$year) - min(fitting_data$year)

# No. of years between each data point from first to last
year_diff <- vector("numeric", length = nrow(fitting_data))
for (i in 1:nrow(fitting_data)-1) {
  year_diff[i] <- fitting_data$year[i+1] - fitting_data$year[i]
}

year_diff <- c(pars$burnin, year_diff[1:(length(year_diff)-1)])  #add burn-in time
year_diff <- cumsum(year_diff)  #cumulative sum to get model sampling times

## Set interpolated parameters from demographic data (see demogdat.R)
# Population size
f       <- spline(age_pop, tot_pop, xout=pars$age)
pars$Na <- f$y*(sum(tot_pop[1:which(age_pop==pars$amax-2.5)])/sum(f$y))  #select correct age specific population size depending on choice of pars$amax
pars$N  <- sum(pars$Na)  #total population size

# Mortality rate
pars$d <- spline(x=age_mort, y=mort, xout=pars$age)$y  #interpolated age specific per capita mortality rates
pars$d[which(pars$d<0)] <- 0   #zero minus elements

# Birth rate
pars$propfert <- approx(age_fert, prop_fert_age, xout=pars$age)$y
pars$propfert[is.na(pars$propfert)] <- 0
pars$propfert <- pars$propfert/sum(pars$propfert)

## Initial values (state variables) ##
S0 <- I0 <- Im0 <- vector("numeric", length=pars$agrps)

# Entire population initially susceptible
S0[1:length(S0)] <- pars$Na
y                <- c(S0, I0, Im0)

## Time ##
if (pars$forecast == 0) {
  time <- seq(1, pars$burnin + pars$tdiff, 1)
  
} else if (pars$forecast == 1) {
  time <- seq(1, pars$burnin + pars$tdiff + pars$years_forecast, 1)
}
