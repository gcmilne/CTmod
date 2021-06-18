#####################################
## SET PARAMETERS & INITIAL VALUES ##
#####################################

pars  <- list(amax=55, 
              log.lambda0 = log(0.02),
              log.shape= log(0.60),
              log.tdecline = log(20))

pars$agrps    <- pars$amax*4
pars$la       <- (pars$amax/pars$agrps)  # no. years in each age group
pars$da       <-  1/pars$la  # ageing rate
pars$age      <- seq(0+pars$la/2, pars$amax-pars$la/2, length.out=pars$agrps) # age at midpoints of age groups
pars$r        <- 1/(21/365)  # maternally-derived IgG half life of 3 weeks (Villard et al., 2016. Diagnostic microbiology and infectious disease;84(1):22-33)
pars$mctr     <- c(0.15, 0.44, 0.71)  # SYROCOT 2007. The Lancet 369
pars$burnin   <- 750

# Seroprevalence data
if (getwd()=="/Users/gregorymilne/Desktop/R Projects/stan"){ #local
  df <- readRDS("data/global_data.rds")
  
} else if (getwd()=="/storage/users/gmilne/test"){ #cluster
  df <- readRDS("global_data.rds")
}

countries <- sort(unique(df$country)) #countries in the dataset

pars$country <- "Brazil"
fitting_data <- subset(df, df$country == pars$country)

# Load demographic data
if (getwd()=="/Users/gregorymilne/Desktop/R Projects/stan"){ #local
  source("R files/demogdat.R")
  
} else if (getwd()=="/storage/users/gmilne/test"){ #cluster
  source("demogdat.R")
}

# TEMPORAL foi decline
# Options: "none", "stepwise" or "linear" 
pars$temporal_foi <-  "linear"

# AGE-related foi change
# Options: "constant", "half" or "double")
pars$age_foi      <- "constant"

# trouble shooting parameter
# Options: 1 (plots graph of age-foi profile); 0 (model runs without age-foi graph plot)
pars$troubleshoot <- 0

# For forecasting: if 1, time interval set to be longer
pars$forecast <- 1

## Number of years between 1st & last data time points
pars$tdiff <- max(fitting_data$year) - min(fitting_data$year)

# no. years between each data point from first to last
year_diff <- vector("numeric", length = nrow(fitting_data))
for(i in 1:nrow(fitting_data)-1){
  year_diff[i] <- fitting_data$year[i+1] - fitting_data$year[i]
}

year_diff <- c(pars$burnin, year_diff[1:(length(year_diff)-1)]) # add burnin time
year_diff <- cumsum(year_diff)  # cumulative sum to get model sampling times

### interpolated parameters from demographic data
# Population size #
f <- spline(age_pop, tot_pop, xout=pars$age)

pars$Na <- f$y*(sum(tot_pop[1:which(age_pop==pars$amax-2.5)])/sum(f$y)) # select correct age specific population size depending on choice of pars$amax
pars$N  <- sum(pars$Na) # total population size

# Mortality rate #
pars$d <- spline(x=age_mort, y=mort, xout=pars$age)$y # interpolated age specific per capita mortality rates
pars$d[which(pars$d<0)] <- 0  # zero minus elements

# Birth rate #
pars$propfert <- approx(age_fert, prop_fert_age, xout=pars$age)$y
pars$propfert[is.na(pars$propfert)] <- 0
pars$propfert <- pars$propfert/sum(pars$propfert)

# rm(f, age_pop, tot_pop, age_mort, mort, age_fert, births_age)

## Initial values (state variables)
S0 <- I0 <- Im0 <- vector("numeric", length=pars$agrps)

# entire population initially susceptible
S0[1:length(S0)] <- pars$Na
y <- c(S0, I0, Im0)

# set time for model running
if (pars$forecast == 0) {
  time <- seq(1, pars$burnin + pars$tdiff, 1)
  
} else if(pars$forecast == 1){
  time <- seq(1, pars$burnin + pars$tdiff + 50, 1)
}
