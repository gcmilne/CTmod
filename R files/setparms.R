#############################
## SET PARAMETERS & INITs ######
#############################
#directory when not using cluster
source("R files/demogdat.R")
# change of directory for cluster
# source("demogdat.R")

pars  <- list(amax=100, 
              log.lambda0 = log(0.02),
              log.lambda1 = log(0), 
              log.gradient = log(0.20),
              log.shape= log(0.60),
              log.tdecline = log(20), 
              se=0.98475, sp=0.98725)

pars$agrps <- pars$amax*4
pars$la <- (pars$amax/pars$agrps)  # no. years in each age group
pars$da <-  1/pars$la  # ageing rate
pars$age <- seq(0+pars$la/2, pars$amax-pars$la/2, length.out=pars$agrps) # age at midpoints of age groups
pars$r <- 1/(21/365)  # maternally-derived IgG half life of 3 weeks (Villard et al., 2016. Diagnostic microbiology and infectious disease;84(1):22-33)
pars$mctr <- c(0.15, 0.44, 0.71)  # SYROCOT 2007. The Lancet 369
pars$burnin <- 100
pars$stepwise <- 0 # controls whether foi decrease is stepwise (1) or linear (0)
pars$constant <- 1 # controls whether foi is constant over age (1) or age-varying (0)

## Number of years between 1st & last data time points (change depending on the data being fit!)
pars$tdiff <- 20
 
### interpolated parameters from demographic data
# Population size #
f <- spline(age_pop, tot_pop, xout=pars$age)

pars$Na <- f$y*(sum(tot_pop[1:which(age_pop==pars$amax-2.5)])/sum(f$y)) # select correct age specific population size depending on choice of pars$amax
pars$N <- sum(pars$Na) # total population size

# Mortality rate #
pars$d <- spline(x=age_mort, y=mort, xout=pars$age)$y # interpolated age specific per capita mortality rates
pars$d[which(pars$d<0)] <- 0  # zero minus elements

# Birth rate #
pars$propfert <- approx(age_fert,prop_fert_age,  xout=pars$age)$y
pars$propfert[is.na(pars$propfert)] <- 0
pars$propfert <- pars$propfert/sum(pars$propfert)

# rm(f, age_pop, tot_pop, age_mort, mort, age_fert, births_age)

## Initial values (state variables)
S0 <- I0 <- Im0 <- vector("numeric", length=pars$agrps)

# entire population initially susceptible
S0[1:length(S0)] <- pars$Na
y <- c(S0, I0, Im0)

# set time for model running (burnin period = 850 years)
time <- seq(1,pars$burnin+pars$tdiff, 1)
