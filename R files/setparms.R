#############################
## SET PARAMETERS & INITs ######
#############################
#directory when not using cluster
source("R files/demogdat.R")
# change of directory for cluster
# source("demogdat.R")


pars  <- list(scale=50, shape= 1, agrps=400, amax=100, lambda0 = 0.05,
              lambda1 = 0, gradient = 0, se=0.98475, sp=0.98725)
pars$la <- (pars$amax/pars$agrps)  # no. years in each age group (here = 1 yr)
pars$da <-  1/pars$la  # ageing rate
pars$age <- seq(0+pars$la/2, pars$amax-pars$la/2, length.out=pars$agrps) # age at midpoints of age groups
pars$r <- 1/(21/365)  # maternally-derived IgG half life of 3 weeks (Villard et al., 2016. Diagnostic microbiology and infectious disease;84(1):22-33)
pars$mctr <- c(0.15, 0.44, 0.71)  # SYROCOT 2007. The Lancet 369
#pars$foi <- NA
# 
# ## Best fit parameters for Netherlands '95/'96 data:
# # lambda0 = 2.235830e+01; lambda1 = 7.570937e-04; gradient = 2.957032e-02; shape = 5.647195e-04
# 
# ### interpolated parameters from demographic data
f <- spline(age_pop, tot_pop, xout=pars$age)

########
#### tot_pop from demogdat.R  -- !! currently length 20 and needs to be length 16 !!
######
pars$Na <- f$y/sum(f$y)*sum(tot_pop) # interpolates age specific population size
pars$N <- pars$Na # total population size
## check how to prevent neg death values
pars$d <- spline(age_mort, mort, xout=pars$age)$y # interpolated age specific per capita mortality rates
# pars$d[length(pars$d)] <- pars$d[length(pars$d)]
pars$propfert <- approx(age_fert,prop_fert_age,  xout=pars$age)$y
pars$propfert[is.na(pars$propfert)] <- 0
pars$propfert <- pars$propfert/sum(pars$propfert)

#rm(f, age_pop, tot_pop, age_mort, mort, age_fert, births_age, pf)

## Initial values (state variables)
S0 <- I0 <- Im0 <- vector("numeric", length=pars$agrps)
# entire population initially susceptible
S0[1:length(S0)] <- pars$Na
y <- c(S0, I0, Im0)
time <- seq(1,10, 1)
#time<-seq(0,250,1)

