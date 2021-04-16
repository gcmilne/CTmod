#######################################
########### Model fitting #############
#######################################

# rm(list = ls()) # clear working environment

#######################
#### Set directory ####
#######################
setwd("/storage/users/gmilne/test/parallel")  #cluster
# setwd("~/Desktop/R Projects/stan")  #local

############
# Set seed #
############
SEED = as.numeric(Sys.getenv("SEED"))  #cluster
# SEED = 1                              #local
set.seed(SEED)

##################
## Load scripts ##
##################
source("demogdat.R")      #cluster
source("setparms.R")      #cluster
source("seroprev_dat.R")  #cluster
source("model.R")         #cluster
# source("R files/demogdat.R")      #local
# source("R files/setparms.R")      #local
# source("R files/seroprev_dat.R")  #local
# source("R files/model.R")         #local

#################
# Load packages #
#################
library(deSolve)
library(lhs)
library(dplyr)

###################
# Additional data #
###################
nsim  <- 100   #no. iterations on each for loop
niter <- 100  #no. iterations on the cluster
npars <- 5    #no. parameters to fit

###############
## Functions ##
###############
## Function to get age profiles at given time
getit <- function(time) {
  row <- which(abs(sol[,"time"]-time)==min(abs(sol[,"time"]-time)))
  df <- sol[row,-1]
  S  <- df[1:pars$agrps]  
  I  <- df[(pars$agrps+1):(2*pars$agrps)]  
  Im <- df[(2*pars$agrps+1):(3*pars$agrps)]
  pI <- df[(3*pars$agrps+1):(4*pars$agrps)]
  obs_pI <- df[(4*pars$agrps+1):(5*pars$agrps)]
  dprev <- df[(5*pars$agrps+1):(6*pars$agrps)]
  seroconv1 <- df[(6*pars$agrps+1):(7*pars$agrps)] 
  seroconv2 <- df[(7*pars$agrps+1):(8*pars$agrps)] 
  seroconv3 <- df[(8*pars$agrps+1):(9*pars$agrps)] 
  matAb1 <- df[(9*pars$agrps+1):(10*pars$agrps)]
  matAb2 <- df[(10*pars$agrps+1):(11*pars$agrps)]
  matAb3 <- df[(11*pars$agrps+1):(12*pars$agrps)]
  ct1 <- df[(12*pars$agrps+1):(13*pars$agrps)]
  ct2 <- df[(13*pars$agrps+1):(14*pars$agrps)]
  ct3 <- df[(14*pars$agrps+1):(15*pars$agrps)]
  Na <- df[(15*pars$agrps+1):(16*pars$agrps)]
  age <- pars$age  
  out <- data.frame(a=age, I=I,Im=Im, pI=pI, obs_pI=obs_pI, dprev=dprev, 
                    seroconv1=seroconv1, seroconv2=seroconv2, seroconv3=seroconv3, 
                    matAb1=matAb1, matAb2=matAb2, matAb3=matAb3, 
                    ct1=ct1, ct2=ct2, ct3=ct3, Na=Na)
  ## remove last age category for aesthetics
  out <- out[-nrow(out),]
}

## Joint likelihood function: for fitting both timepoints simultaneously
loglik <- function(k1, n1, prev1, k2, n2, prev2){
  dbinom(k1, n1, prev1, log=T) + dbinom(k2, n2, prev2, log=T)
}

##############################
## Latin hypercube sampling ##
##############################
par_arr <- randomLHS(nsim, npars) #create parameter array
par_arr[,1] <- log(qunif(par_arr[,1], min=0, max=0.2))        #log lambda0
par_arr[,2] <- log(qunif(par_arr[,2], min=0, max=0.2))        #log lambda1
par_arr[,3] <- log(qunif(par_arr[,3], min=0, max=2))          #log gradient
par_arr[,4] <- log(qunif(par_arr[,4], min=0, max=0.8))        #log shape
par_arr[,5] <- log(qunif(par_arr[,5], min=0, max=100))        #log tdecline

# par(mfrow=c(3,2))
# dummy <- apply(exp(par_arr), 2, hist, main = "") #plot prior distributions

#####################################
## Run model and store likelihoods ##
#####################################
lik_arr <- vector(mode="numeric", length=nsim) #create likelihood array

system.time(
  for(i in 1:nrow(par_arr)){
    pars$log.lambda0  <- par_arr[i,1]
    pars$log.lambda1  <- par_arr[i,2]
    pars$log.gradient <- par_arr[i,3]
    pars$log.shape    <- par_arr[i,4]
    pars$log.tdecline <- par_arr[i,5]
    sol           <- ode(y = y, times = time, parms = pars,  func = age_si)  #save model solution
    store_sim     <- getit(pars$burnin)  #store age profile after burnin period (TIMEPOINT 1, 95/96)
    matched_prev  <- store_sim[,"obs_pI"][matched_indices]  #select observed prevalence from relevant age categories
    store_sim2    <- getit(pars$burnin+11)  #store age profile 11 years after burnin period (TIMEPOINT 2, 06/07)
    matched_prev2 <- store_sim2[,"obs_pI"][matched_indices]  #select observed prevalence from relevant age categories
    logliks       <- loglik(k1 = neth_95$k, n1 = neth_95$n, prev1 = matched_prev, 
                            k2 = neth_06$k, n2 = neth_06$n, prev2 = matched_prev2)
    lik_arr[i]    <- sum(-logliks)
  }
)

likpar_arr <- data.frame(par_arr, lik_arr)
names(likpar_arr) <- c("log.lambda0", "log.lambda1", "log.gradient", "log.shape", "log.tdecline", "log.lik")
saveRDS(likpar_arr, file = paste("parliks5_linear_", SEED, ".Rdata", sep = ""))
