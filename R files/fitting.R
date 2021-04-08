########### Model fitting
rm(list = ls()) # clear working environment

#directory when using cluster
# setwd("/storage/users/gmilne/test/parallel")

#directory when not using the cluster
setwd("~/Desktop/R Projects/stan")

# Set seed
# SEED = as.numeric(Sys.getenv("SEED")) # cluster
SEED = 1                              # local
set.seed(SEED)

# Load scripts #
# source("demogdat.R")
# source("setparms.R")
# source("seroprev_dat.R")
# source("model.R")
source("R files/demogdat.R")
source("R files/setparms.R")
source("R files/seroprev_dat.R")
source("R files/model.R")

# Load packages
library(deSolve)
library(lhs)
library(dplyr)

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

## Latin hypercube sampling
nsim  <- 100
npars <- 4

par_arr <- randomLHS(nsim, npars) #create parameter array
par_arr[,1] <- log(qunif(par_arr[,1], min=0, max=0.2))        #log lambda0
par_arr[,2] <- log(rbeta(par_arr[,2], shape1 = 2, shape2=80)) #log lambda1
par_arr[,3] <- log(runif(par_arr[,3], min=0, max=0.2))        #log gradient
par_arr[,4] <- log(runif(par_arr[,4], min=0.2, max=0.8))      #log shape

# par(mfrow=c(2,2))
# dummy <- apply(exp(par_arr), 2, hist, main = "") #plot prior distributions

# Run model and store likelihoods
lik_arr <- vector(mode="numeric", length=nsim) #create likelihood array

system.time(
  for(i in 1:nrow(par_arr)){
    pars$log.lambda0  <- par_arr[i,1]
    pars$log.lambda1  <- par_arr[i,2]
    pars$log.gradient <- par_arr[i,3]
    pars$log.shape    <- par_arr[i,4]
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

## save parameter values & likelihoods in one dataframe
# likpar_arr <- data.frame(par_arr, lik_arr)
# names(likpar_arr) <- c("log.lambda0", "log.lambda1", "log.gradient", "log.shape", "log.lik")
# saveRDS(likpar_arr, file = paste("parliks_", SEED, ".Rdata", sep = ""))

## Read in and combine all parameter/likelihood values ##
niter <- 20  #no. iterations on the cluster
out_likpar <- data.frame(matrix(ncol=npars+1, nrow=niter*nsim))
names(out_likpar) <- c("log.lambda0", "log.lambda1", "log.gradient", "log.shape", "likelihood")
counter <- seq(1, (nsim*niter), by = 1/nsim)  # used in for loop to pick correct parliks_ file

#works for multiples of 10
for(i in seq(1, niter*nsim, by=nsim)){
  if(i==1){
    out_likpar[i:(i+nsim-1),] <- readRDS(file = paste("mod_output/parliks_", i, ".RData", sep = ""))
  } else if (i > 1 & i < (niter*nsim)){
    out_likpar[i:(i+nsim-1),] <- readRDS(file = paste("mod_output/parliks_", i-(i-counter[i]), ".RData", sep = ""))
  }
}

## plot parameter values vs. likelihood
par(mfrow=c(2,2))
plot(exp(out_likpar$log.lambda0), out_likpar$likelihood)
abline(lm(out_likpar$likelihood ~ exp(out_likpar$log.lambda0)), col="red")

plot(exp(out_likpar$log.lambda1), out_likpar$likelihood)
abline(lm(out_likpar$likelihood ~ exp(out_likpar$log.lambda1)), col="red")

plot(exp(out_likpar$log.gradient), out_likpar$likelihood)
abline(lm(out_likpar$likelihood ~ exp(out_likpar$log.gradient)), col="red")

plot(exp(out_likpar$log.shape), out_likpar$likelihood)
abline(lm(out_likpar$log.shape ~ exp(out_likpar$log.gradient)), col="red")

# See what the best fitting foi form looks like
x<-exp(out_likpar[which.min(out_likpar$likelihood),])
lambda0 <-  x[1]$log.lambda0
lambda1 <-  x[2]$log.lambda1
gradient <- x[3]$log.gradient
shape <-    x[4]$log.shape
par(mfrow=c(2,1))
foi <- lambda0 + lambda1 * (pars$age * exp(-gradient*pars$age))
plot(pars$age, foi, type='l')
foi <- (lambda0 + lambda1 * (pars$age * exp(-gradient*pars$age)))*shape
plot(pars$age, foi, type='l', lty=2)


# Get x number of best-fitting par sets & their indices
bfit <- sort(out_likpar$likelihood, index.return=T)
top_pars <- out_likpar[head(bfit$ix, 10),]

# see what range the best-fit pars span
par(mfrow=c(2,2))
plot(exp(top_pars$log.lambda0), top_pars$likelihood)
plot(exp(top_pars$log.lambda1), top_pars$likelihood)
plot(exp(top_pars$log.gradient), top_pars$likelihood)
plot(exp(top_pars$log.shape), top_pars$likelihood)


## simulate model & estimate cases of CT for best fitting parameter set
bfit <- out_likpar[which.min(out_likpar$likelihood),]
pars$log.lambda0  <- bfit$log.lambda0
pars$log.lambda1  <- bfit$log.lambda1
pars$log.gradient <- bfit$log.gradient
pars$log.shape    <- bfit$log.shape
sol <- ode(y = y, times = time, parms = pars,  func = age_si)  #save model solution

#Plot first time point vs. data
df <- getit(850)
par(mfrow=c(1,1))
plot(df[,"obs_pI"]~df[,"a"], type='l', xlab="Age (years)", ylab="Prevalence")
points(neth_95$age_mid, neth_95$prevalence)
# how many cases of CT that year?
sol[850, "ctt"]

#Plot second time point vs. data
df2 <- getit(850+11)
par(mfrow=c(1,1))
plot(df2[,"obs_pI"]~df2[,"a"], type='l', xlab="Age (years)", ylab="Prevalence")
points(neth_06$age_mid, neth_06$prevalence)
# how many cases of CT that year?
sol[(850+11), "ctt"]



# plot(df[,"Na"]~df[,"a"], type='l', ylim=c(0,70000))
# points(pars$age, pars$Na)
# lines(rep(15, 401), seq(0,60000, by=60000/400), col="red", lty=2)
# lines(rep(40, 401), seq(0,60000, by=60000/400), col="red", lty=2)


## Plot priors ##
# nsim <- 1000
# set prior on lambda0 #
# set.seed(1001)
# plot(density(par_arr[,1]), main = "lambda0 prior")
# polygon(density(par_arr[,1]), col = "lightblue")

# set prior on lambda1 #
# plot(density(par_arr[,2]), main = "lambda1 prior")
# polygon(density(par_arr[,2]), col = "lightblue")

# set prior on gradient #
# plot(density(par_arr[,3]), main = "gradient prior")
# polygon(density(par_arr[,3]), col = "lightblue")
