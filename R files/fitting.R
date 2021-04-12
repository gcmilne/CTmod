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
par_arr[,2] <- log(qunif(par_arr[,2], min=0, max=0.1))        #log lambda1
par_arr[,3] <- log(qunif(par_arr[,3], min=0, max=1))          #log gradient
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
saveRDS(likpar_arr, file = paste("parliks5-r2_", SEED, ".Rdata", sep = ""))

##############################################################
## POST-CLUSTER: Read in parameter sets & likelihood values ##
##############################################################
out_likpar <- data.frame(matrix(ncol=npars+1, nrow=niter*nsim))
names(out_likpar) <- c("log.lambda0", "log.lambda1", "log.gradient", "log.shape", "log.tdecline", "likelihood")
counter <- seq(1, (nsim*niter), by = 1/nsim)  # used in for loop to pick correct parliks_ file

#works for multiples of 10
i_seq <- seq(1, niter*nsim, by=nsim)
for(i in i_seq[1:26]){
  if(i==1){
    out_likpar[i:(i+nsim-1),] <- readRDS(file = paste("mod_output/5pars_fit2/parliks5-r2_", i, ".RData", sep = ""))
  } else if (i > 1 & i < (niter*nsim)){
    out_likpar[i:(i+nsim-1),] <- readRDS(file = paste("mod_output/5pars_fit2/parliks5-r2_", i-(i-counter[i]), ".RData", sep = ""))
  }
}

for(i in i_seq[28:100]){
  if(i==1){
    out_likpar[i:(i+nsim-1),] <- readRDS(file = paste("mod_output/5pars_fit2/parliks5-r2_", i, ".RData", sep = ""))
  } else if (i > 1 & i < (niter*nsim)){
    out_likpar[i:(i+nsim-1),] <- readRDS(file = paste("mod_output/5pars_fit2/parliks5-r2_", i-(i-counter[i]), ".RData", sep = ""))
  }
}

##########################
## Exploratory plotting ##
##########################
## Plot parameter values vs. likelihood
par(mfrow=c(2,2))
plot(exp(out_likpar$log.lambda0), out_likpar$likelihood)
abline(lm(out_likpar$likelihood ~ exp(out_likpar$log.lambda0)), col="red")

plot(exp(out_likpar$log.lambda1), out_likpar$likelihood)
abline(lm(out_likpar$likelihood ~ exp(out_likpar$log.lambda1)), col="red")

plot(exp(out_likpar$log.gradient), out_likpar$likelihood)
abline(lm(out_likpar$likelihood ~ exp(out_likpar$log.gradient)), col="red")

plot(exp(out_likpar$log.shape), out_likpar$likelihood)
abline(lm(out_likpar$likelihood ~ exp(out_likpar$log.shape)), col="red")

plot(exp(out_likpar$log.tdecline), out_likpar$likelihood)
abline(lm(out_likpar$likelihood ~ exp(out_likpar$log.tdecline)), col="red")

## Plot the best fitting foi form
x<-exp(out_likpar[which.min(out_likpar$likelihood),])
lambda0 <-  x[1]$log.lambda0
lambda1 <-  x[2]$log.lambda1
gradient <- x[3]$log.gradient
shape <-    x[4]$log.shape
par(mfrow=c(1,1))
foi <- lambda0 + lambda1 * (pars$age * exp(-gradient*pars$age))
plot(pars$age, foi, type='l', ylim=c(0,0.08))
foi <- (lambda0 + lambda1 * (pars$age * exp(-gradient*pars$age)))*shape
lines(pars$age, foi, type='l', lty=2)

## Get x number of best-fitting par sets
top_pars <- head(out_likpar[order(out_likpar$likelihood),], 1000)

## See what range the best-fit pars span
par(mfrow=c(3,2))
plot(exp(top_pars$log.lambda0), top_pars$likelihood)
abline(lm(top_pars$likelihood ~ exp(top_pars$log.lambda0)), col="red")

plot(exp(top_pars$log.lambda1), top_pars$likelihood)
abline(lm(top_pars$likelihood ~ exp(top_pars$log.lambda1)), col="red")

plot(exp(top_pars$log.gradient), top_pars$likelihood)
abline(lm(top_pars$likelihood ~ exp(top_pars$log.gradient)), col="red")

plot(exp(top_pars$log.shape), top_pars$likelihood)
abline(lm(top_pars$likelihood ~ exp(top_pars$log.shape)), col="red")

plot(exp(top_pars$log.tdecline), top_pars$likelihood)
abline(lm(top_pars$likelihood ~ exp(top_pars$log.tdecline)), col="red")


## Plot the distribution of 1,000 best-fit par sets (near prior boundaries?)
top_pars <- head(out_likpar[order(out_likpar$likelihood),], 1000)

par(mfrow=c(3,2))
hist(exp(top_pars$log.lambda0),  main = "lambda0 (1,000 best fit)", xlab = "", ylab = "")
hist(exp(top_pars$log.lambda1),  main = "lambda1 (1,000 best fit)", xlab = "", ylab = "")
hist(exp(top_pars$log.gradient), main = "gradient (1,000 best fit)", xlab = "", ylab = "")
hist(exp(top_pars$log.shape),    main = "shape (1,000 best fit)", xlab = "", ylab = "")
hist(exp(top_pars$log.tdecline), main = "tdecline (1,000 best fit)", xlab = "", ylab = "")


#####################################################################
## simulate model & show fit to data for a range of parameter sets ##
#####################################################################
params <- data.frame(matrix(nrow=3, ncol=npars+1))
names(params) <- c("log.lambda0", "log.lambda1", "log.gradient", "log.shape", "log.tdecline", "likelihood")
sorted_lik    <- sort(out_likpar$likelihood, index.return=T)$ix #sort by likelihood value
params[1:3,]  <- out_likpar[head(sorted_lik, 3),]               #best fit
# params[3,]    <- out_likpar[tail(sorted_lik, 1),]               #worst fit
# params[4,]    <- out_likpar[sorted_lik[nrow(out_likpar)/2],]    #middling fit

prev_list1 <- rep(list(matrix(nrow=length(neth_95$age_mid), ncol=1)), nrow(params)) #for timepoint 1
prev_list2 <- rep(list(matrix(nrow=length(neth_95$age_mid), ncol=1)), nrow(params)) #for timepoint 2

for(i in 1:nrow(params)){
  pars$log.lambda0  <- params[i, "log.lambda0"]
  pars$log.lambda1  <- params[i, "log.lambda1"]
  pars$log.gradient <- params[i, "log.gradient"]
  pars$log.shape    <- params[i, "log.shape"]
  pars$log.tdecline <- params[i, "log.tdecline"]
  sol <- ode(y = y, times = time, parms = pars,  func = age_si)  #save model solution
  df  <- getit(pars$burnin)
  prev_list1[[i]]  <- df[,"obs_pI"][matched_indices]  #select observed prevalence from relevant age categories
  df2 <- getit(pars$burnin+11)
  prev_list2[[i]]  <- df2[,"obs_pI"][matched_indices]  #select observed prevalence from relevant age categories
}

#plot first time point vs. data
par(mfrow=c(2,1))
plot(pars$age[matched_indices], prev_list1[[1]], type='l', ylim=c(0,1), ylab="Prevalence",
     xlab="Age (years)", main = "Netherlands 1995")
points(neth_95$age_mid, neth_95$prevalence)
lines(pars$age[matched_indices], prev_list1[[2]], col="grey")
lines(pars$age[matched_indices], prev_list1[[3]], col="grey")
# lines(pars$age[matched_indices], prev_list1[[4]], col="lightgrey")

# plot second time point vs. data
plot(pars$age[matched_indices], prev_list2[[1]], type='l', ylim=c(0,1), ylab="Prevalence",
     xlab="Age (years)", main = "Netherlands 2006")
points(neth_06$age_mid, neth_06$prevalence)
lines(pars$age[matched_indices], prev_list2[[2]], col="grey")
lines(pars$age[matched_indices], prev_list2[[3]], col="grey")
# lines(pars$age[matched_indices], prev_list2[[4]], col="lightgrey")
plot(sol[,"time"], sol[,"ctt"])
sol[850+11,"ctt"]

# head(sort(exp(out_likpar$log.lambda1)))

#################
## Plot the range of prevalence curves over the course of temporal foi decrease
#################
par(mfrow=c(1,1))
df<-getit(pars$burnin)
plot(df[,"a"], df[,"obs_pI"], type='l')
for(i in exp(pars$log.tdecline):-11){
  df  <- getit(pars$burnin-i)
  lines(df[,"a"], df[,"obs_pI"])
}
df<-getit(pars$burnin)
lines(df[,"a"], df[,"obs_pI"], type='l', col="red")  #first timepoint
df<-getit(pars$burnin+11)
lines(df[,"a"], df[,"obs_pI"], type='l', col="red")  #second timepoint

for(i in 1:11){
  df<-getit(pars$burnin+i)
  print(df[,"ct1"])
  
}

sol[(pars$burnin):(pars$burnin+11), "ctt"] ## CT cases over fitting time period

########################################
### sample 1000 foi profiles & plot ####
########################################
par(mfrow=c(1,1))
foi <- (exp(out_likpar$log.lambda0[1]) + exp(out_likpar$log.lambda1[1]) * (pars$age * exp(-exp(out_likpar$log.gradient[1])*pars$age)))*exp(out_likpar$log.shape[1]) #decrease foi after burnin
plot(pars$age, foi, type='l')
foi <- vector("list", length=1000)
for(i in 1:1000){
  foi[[i]] <- (exp(out_likpar$log.lambda0[i]) + exp(out_likpar$log.lambda1[i]) * (pars$age * exp(-exp(out_likpar$log.gradient[i])*pars$age)))*exp(out_likpar$log.shape[i]) #decrease foi after burnin
  lines(pars$age, foi[[i]])
}

#################
## Plot priors ##
#################
# nsim <- 1000
# set prior on lambda0 #
# set.seed(1001)
# par(mfrow=c(2,2))
# plot(density(exp(out_likpar$log.lambda0)), main = "lambda0 prior", xlab="")
# polygon(density(exp(out_likpar$log.lambda0)), col = "lightblue")

# set prior on lambda1 #
# plot(density(exp(out_likpar$log.lambda1)), main = "lambda1 prior", xlab="")
# polygon(density(exp(out_likpar$log.lambda1)), col = "lightblue")

# set prior on gradient #
# plot(density(exp(out_likpar$log.gradient)), main = "gradient prior", xlab="")
# polygon(density(exp(out_likpar$log.gradient)), col = "lightblue")

# set prior on shape #
# plot(density(exp(out_likpar$log.shape)), main = "shape prior", xlab="")
# polygon(density(exp(out_likpar$log.shape)), col = "lightblue")

# set prior on tdecline #
# plot(density(exp(out_likpar$log.tdecline)), main = "tdecline prior", xlab="")
# polygon(density(exp(out_likpar$log.tdecline)), col = "lightblue")
