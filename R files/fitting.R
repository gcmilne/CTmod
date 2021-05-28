#######################################
########### Model fitting #############
#######################################
rm(list = ls()) # clear working environment

#######################
#### Set directory ####
#######################
# setwd("/storage/users/gmilne/test/parallel")  #cluster
setwd("~/Desktop/R Projects/stan")  #local

############
# Set seed #
############
if(getwd()=="/Users/gregorymilne/Desktop/R Projects/stan"){ #local
  SEED = 1
  
} else if (getwd()=="/storage/users/gmilne/test/parallel"){ #cluster
  SEED = as.numeric(Sys.getenv("SEED"))
}

set.seed(SEED)

##################
## Load scripts ##
##################
if (getwd()=="/Users/gregorymilne/Desktop/R Projects/stan"){ #local
  source("R files/demogdat.R")
  source("R files/setparms.R")
  source("R files/model.R")
  
} else if (getwd()=="/storage/users/gmilne/test/parallel"){ #cluster
  source("demogdat.R")
  source("setparms.R")
  source("model.R")
}

#################
# Load packages #
#################
library(deSolve)
library(lhs)
library(dplyr)

###################
# Additional data #
###################
## no. iterations on each for loop
nsim  <- 100  

## no. parameters to fit
npars <- 3

## Change objects according to dataset being fit:
matched_indices <- which(pars$propfert !=0)  #child-bearing ages

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
  # out <- out[-nrow(out),]
}

## Joint likelihood function: for fitting both timepoints simultaneously

## To do: Update function to adjust number of ks and ns depending on rows of dataset being fit
loglik <- function(k1, n1, prev1, k2, n2, prev2){
  dbinom(k1, n1, prev1, log=T) + dbinom(k2, n2, prev2, log=T)
}

##############################
## Latin hypercube sampling ##
##############################

par_arr <- randomLHS(nsim, npars) #create parameter array
par_arr[,1] <- log(qunif(par_arr[,1], min=0, max=0.2))        #log lambda0
par_arr[,2] <- log(qunif(par_arr[,4], min=0, max=0.8))        #log shape
par_arr[,3] <- log(qunif(par_arr[,5], min=0, max=100))        #log tdecline

# par(mfrow=c(3,2))
# dummy <- apply(exp(par_arr), 2, hist, main = "") #plot prior distributions

#####################################
## Run model and store likelihoods ##
#####################################

lik_arr <- vector(mode="numeric", length=nsim) #create likelihood array


# Create list to store model output for each of the timepoints being fit
matched_prev <- vector("list", length = nrow(fitting_data))


system.time (
  
  for (i in 1:nrow(par_arr)) {
    
    for (j in 1:nrow(fitting_data)) {
      
      pars$log.lambda0  <- par_arr[i,"log.lambda0"]
      pars$log.shape    <- par_arr[i,"log.shape"]
      pars$log.tdecline <- par_arr[i,"log.tdecline"]
      
      
      ## First time point ##
      if (j == 1){ 
        
      sol                <- ode(y = y, times = time, parms = pars, func = age_si)  #save model solution
      store_sim          <- getit(pars$burnin)                                     #store age profile
      matched_prev[[j]]  <- store_sim[,"pI"][matched_indices]                      #select true prevalence from relevant age categories
      
      
      # select relevant sensitivity & specificity values
      if (is.na(fitting_data$method2[j]) & is.na(fitting_data$method3[j])) {  # if timepoint used 1 immunoassay type
        
        pars$se <- se_values$se [ which(fitting_data$method[j] == se_values$method) ]  #se
        pars$sp <- sp_values$sp [ which(fitting_data$method[j] == sp_values$method) ]  #sp
        
        #adjust true prevalence by relevant diagnostic se and sp
        matched_prev[[j]] <- matched_prev[[j]] * (pars$se + pars$sp - 1) + (1 - pars$sp)
        
        
      } else if (!is.na(fitting_data$method2[j]) & is.na(fitting_data$method3[j])) {  # if timepoint used 2 immunoassays
        
        se_method1 <- se_values$se [ which(fitting_data$method [j] == se_values$method) ]  #se of first immunoassay
        se_method2 <- se_values$se [ which(fitting_data$method2[j] == se_values$method) ]  #se of second immunoassay
        
        sp_method1 <- sp_values$sp [ which(fitting_data$method [j] == sp_values$method) ]  #sp of first immunoassay
        sp_method2 <- sp_values$sp [ which(fitting_data$method2[j] == sp_values$method) ]  #sp of second immunoassay
        
        pars$se <- (se_method1 * fitting_data$prop_method1[j]) + (se_method2 * fitting_data$prop_method2[j]) #weighted se
        pars$sp <- (sp_method1 * fitting_data$prop_method1[j]) + (sp_method2 * fitting_data$prop_method2[j]) #weighted sp
        
        #adjust true prevalence by relevant diagnostic se and sp
        matched_prev[[j]] <- matched_prev[[j]] * (pars$se + pars$sp - 1) + (1 - pars$sp)
          
        
      } else if (!is.na(fitting_data$method2[j]) & !is.na(fitting_data$method3[j])) {  # if timepoint used 3 immunoassays
        
        se_method1 <- se_values$se [ which(fitting_data$method [j] == se_values$method) ]  #se of first immunoassay
        se_method2 <- se_values$se [ which(fitting_data$method2[j] == se_values$method) ]  #se of second immunoassay
        se_method3 <- se_values$se [ which(fitting_data$method3[j] == se_values$method) ]  #se of third immunoassay
        
        sp_method1 <- sp_values$sp [ which(fitting_data$method [j] == sp_values$method) ]  #sp of first immunoassay
        sp_method2 <- sp_values$sp [ which(fitting_data$method2[j] == sp_values$method) ]  #sp of second immunoassay
        sp_method3 <- sp_values$sp [ which(fitting_data$method3[j] == sp_values$method) ]  #sp of third immunoassay
        
        pars$se <- (se_method1 * fitting_data$prop_method1[j]) + (se_method2 * fitting_data$prop_method2[j]) + (se_method3 * fitting_data$prop_method3[j]) #weighted se
        pars$sp <- (sp_method1 * fitting_data$prop_method1[j]) + (sp_method2 * fitting_data$prop_method2[j]) + (sp_method3 * fitting_data$prop_method3[j]) #weighted sp
        
        #adjust true prevalence by relevant diagnostic se and sp
        matched_prev[[j]] <- matched_prev[[j]] * (pars$se + pars$sp - 1) + (1 - pars$sp)
        
      }        
      
      
      ## After first time point ##
      } else if (j > 1 ){ 
        
        sol                <- ode(y = y, times = time, parms = pars,  func = age_si)  #save model solution
        store_sim          <- getit(pars$burnin + year_diff[j])  #store age profile at ith timepoint
        matched_prev[[j]]  <- store_sim[,"pI"][matched_indices]  #select true prevalence from relevant age categories
        
        
        # select relevant sensitivity & specificity values
        if (is.na(fitting_data$method2[j]) & is.na(fitting_data$method3[j])) {  # if timepoint used 1 immunoassay type
          
          pars$se <- se_values$se [ which(fitting_data$method[j] == se_values$method) ]  #se
          pars$sp <- sp_values$sp [ which(fitting_data$method[j] == sp_values$method) ]  #sp
          
          #adjust true prevalence by relevant diagnostic se and sp
          matched_prev[[j]] <- matched_prev[[j]] * (pars$se + pars$sp - 1) + (1 - pars$sp)
          
          
        } else if (!is.na(fitting_data$method2[j]) & is.na(fitting_data$method3[j])) {  # if timepoint used 2 immunoassays
          
          se_method1 <- se_values$se [ which(fitting_data$method [j] == se_values$method) ]  #se of first immunoassay
          se_method2 <- se_values$se [ which(fitting_data$method2[j] == se_values$method) ]  #se of second immunoassay
          
          sp_method1 <- sp_values$sp [ which(fitting_data$method [j] == sp_values$method) ]  #sp of first immunoassay
          sp_method2 <- sp_values$sp [ which(fitting_data$method2[j] == sp_values$method) ]  #sp of second immunoassay
          
          pars$se <- (se_method1 * fitting_data$prop_method1[j]) + (se_method2 * fitting_data$prop_method2[j]) #weighted se
          pars$sp <- (sp_method1 * fitting_data$prop_method1[j]) + (sp_method2 * fitting_data$prop_method2[j]) #weighted sp
          
          #adjust true prevalence by relevant diagnostic se and sp
          matched_prev[[j]] <- matched_prev[[j]] * (pars$se + pars$sp - 1) + (1 - pars$sp)
          
          
        } else if (!is.na(fitting_data$method2[j]) & !is.na(fitting_data$method3[j])) {  # if timepoint used 3 immunoassays
          
          se_method1 <- se_values$se [ which(fitting_data$method [j] == se_values$method) ]  #se of first immunoassay
          se_method2 <- se_values$se [ which(fitting_data$method2[j] == se_values$method) ]  #se of second immunoassay
          se_method3 <- se_values$se [ which(fitting_data$method3[j] == se_values$method) ]  #se of third immunoassay
          
          sp_method1 <- sp_values$sp [ which(fitting_data$method [j] == sp_values$method) ]  #sp of first immunoassay
          sp_method2 <- sp_values$sp [ which(fitting_data$method2[j] == sp_values$method) ]  #sp of second immunoassay
          sp_method3 <- sp_values$sp [ which(fitting_data$method3[j] == sp_values$method) ]  #sp of third immunoassay
          
          pars$se <- (se_method1 * fitting_data$prop_method1[j]) + (se_method2 * fitting_data$prop_method2[j]) + (se_method3 * fitting_data$prop_method3[j]) #weighted se
          pars$sp <- (sp_method1 * fitting_data$prop_method1[j]) + (sp_method2 * fitting_data$prop_method2[j]) + (sp_method3 * fitting_data$prop_method3[j]) #weighted sp
          
          #adjust true prevalence by relevant diagnostic se and sp
          matched_prev[[j]] <- matched_prev[[j]] * (pars$se + pars$sp - 1) + (1 - pars$sp)
          
        }
        
      }    
      
      
      # simple mean of modelled adjusted seroprevalence
      for (j in 1:length(matched_prev)) {
        mean_matched_prev[j] <- mean(matched_prev[[j]])
      }
      
      ## To do: Update function to adjust number of ks and ns depending on rows of dataset being fit
      logliks       <- loglik(k1 = nz1$k, n1 = nz1$n, prev1 = matched_prev, 
                              k2 = nz2$k, n2 = nz2$n, prev2 = matched_prev2)
      lik_arr[i]    <- sum(-logliks)
      
    }
  }
)

## save output ##
# create dataframe
likpar_arr <- data.frame(par_arr, lik_arr)

# name dataframe
if (pars$constant==1){
  names(likpar_arr) <- c("log.lambda0", "log.shape", "log.tdecline", "log.lik")
} else if (pars$constant==0){
  names(likpar_arr) <- c("log.lambda0", "log.lambda1", "log.gradient", "log.shape", "log.tdecline", "log.lik")
}

# save dataframe
if (pars$constant==1){
  
  if (pars$stepwise==1){
    saveRDS(likpar_arr, file = paste("parliks_consant_stepwise_", SEED, ".Rdata", sep = ""))
  } else if (pars$stepwise==0){
    saveRDS(likpar_arr, file = paste("parliks_consant_linear_", SEED, ".Rdata", sep = ""))
  }
  
} else if (pars$constant==0){
  
  if (pars$stepwise==1){
    saveRDS(likpar_arr, file = paste("parliks_varying_stepwise_", SEED, ".Rdata", sep = ""))
  } else if (pars$stepwise==0){
    saveRDS(likpar_arr, file = paste("parliks_varying_linear_", SEED, ".Rdata", sep = ""))
  }
}
