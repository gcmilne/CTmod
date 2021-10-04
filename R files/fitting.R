#######################################
########### Model fitting #############
#######################################
rm(list = ls()) # clear working environment

###################
## Load packages ##
###################
library(deSolve)
library(lhs)
library(dplyr)

###################
## Set directory ##
###################
# Sets working directory based on environment (Cluster (RVC or UCL), or local). 
# Options: "RVC", "UCL", "none"
cluster <- "UCL"

if (cluster == "RVC") {
  setwd("/storage/users/gmilne/test")  #cluster (RVC)
  
} else if (cluster == "UCL") { 
  setwd("/lustre/scratch/scratch/ucbtgmi")  #cluster (UCL)
  
} else if (cluster == "none") { 
  setwd("~/GitHub/toxCTmod")  #local
  
}

##################
## Load scripts ##
##################
if (cluster == "none") { #local
  source("R files/setparms.R")
  source("R files/diagnostics.R")
  source("R files/model.R")
  source("R files/funcs.R")
  
} else if (cluster == "RVC" | cluster == "UCL"){ #cluster
  source("setparms.R")
  source("diagnostics.R")
  source("model.R")
  source("funcs.R")
}

##############
## Set seed ##
##############
if (cluster == "none") { #local
  SEED = 1
  
} else if (cluster == "RVC" | cluster == "UCL"){ #cluster
  SEED = as.numeric(Sys.getenv("SEED"))
}

set.seed(SEED)

#####################
## Additional data ##
#####################
if (cluster == "UCL") {
  niter <- 1000                 # no. jobs to submit
  tot_iter <- 60000             # total number of samples to run

} else if (cluster == "RVC") {
  niter <- 68                   # no. jobs to submit
  tot_iter <- 60000             # total number of samples to run

}

nsim <- ceiling(tot_iter / niter) # no. iterations on each for loop

## no. parameters to fit
npars <- 3

## Change objects according to dataset being fit:
matched_indices <- which(pars$propfert !=0)  #child-bearing ages

logliks <- vector("numeric", length=nrow(fitting_data))  #initialise log likelihood vector

##############################
## Latin hypercube sampling ##
##############################
par_arr <- randomLHS(nsim, npars) #create parameter array
par_arr[,1] <- log(qunif(par_arr[,1], min=0, max=0.2))  #log lambda0
par_arr[,2] <- log(qunif(par_arr[,2], min=0, max=2))    #log beta
par_arr[,3] <- log(qunif(par_arr[,3], min=0, max=100))  #log tau

colnames(par_arr) <- c("log.lambda0", "log.beta", "log.tau")

# par(mfrow=c(3,2))
# dummy <- apply(exp(par_arr), 2, hist, main = "") #plot prior distributions

#####################################
## Run model and store likelihoods ##
#####################################
lik_arr <- vector(mode="numeric", length=nsim) #create likelihood array

# Create list to store model output for each of the timepoints being fit
matched_prev      <- vector("list", length = nrow(fitting_data))
mean_matched_prev <- vector("list", length = nrow(fitting_data))


system.time (
  
  for (i in 1:nrow(par_arr)) {
    
    # extract parameter values from priors
    pars$log.lambda0  <- par_arr[i,"log.lambda0"]
    pars$log.beta    <- par_arr[i,"log.beta"]
    pars$log.tau <- par_arr[i,"log.tau"]
    
    sol                <- ode(y = y, times = time, parms = pars,  func = age_si)  #save model solution
    
    for (j in 1:nrow(fitting_data)) {
      
      store_sim          <- getit(year_diff[j])  #store age profile at jth timepoint
      
      matched_prev[[j]]  <- store_sim[,"pI"][matched_indices]  #select true prevalence from relevant age categories
      
      
      ## Select relevant sensitivity & specificity values ##
        if (is.na(fitting_data$method2[j]) & is.na(fitting_data$method3[j])) {  #if timepoint used 1 immunoassay type
          
          pars$se <- assays$se [ which(fitting_data$method[j] == assays$method) ]  #se
          pars$sp <- assays$sp [ which(fitting_data$method[j] == assays$method) ]  #sp
          
          
        } else if (!is.na(fitting_data$method2[j]) & is.na(fitting_data$method3[j])) {  # if timepoint used 2 immunoassays
          
          se_method1 <- assays$se [ which(fitting_data$method [j] == assays$method) ]  #se of first immunoassay
          se_method2 <- assays$se [ which(fitting_data$method2[j] == assays$method) ]  #se of second immunoassay
          
          sp_method1 <- assays$sp [ which(fitting_data$method [j] == assays$method) ]  #sp of first immunoassay
          sp_method2 <- assays$sp [ which(fitting_data$method2[j] == assays$method) ]  #sp of second immunoassay
          
          pars$se <- weighted.mean(x = c(se_method1, se_method2), w = c(fitting_data$prop_method1[j], fitting_data$prop_method2[j])) #weighted mean se
          pars$sp <- weighted.mean(x = c(sp_method1, sp_method2), w = c(fitting_data$prop_method1[j], fitting_data$prop_method2[j])) #weighted mean sp
          
          
        } else if (!is.na(fitting_data$method2[j]) & !is.na(fitting_data$method3[j])) {  # if timepoint used 3 immunoassays
          
          se_method1 <- assays$se [ which(fitting_data$method [j] == assays$method) ]  #se of first immunoassay
          se_method2 <- assays$se [ which(fitting_data$method2[j] == assays$method) ]  #se of second immunoassay
          se_method3 <- assays$se [ which(fitting_data$method3[j] == assays$method) ]  #se of third immunoassay
          
          sp_method1 <- assays$sp [ which(fitting_data$method [j] == assays$method) ]  #sp of first immunoassay
          sp_method2 <- assays$sp [ which(fitting_data$method2[j] == assays$method) ]  #sp of second immunoassay
          sp_method3 <- assays$sp [ which(fitting_data$method3[j] == assays$method) ]  #sp of third immunoassay
          
          pars$se <- weighted.mean(x = c(se_method1, se_method2, se_method3), w = c(fitting_data$prop_method1[j], fitting_data$prop_method2[j], fitting_data$prop_method3[j])) #weighted mean se
          pars$sp <- weighted.mean(x = c(sp_method1, sp_method2, sp_method3), w = c(fitting_data$prop_method1[j], fitting_data$prop_method2[j], fitting_data$prop_method3[j])) #weighted mean sp
          
        }
    
    ## Adjust true prevalence by relevant diagnostic se and sp ##
    matched_prev[[j]] <- matched_prev[[j]] * (pars$se + pars$sp - 1) + (1 - pars$sp)
    
    ## Demographically weighted mean of modelled adjusted seroprevalence ##
    mean_matched_prev[j] <- weighted.mean(x = matched_prev[[j]], w = pars$propfert[matched_indices] * pars$Na[matched_indices])
    
    ## Calculate log likelihood ##
    logliks[j] <- loglik(k = fitting_data$k[j], n = fitting_data$n[j], prev = unlist(mean_matched_prev[j]))
    
    }
    
    lik_arr[i]    <- sum(-logliks)  #sum the log likelihoods
    
  }
  
)


## save output ##
likpar_arr <- data.frame(par_arr, lik_arr)

names(likpar_arr)[4] <- "log.lik"

saveRDS(likpar_arr, file = paste("parliks_", pars$country, "_t", pars$temporal_foi, "_a", 
                                 pars$age_foi, "_", SEED, ".Rdata", sep = ""))
