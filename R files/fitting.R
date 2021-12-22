#########################################
## Script to fit the model to the data ##
#########################################
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
# Sets working directory based on environment (Cluster (RVC or UCL), or local)
cluster <- "UCL" #options: "RVC", "UCL", "none"

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

# No. of parameters to fit
npars <- 3

# Get indices for child-bearing ages (changes  according to dataset being fit)
matched_indices <- which(pars$propfert !=0)

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

# Call function to store relevant immunoassay sensitivity & specificity values 
diagnostic_values <- find_diagnostic_values(fitting_data, assays)


system.time (
  
  for (i in 1:nrow(par_arr)) {
    
    # Extract parameter values from priors
    pars$log.lambda0 <- par_arr[i,"log.lambda0"]
    pars$log.beta    <- par_arr[i,"log.beta"]
    pars$log.tau     <- par_arr[i,"log.tau"]
    
    # Run model & store solution
    sol              <- ode(y = y, times = time, parms = pars,  func = age_si)
    
    
    for (j in 1:nrow(fitting_data)) {
      
      store_sim          <- getit(year_diff[j])  #store age profile at jth timepoint
      
      matched_prev[[j]]  <- store_sim[,"pI"][matched_indices]  #select true prevalence from relevant age categories
      
      # Calculate observed prevalence using relevant se and sp values
      matched_prev[[j]] <- matched_prev[[j]] * (diagnostic_values$se[j] + diagnostic_values$sp[j] - 1) + (1 - diagnostic_values$sp[j])
      
      # Calculate demographically weighted mean of observed seroprevalence
      mean_matched_prev[j] <- weighted.mean(x = matched_prev[[j]], w = pars$propfert[matched_indices] * pars$Na[matched_indices])
      
      # Calculate log likelihood
      logliks[j] <- loglik(k = fitting_data$k[j], n = fitting_data$n[j], prev = unlist(mean_matched_prev[j]))
      
    }
    
    lik_arr[i]    <- sum(-logliks)  #sum the log likelihoods
    
  }
  
)


## Save output ##
likpar_arr <- data.frame(par_arr, lik_arr)

names(likpar_arr)[4] <- "log.lik"

saveRDS(likpar_arr, file = paste("parliks_", pars$country, "_t", pars$temporal_foi, "_a", 
                                 pars$age_foi, "_", SEED, ".Rdata", sep = ""))
