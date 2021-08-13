###########################
## Posterior Predictions ##
###########################

#######################
#### Set directory ####
#######################
cluster <- "UCL"

if (cluster == "RVC") {
  setwd("/storage/users/gmilne/test")  #cluster (RVC)
  
} else if (cluster == "UCL") { 
  setwd("/lustre/scratch/scratch/ucbtgmi")        #cluster (UCL)
  
} else if (cluster == "none") { 
  setwd("~/GitHub/stan")  #local
  
}

############
# Set seed #
############
if(cluster == "none"){ #local
  SEED = 1
  
} else if (cluster == "RVC" | cluster == "UCL"){ #cluster
  SEED = as.numeric(Sys.getenv("SEED"))
  
}

set.seed(SEED)

##################
## Load scripts ##
##################
if (cluster == "none") { #local
  source("R files/setparms.R")
  source("R files/diagnostics.R")
  source("R files/model.R")
  
} else if (cluster == "RVC" | cluster == "UCL"){ #cluster
  source("setparms.R")
  source("diagnostics.R")
  source("model.R")
}

  #################
  # Load packages #
  #################
  library(deSolve)
  library(lhs)
  library(dplyr)

###############
## Load data ##
###############
## Posterior distributions for given country
if (cluster == "none") { #local
  post <- readRDS(file = paste("posteriors/", pars$country, "/", "posteriors_", pars$country, "_t", pars$temporal_foi, "_", "a", pars$age_foi, ".RDS", sep=""))

} else if (cluster == "RVC" | cluster == "UCL"){ #cluster
  post <- readRDS(file = paste("posteriors_", pars$country, "_t", pars$temporal_foi, "_", "a", pars$age_foi, ".RDS", sep=""))
  
}

attach(post)

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

###################################
## Perform posterior predictions ##
###################################
## No. samples for each cluster job (whole posterior is 600 samples)
n_samples <- 1

if (n_samples > 1) {
  
  if (SEED == 1) {
    lambda.vec   <- post.lambda  [SEED : n_samples]
    shape.vec    <- post.shape   [SEED : n_samples]
    tdecline.vec <- post.tdecline[SEED : n_samples]
    
  } else if (SEED > 1) {
    lambda.vec   <- post.lambda  [ (n_samples*(SEED-1)+1) : (n_samples*SEED) ]
    shape.vec    <- post.shape   [ (n_samples*(SEED-1)+1) : (n_samples*SEED) ]
    tdecline.vec <- post.tdecline[ (n_samples*(SEED-1)+1) : (n_samples*SEED) ]
    
    
  }
  
} else if (n_samples == 1) {
  
  lambda.vec   <- post.lambda  [SEED]
  shape.vec    <- post.shape   [SEED]
  tdecline.vec <- post.tdecline[SEED]
}


## Run model for each parameter set, store estimated prevalence & CT incidence
ct_cases        <- mean_matched_prev <- vector("list", length = n_samples)
matched_prev    <- vector("list", length=max(time))
matched_indices <- which(pars$propfert !=0)  #child-bearing ages

for(i in 1:length(lambda.vec)){
  
  pars$log.lambda0  <- lambda.vec[i]
  pars$log.shape    <- shape.vec[i]
  pars$log.tdecline <- tdecline.vec[i]
  
  sol <- ode(y = y, times = time, parms = pars,  func = age_si)  #save model solution
  
  for(j in 1:max(time)){
    
    ct_cases[[i]][j] <- sol[j, "ctt"] #store CT incidence
    store_sim <- getit(j) #get age profile
    matched_prev[[j]]  <- store_sim[,"pI"][matched_indices]  #select true prevalence from relevant age categories
    mean_matched_prev[[i]][j] <- weighted.mean(x = matched_prev[[j]], w = pars$propfert[matched_indices] * pars$Na[matched_indices]) #calculate demographically weighted mean prevalence
    
  }
}


# Store estimates in matrix (each row = estimates from a different parameter set)
prev_mat <- t(sapply(mean_matched_prev,c))
ct_mat   <- t(sapply(ct_cases,c))

# Save the output
if (cluster == "RVC") {
  
  saveRDS(prev_mat, file = paste("posteriors/", pars$country, "/", "prev_predictions_", pars$country, "_t", pars$temporal_foi, "_a", 
                                 pars$age_foi, "_", SEED, ".Rdata", sep = ""))
  
  saveRDS(ct_mat, file = paste("posteriors/", pars$country, "/", "ct_predictions_", pars$country, "_t", pars$temporal_foi, "_a", 
                               pars$age_foi, "_", SEED, ".Rdata", sep = ""))
  
} else if (cluster == "UCL") { 
  
  saveRDS(prev_mat, file = paste("prev_predictions_", pars$country, "_t", pars$temporal_foi, "_a", 
                                 pars$age_foi, "_", SEED, ".Rdata", sep = ""))
  
  saveRDS(ct_mat, file = paste("ct_predictions_", pars$country, "_t", pars$temporal_foi, "_a", 
                               pars$age_foi, "_", SEED, ".Rdata", sep = ""))
  
}
