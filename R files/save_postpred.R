###########################################
## Script to save posterior predictions  ##
###########################################

#################
# Load packages #
#################
library(binom)

################
# Load scripts #
################
source("R files/diagnostics.R")
source("R files/setparms.R")

#####################################
## Read in  posterior distribution ##
#####################################
new_fit <- 1 #fitting with shape prior 0-2 (1) or 0-0.8 (0)

rm(post.lambda, post.shape, post.tdecline)

if (new_fit == 0) {
  post <- readRDS(file = paste("posteriors/", pars$country, "/", "posteriors_", pars$country, "_t", pars$temporal_foi, "_", "a", pars$age_foi, ".RDS", sep=""))
} else if (new_fit == 1){
  post <- readRDS(file = paste("posteriors/", pars$country, "/new_fit/", "posteriors_", pars$country, "_t", pars$temporal_foi, "_", "a", pars$age_foi, ".RDS", sep=""))
}

attach(post)

###################################
## Read in posterior predictions ##
###################################

if (new_fit == 0) {
  
  if (pars$country == "Brazil" | pars$country == "Burkina Faso") {
    
    n_samples  <- 6    #no. parameter sets run on each cluster job
    niter      <- 100  #overall no. cluster jobs
    
  } else if (pars$country == "Iran (Islamic Republic of)") {
    
    n_samples  <- 1    #no. parameter sets run on each cluster job
    niter      <- 599  #overall no. cluster jobs
    
  } else if (pars$country != "Brazil" | pars$country != "Burkina Faso" | pars$country != "Iran (Islamic Republic of)" )  { #for the rest of the countries 
    
    n_samples  <- 1    #no. parameter sets run on each cluster job
    niter      <- 600  #overall no. cluster jobs
    
  }
  
} else if (new_fit == 1) {
  n_samples <- 1
  niter     <- 600
}


prev_list <- ct_list <- vector("list", length=niter)

if (new_fit == 0) {
  
  for(i in 1:niter){
    
    prev_list[[i]] <- readRDS(file = paste("posteriors/", pars$country, "/", "prev_predictions/", "prev_predictions_", pars$country, "_t", pars$temporal_foi, "_a",
                                           pars$age_foi, "_", i, ".Rdata", sep = ""))
    
    ct_list[[i]] <- readRDS(file = paste("posteriors/", pars$country, "/", "ct_predictions/", "ct_predictions_", pars$country, "_t", pars$temporal_foi, "_a", 
                                         pars$age_foi, "_", i, ".Rdata", sep = ""))
  }
  
} else if (new_fit == 1) {
  
  for(i in 1:niter){
    
    prev_list[[i]] <- readRDS(file = paste("posteriors/", pars$country, "/new_fit/", "prev_predictions/", "prev_predictions_", pars$country, "_t", pars$temporal_foi, "_a",
                                           pars$age_foi, "_", i, ".Rdata", sep = ""))
    
    ct_list[[i]] <- readRDS(file = paste("posteriors/", pars$country, "/new_fit/", "ct_predictions/", "ct_predictions_", pars$country, "_t", pars$temporal_foi, "_a", 
                                         pars$age_foi, "_", i, ".Rdata", sep = ""))
    
  }
}

####################################################################
## Calculate median, upper & lower percentiles of model estimates ##
####################################################################

# put data lists into dataframes
prev_mat <- do.call(rbind.data.frame, prev_list)
ct_mat <- do.call(rbind.data.frame, ct_list)
ct_med <- ct_lower <- ct_upper <- prev_med <- prev_lower <- prev_upper <- vector("numeric", length=max(time))

# calculate percentiles
for(i in 1:max(time)){
  
  ct_med[i]   <- median(ct_mat[,i])
  ct_lower[i] <- quantile(ct_mat[,i], probs = 0.025)
  ct_upper[i] <- quantile(ct_mat[,i], probs = 0.975)
  
  prev_med[i]   <- median(prev_mat[,i])
  prev_lower[i] <- quantile(prev_mat[,i], probs = 0.025)
  prev_upper[i] <- quantile(prev_mat[,i], probs = 0.975)
}


####################################
## Adjust data to true prevalence ##
####################################
tp <- vector("numeric", length=nrow(fitting_data))

## Extract correct se & sp, & adjust observed prevalence to true prevalence
for(j in 1:nrow(fitting_data)){
  
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
  
  
  tp[j] <- (fitting_data$prev[j] + (pars$sp - 1)) / (pars$sp + (pars$se - 1))
  
  
}

############################################################
## Calculate 95% CIs for data adjusted to true prevalence ##
############################################################

### Estimate new k for true prevalence ###

## find value of k that leads to smallest difference in the calculated true prev estimate
something <- vector("list", length=nrow(fitting_data))
est_prev  <- vector("list", length=nrow(fitting_data))
for(i in 1:length(something)){
  something[[i]] <- seq(-1000,1000)
  est_prev[[i]] <- rep(0, length=nrow(fitting_data))
}

for(i in 1:nrow(fitting_data)){
  for(j in 1:length(something[[i]])){
    est_prev[[i]][j] <- (fitting_data$k[i] + something[[i]][j]) / fitting_data$n[i]
  }
}  

## find best value of k for true prevalence, given the above calculations
k_tp <- vector("numeric", length=nrow(fitting_data))

for(i in 1:nrow(fitting_data)){
  k_tp[i] <- something[[i]] [which.min(abs(est_prev[[i]] - tp[i]))]
}

## recalculate true prevalence given updated k
new_tp <- (fitting_data$k + k_tp) / fitting_data$n
true_prev <- data.frame(prev = new_tp, k = fitting_data$k + k_tp, n = fitting_data$n)

## make sure can't get true prevalence < 0 
for(i in 1:length(true_prev$k)){
  if(true_prev$k[i]   < 0){
    true_prev$k[i]    <- 0
    true_prev$prev[i] <- true_prev$k[i]/true_prev$n[i]
    
  }
}

# plot(true_prev$prev, tp)  #check concordance

## Calculate data 95% CIs
cis <- binom.confint(x=true_prev$k, n=true_prev$n, conf.level=0.95, methods="exact")

# Store lower and upper estimates
true_prev$ci_lo <- cis$lower
true_prev$ci_up <- cis$upper

#####################
## Create datasets ##
#####################

################
## Prevalence ##
################

## Prevalence across years with data ##
if (!exists("prev_fit")) {  #only create if list not in existence
  prev_fit <- vector("list", length=length(countries))
}

prev_fit[[which(countries == pars$country)]] <- 
  data.frame(
    time = fitting_data$year,             #sampling year
    dat_prev = true_prev$prev,            #data adjusted to true prevalence
    dat_low = true_prev$ci_lo,            #data lower interval
    dat_up = true_prev$ci_up,             #data upper interval
    mod_prev = prev_med[year_diff],       #modelled median prevalence
    mod_prev_low = prev_lower[year_diff], #modelled prevalence lower interval
    mod_prev_up = prev_upper[year_diff]   #modelled prevalence upper interval
  )


## Prevalence, past & forecasting ##
# time points across which foi is declining
timepoints <- (pars$burnin - (round(max(exp(post$post.tdecline)), 0))-10) : max(time)

# equivalent years
years <- (min(fitting_data$year) - (round(max(exp(post$post.tdecline)), 0))-10) : (max(fitting_data$year)+pars$years_forecast)

# period from beginning of foi decrease to final data timepoint
sampling_period <- (min(fitting_data$year) - round(exp(pars$log.tdecline), 0)) : max(fitting_data$year)

## Create df for model estimates
# make list to store model output from different countries
if (!exists("prev_all")) {  #only create if list not in existence
  prev_all <- vector("list", length=length(countries))
}

prev_all[[which(countries == pars$country)]] <- data.frame(
  time         = years,                  #years of foi decrease
  mod_prev     = prev_med[timepoints],   #modelled median prevalence
  mod_prev_low = prev_lower[timepoints], #modelled prevalence lower interval
  mod_prev_up  = prev_upper[timepoints]  #modelled prevalence upper interval
)


##################
## CT incidence ##
##################

## CT incidence, past & forecasting ##

# Make CT incidence relative to no. births
fitting_timepoints <- pars$burnin : (pars$burnin + pars$tdiff)
fitting_years <- min(fitting_data$year) : max(fitting_data$year)

#Calculate annual no. births
deaths     <- sum((pars$d) * y[1:pars$agrps])
byebye     <- y[pars$agrps] * pars$da
births_age <- (deaths + byebye) * pars$propfert
births     <- sum(births_age)

ct_rel     <- ct_med   / (births / 10000)  #CT incidence per 10,000 live births
ct_rel_low <- ct_lower / (births / 10000)
ct_rel_up  <- ct_upper / (births / 10000)

if (pars$forecast == 1){
  
  # time points across which foi is declining
  timepoints <- (pars$burnin - (round(max(exp(post$post.tdecline)), 0))-10) : max(time)
  
  # equivalent years
  years <- (min(fitting_data$year) - (round(max(exp(post$post.tdecline)), 0))-10) : (max(fitting_data$year)+pars$years_forecast)
  
  # Create df
  if (!exists("ct_all")) {  #only create if list not in existence
    ct_all <- vector("list", length=length(countries))
  }
  
  ct_all[[which(countries == pars$country)]] <- data.frame(
    time = years,                         #years of foi decrease
    ct_rel = ct_rel[timepoints],          #modelled ct incidence
    ct_rel_low = ct_rel_low[timepoints],  #modelled ct incidence lower interval
    ct_rel_up = ct_rel_up[timepoints]     #modelled ct incidence upper interval
  )
  
  # add sequelae (Torgerson et al., 2013. Bull WHO, 91, pp. 501-508)
  ct_all[[which(countries == pars$country)]]$foetal_loss                 <- 0.024 * ct_all[[which(countries == pars$country)]]$ct_rel
  ct_all[[which(countries == pars$country)]]$neonatal_death              <- 0.007 * ct_all[[which(countries == pars$country)]]$ct_rel
  ct_all[[which(countries == pars$country)]]$intracranial_calcifications <- 0.11  * ct_all[[which(countries == pars$country)]]$ct_rel
  ct_all[[which(countries == pars$country)]]$hydrocephalus               <- 0.02  * ct_all[[which(countries == pars$country)]]$ct_rel
  ct_all[[which(countries == pars$country)]]$cns_abnormalities           <- 0.029 * ct_all[[which(countries == pars$country)]]$ct_rel
  
  if (pars$country != "Brazil") {
    ct_all[[which(countries == pars$country)]]$chorioretinitis_first     <- 0.13  * ct_all[[which(countries == pars$country)]]$ct_rel
    ct_all[[which(countries == pars$country)]]$chorioretinitis_later     <- 0.16  * ct_all[[which(countries == pars$country)]]$ct_rel
    
  }  else if (pars$country == "Brazil") {
    residual = (1 - (0.024 + 0.007  + 0.02 + 0.029)) #to avoid incidence > total cases, apply to residual incidence after sequelae with higher disability weights extracted
    
    ct_all[[which(countries == pars$country)]]$chorioretinitis_first   <- (0.80  * residual) * ct_all[[which(countries == pars$country)]]$ct_rel
    ct_all[[which(countries == pars$country)]]$chorioretinitis_later   <- (0.10  * residual) * ct_all[[which(countries == pars$country)]]$ct_rel
    
  }
  
}

## Save data for plotting
saveRDS(ct_all, "data/ct_predictions.rds")
saveRDS(prev_all, "data/prev_predictions.rds")
