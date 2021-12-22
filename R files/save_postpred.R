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
post <- readRDS(file = paste("posteriors/", pars$country, "/posteriors_",
                             pars$country, "_t", pars$temporal_foi, "_", "a",
                             pars$age_foi, ".RDS", sep=""))

attach(post)

###################################
## Read in posterior predictions ##
###################################
n_samples <- 1
niter     <- 600

prev_list <- ct_list <- vector("list", length=niter)

for(i in 1:niter){
  
  prev_list[[i]] <- readRDS(file = paste("posteriors/", pars$country,
                                         "/prev_predictions/", "prev_predictions_", 
                                         pars$country, "_t", pars$temporal_foi, "_a",
                                         pars$age_foi, "_", i, ".Rdata", sep = ""))
  
  ct_list[[i]] <- readRDS(file = paste("posteriors/", pars$country, 
                                       "/ct_predictions/", "ct_predictions_", 
                                       pars$country, "_t", pars$temporal_foi, "_a", 
                                       pars$age_foi, "_", i, ".Rdata", sep = ""))
  
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

# call function to find relevant immunoassay sensitivity & specificity values 
diagnostic_values <- find_diagnostic_values(fitting_data, assays)

# calculate true prevalence
tp <- vector("numeric", length=nrow(fitting_data))
for(j in 1:nrow(fitting_data)){
  tp[j] <- (fitting_data$prev[j] + (diagnostic_values$sp[j] - 1)) / (diagnostic_values$sp[j] + (diagnostic_values$se[j] - 1))
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

## to make sure can't get true prevalence < 0 
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
    time = fitting_data$year,      #sampling year
    dat_prev = true_prev$prev*100, #data adjusted to true prevalence (%)
    dat_low = true_prev$ci_lo*100, #data lower interval (%)
    dat_up = true_prev$ci_up *100  #data upper interval (%)
  )


## Prevalence, past & forecasting ##
# time points across which foi is declining
timepoints <- (pars$burnin - (round(max(post$tau), 0))) : max(time)

# equivalent years
years <- (min(fitting_data$year) - (round(max(post$tau), 0))) : (max(fitting_data$year)+pars$years_forecast)


##############################################################################################
## ! IN PROGRESS ! : MAKE CODE THAT DEFINES TIMEPOINTS MORE INTUITIVELY

## This works -- now need to change whole script to use this and check results/plots vs old ones
# 
# # Define years wanted for plotting
# years         <- 1980:2030
# 
# # Calculate equivalent model timepoints
# min_timepoint <- pars$burnin + min(years) - min(fitting_data$year)
# timepoints    <- min_timepoint: (min_timepoint + (max(years) - min(years)))
# 
# 
# # Equivalent model timepoints
# timepoints <- (pars$burnin - (min(fitting_data$year) - min(years))) : 
#   (pars$burnin - (min(fitting_data$year) - min(years)) + (max(years) - min(years)))
# 
# 
# test <- data.frame(
#   time = years,
#   mod_prev     = prev_med[timepoints]  ,  #modelled median prevalence (%)
#   mod_prev_low = prev_lower[timepoints], #modelled prevalence lower interval (%)
#   mod_prev_up  = prev_upper[timepoints]   #modelled prevalence upper interval (%)
# )
##############################################################################################

## Create df for model estimates
# make list to store model output from different countries
if (!exists("prev_all")) {  #only create if list not in existence
  prev_all <- vector("list", length=length(countries))
}

prev_all[[which(countries == pars$country)]] <- data.frame(
  time         = years,                       #years of foi decrease
  mod_prev     = prev_med[timepoints]  *100,  #modelled median prevalence (%)
  mod_prev_low = prev_lower[timepoints]*100, #modelled prevalence lower interval (%)
  mod_prev_up  = prev_upper[timepoints]*100   #modelled prevalence upper interval (%)
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

## Save model estimates & data for plotting
saveRDS(ct_all,   "data/ct_predictions.rds")
saveRDS(prev_all, "data/prev_predictions.rds")
saveRDS(prev_fit, "data/prev_data.rds")
