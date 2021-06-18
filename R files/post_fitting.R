######################################
## POST-CLUSTER analysis & plotting ##
######################################

##################
## Load scripts ##
##################
source("R files/seroprev_dat.R")
source("R files/demogdat.R")
source("R files/setparms.R")
source("R files/model.R")   
source("R files/diagnostics.R")   

###################
## Load packages ##
###################
library(deSolve)
library(ggplot2)
library(bayestestR)
library(ggpubr)
library(Cairo)

###############################
## Load age profile function ##
###############################
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

################################################
## Read in parameter sets & likelihood values ##
################################################
nsim  <- 600
niter <- 100
npars <- 3

out_likpar <- data.frame(matrix(ncol=npars+1, nrow=niter*nsim))
names(out_likpar) <- c("log.lambda0", "log.shape", "log.tdecline", "likelihood")
counter <- seq(1, (nsim*niter), by = 1/nsim)  # used in for loop to pick correct parliks_ file

#works for multiples of 10
i_seq <- seq(1, niter*nsim, by=nsim)
for(i in i_seq){
  if(i==1){
    out_likpar[i:(i+nsim-1),] <- readRDS(file = paste("mod_output/", pars$country, "/parliks_", pars$country, "_t", pars$temporal_foi, "_a", pars$age_foi, "_", i, ".Rdata", sep = ""))
  } else if (i > 1 & i < (niter*nsim)){
    out_likpar[i:(i+nsim-1),] <- readRDS(file = paste("mod_output/", pars$country, "/parliks_", pars$country, "_t", pars$temporal_foi, "_a", pars$age_foi, "_", i-(i-counter[i]), ".Rdata", sep = ""))
  }
}

colnames(out_likpar) <- c("log.lambda0", "log.shape",  "log.tdecline", "likelihood")

##########################
## Exploratory plotting ##
##########################
## Plot parameter values vs. likelihood
par(mfrow=c(2,2))
plot(exp(out_likpar$log.lambda0), out_likpar$likelihood)
plot(exp(out_likpar$log.shape), out_likpar$likelihood)
plot(exp(out_likpar$log.tdecline), out_likpar$likelihood)

##################################################################################
## SENSE CHECK: simulate model & show fit to data for a range of parameter sets ##
##################################################################################
params <- data.frame(matrix(nrow=3, ncol=npars+1))
names(params) <- c("log.lambda0", "log.shape", "log.tdecline", "likelihood")
sorted_lik    <- out_likpar[order(out_likpar$likelihood),] #sort by likelihood value
sorted_lik <- na.omit(sorted_lik)
params[1,]    <- head(sorted_lik, 1)                #best fit
params[2,]    <- sorted_lik[nrow(out_likpar)/2,]    #middling fit
params[3,]    <- tail(sorted_lik, 1)[1,]            #worst fit

# Create list to store model output for each of the timepoints being fit
ct_cases          <- vector("list", length = nrow(fitting_data))
matched_prev      <- vector("list", length = nrow(fitting_data))
mean_matched_prev <- vector("list", length = nrow(fitting_data))

for(i in 1:nrow(params)){
  
  pars$log.lambda0  <- params[i, "log.lambda0"]
  pars$log.shape    <- params[i, "log.shape"]
  pars$log.tdecline <- params[i, "log.tdecline"]
  sol <- ode(y = y, times = time, parms = pars,  func = age_si)  #save model solution
  
  for(j in 1:nrow(fitting_data)){
    
    store_sim          <- getit(year_diff[j])  #store age profile at jth timepoint
    ct_cases[[i]][j]   <- sum(store_sim[,"ct1"] + store_sim[,"ct2"] + store_sim[,"ct3"])
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
    mean_matched_prev[[i]][j] <- weighted.mean(x = matched_prev[[j]], w = pars$propfert[matched_indices] * pars$Na[matched_indices])
    
    
    }
}

# plot the fits
par(mfrow=c(1,1))
plot(fitting_data$year, fitting_data$prev, ylab = "Seroprevalence", xlab = "Year", ylim = c(0,1))

points(fitting_data$year, mean_matched_prev[[1]], pch=10)
points(fitting_data$year, mean_matched_prev[[2]], col ="blue", pch=10)
points(fitting_data$year, mean_matched_prev[[3]], col = "red", pch=10)


###########################################
## Summarise the posterior distributions ##
###########################################

#priors
prior.lambda   <- out_likpar$log.lambda0
prior.shape    <- out_likpar$log.shape
prior.tdecline <- out_likpar$log.tdecline

#posteriors
sorted_lik    <- head(out_likpar[order(out_likpar$likelihood),], nrow(out_likpar)*0.01) #take top 1%
post.lambda   <- sorted_lik$log.lambda0
post.shape    <- sorted_lik$log.shape
post.tdecline <- sorted_lik$log.tdecline


## Medians & 95% CIs
# lambda0
x<-describe_posterior(post.lambda, centrality = "median")
post.median.lambda <- exp(x$Median)
ci.low.lambda    <- exp(x$CI_low)
ci.high.lambda   <- exp(x$CI_high)

# shape
x<-describe_posterior(post.shape, centrality = "median")
post.median.shape <- exp(x$Median)
ci.low.shape    <- exp(x$CI_low)
ci.high.shape   <- exp(x$CI_high)

# tdecline
x<-describe_posterior(post.tdecline, centrality = "median")
post.median.tdecline <- exp(x$Median)
ci.low.tdecline    <- exp(x$CI_low)
ci.high.tdecline   <- exp(x$CI_high)


############################################
## Plot prior and posterior distributions ##
############################################

## Lambda0

#dataframe
dat <- data.frame(dens = c(prior.lambda, post.lambda), 
                  ci_low = ci.low.lambda,
                  ci_high = ci.high.lambda,
                  lines = c(rep("a", length(prior.lambda)), rep("b", length(post.lambda))))

#plot
lambda0 <- ggplot(dat, aes(x = exp(dens), fill = lines)) + 
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = post.median.lambda) +                   # mode
  geom_vline(xintercept = ci.low.lambda, linetype='dotted') +   # lower ci
  geom_vline(xintercept = ci.high.lambda, linetype='dotted') +  # upper ci
  theme_light(base_size = 12, base_line_size = 0, base_family = "Times") + 
  theme(legend.position = "none", panel.grid = element_blank(), axis.ticks = element_blank()) +
  labs(fill = "") + 
  xlab(expression(lambda)) + 
  ylab("Density") + 
  theme(plot.margin=unit(c(rep(1,4)),"cm"))


## shape

#dataframe
dat <- data.frame(dens = c(prior.shape, post.shape), 
                  ci_low = ci.low.shape,
                  ci_high = ci.high.shape,
                  lines = c(rep("a", length(prior.shape)), rep("b", length(post.shape))))


#plot
shape <- ggplot(dat, aes(x = exp(dens), fill = lines)) + geom_density(alpha = 0.5) + 
  geom_vline(xintercept = post.median.shape) +                   # mode
  geom_vline(xintercept = ci.low.shape, linetype='dotted') +   # lower ci
  geom_vline(xintercept = ci.high.shape, linetype='dotted') +  # upper ci
  theme_light(base_size = 12, base_line_size = 0, base_family = "Times") + 
  theme(legend.position = "none", panel.grid = element_blank(), axis.ticks = element_blank()) +
  labs(fill = "") + 
  # scale_x_continuous(breaks = seq(0, 0.8, 0.4), limits = c(0, 0.8), expand = c(0, 0)) +
  # scale_y_continuous(breaks = seq(0, 1.6, 0.75), limits = c(0, 1.6), expand = c(0, 0)) +
  # xlab("Shape") + 
  xlab(expression(mu)) + 
  ylab("Density") + 
  theme(plot.margin=unit(c(rep(1,4)),"cm"))


## tdecline

#dataframe
dat <- data.frame(dens = c(prior.tdecline, post.tdecline), 
                  ci_low = ci.low.tdecline,
                  ci_high = ci.high.tdecline,
                  lines = c(rep("a", length(prior.tdecline)), rep("b", length(post.tdecline))))

#plot
tdecline <-ggplot(dat, aes(x = exp(dens), fill = lines)) + geom_density(alpha = 0.5) + 
  geom_vline(xintercept = post.median.tdecline) +                   # mode
  geom_vline(xintercept = ci.low.tdecline, linetype='dotted') +   # lower ci
  geom_vline(xintercept = ci.high.tdecline, linetype='dotted') +  # upper ci 
  theme_light(base_size = 12, base_line_size = 0, base_family = "Times") + 
  theme(legend.position = "none", panel.grid = element_blank(), axis.ticks = element_blank()) +
  labs(fill = "") + 
  xlab(expression(tau)) + 
  ylab("Density") + 
  theme(plot.margin=unit(c(rep(1,4)),"cm"))

## Multipanel plot of priors vs. posteriors
ggarrange(lambda0, shape, tdecline, ncol=2, nrow=2, 
          labels=c("(a)", "(b)", "(c)"), font.label=list(size=12, family="Times"), hjust = -2)

# ggsave(filename = paste("plots/", pars$country, "_", pars$temporal_foi, "_posteriors.pdf", sep=""), width = 6, height = 6,
       # units = "in", family = "Times")


####################################
#### Summary of the posteriors #####
####################################
posteriors <- data.frame(par = c("lambda0", "shape", "tdecline"), 
                         med = c(post.median.lambda, post.median.shape, post.median.tdecline),
                         lower_ci = c(ci.low.lambda, ci.low.shape, ci.low.tdecline), 
                         upper_ci = c(ci.high.lambda, ci.high.shape, ci.high.tdecline)
)


###############################
#### Posterior predictions ####
###############################

## Randomly sample n times from posterior & store
n_samples <- 3
lambda.vec   <- sample(post.lambda, n_samples)
shape.vec    <- sample(post.shape, n_samples)
tdecline.vec <- sample(post.tdecline, n_samples)

## Run model for each parameter set, store estimated prevalence & CT incidence
ct_cases <- mean_matched_prev <- vector("list", length = n_samples)
matched_prev <- vector("list", length=max(time))


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


## Calculate median, upper & lower percentiles of model estimates

# Store estimates in matrix (each row = estimates from a different parameter set)
ct_mat   <- t(sapply(ct_cases,c))
prev_mat <- t(sapply(mean_matched_prev,c))


# Column-wise medians & percentiles
ct_med <- ct_lower <- ct_upper <- prev_med <- prev_lower <- prev_upper <- vector("numeric", length=max(time))

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

# Check concordance
plot(true_prev$prev, tp)

## Calculate data 95% CIs
cis <- binom.confint(x=true_prev$k, n=true_prev$n, conf.level=0.95, methods="exact")

# Store lower and upper estimates
true_prev$ci_lo <- cis$lower
true_prev$ci_up <- cis$upper

#####################################
## Plot modelled estimates vs data ##
#####################################

################
## Prevalence ##
################

## Prevalence across years with data ##

# Create df
dat <- data.frame(
  time = fitting_data$year,             #sampling year
  dat_prev = true_prev$prev,            #data adjusted to true prevalence
  dat_low = true_prev$ci_lo,            #data lower interval
  dat_up = true_prev$ci_up,             #data upper interval
  mod_prev = prev_med[year_diff],       #modelled median prevalence
  mod_prev_low = prev_lower[year_diff], #modelled prevalence lower interval
  mod_prev_up = prev_upper[year_diff]   #modelled prevalence upper interval
)

# Plot
prev_fittingyears <- ggplot(data=dat, aes(x=time, y=mod_prev)) + 
  geom_line() + 
  geom_ribbon(aes(ymin = mod_prev_low, ymax = mod_prev_up), alpha=0.2) + 
  geom_point(aes(y=dat_prev), col="grey", size=.3) + 
  geom_errorbar(aes(ymin = dat_low, ymax = dat_up),  col="grey", size=.25) + 
  xlab("Year") + 
  ylab("Seroprevalence") + 
  scale_x_continuous(expand = c(0,0)) + 
  scale_y_continuous(limits=c(0,1), expand = c(0,0)) + 
  theme(plot.margin=unit(c(rep(.5,4)),"cm")) + 
  theme_light(base_size = 12, base_line_size = 0, base_family = "Times") + 
  theme(legend.position = "none")
  
# Save
ggsave(filename = paste("plots/", pars$country, "/", pars$temporal_foi, "_prev_fittingyears", ".pdf", sep=""), 
       device = cairo_pdf, height = 6, width = 6, units = "in")


## Prevalence across years of foi decrease (& forecasting) ##

# time points across which foi is declining
timepoints <- (pars$burnin - round(exp(pars$log.tdecline), 0)) : max(time)

# equivalent years
years <- (min(fitting_data$year) - round(exp(pars$log.tdecline), 0)) : max(fitting_data$year)


# Create df
dat <- data.frame(
  time = years,  #years of foi decrease
  dat_prev = c(rep(NA, (length(years) - length(true_prev$prev))), true_prev$prev), #data adjusted to true prevalence
  dat_low = c(rep(NA, (length(years) - length(true_prev$prev))), true_prev$ci_lo), #data lower interval
  dat_up = c(rep(NA, (length(years) - length(true_prev$prev))), true_prev$ci_up),  #data upper interval
  mod_prev = prev_med[timepoints],       #modelled median prevalence
  mod_prev_low = prev_lower[timepoints], #modelled prevalence lower interval
  mod_prev_up = prev_upper[timepoints]  #modelled prevalence upper interval
)

# plot
prev_allyears <- ggplot(data=dat, aes(x=time, y=mod_prev)) + 
  geom_line() + 
  geom_ribbon(aes(ymin = mod_prev_low, ymax = mod_prev_up), alpha=0.2) + 
  geom_point(aes(y=dat_prev), col="grey", size=.25) + 
  geom_errorbar(aes(ymin = dat_low, ymax = dat_up),  col="grey", size=.2) + 
  xlab("Year") + 
  ylab("Seroprevalence") + 
  scale_x_continuous(expand = c(0,0)) + 
  scale_y_continuous(limits=c(0,1), expand = c(0,0)) + 
  theme(plot.margin=unit(c(rep(.5,4)),"cm")) + 
  theme_light(base_size = 12, base_line_size = 0, base_family = "Times") + 
  theme(legend.position = "none")

# Save
ggsave(filename = paste("plots/", pars$country, "/", pars$temporal_foi, "_prev_allyears", ".pdf", sep=""), 
       device = cairo_pdf, height = 6, width = 6, units = "in")

##################
## CT incidence ##
##################

## CT incidence across years with data ##

# Make CT incidence relative to no. births
fitting_timepoints <- pars$burnin : (pars$burnin + pars$tdiff)
fitting_years <- min(fitting_data$year) : max(fitting_data$year)

births <- sum(pars$propfert*pars$Na)       #annual no. births
ct_rel     <- ct_med   / (births / 10000)  #CT incidence per 10,000 live births
ct_rel_low <- ct_lower / (births / 10000)
ct_rel_up  <- ct_upper / (births / 10000)

# Create df
dat <- data.frame(
  time = fitting_years,                         #sampling year
  ct_rel = ct_rel[fitting_timepoints],          #modelled ct incidence
  ct_rel_low = ct_rel_low[fitting_timepoints],  #modelled ct incidence lower interval
  ct_rel_up = ct_rel_up[fitting_timepoints]     #modelled ct incidence upper interval
)

# Plot
ct_fittingyears <- ggplot(data=dat, aes(x=time, y=ct_rel)) + 
  geom_line() + 
  geom_ribbon(aes(ymin = ct_rel_low, ymax = ct_rel_up), alpha=0.2) + 
  xlab("Year") + 
  ylab("Incidence per 10 000 live births") + 
  scale_x_continuous(expand = c(0,0)) + 
  scale_y_continuous(expand = c(0,0)) + 
  theme(plot.margin=unit(c(rep(.5,4)),"cm")) + 
  theme_light(base_size = 12, base_line_size = 0, base_family = "Times") + 
  theme(legend.position = "none")

# Save
ggsave(filename = paste("plots/", pars$country, "/", pars$temporal_foi, "_CT_fittingyears", ".pdf", sep=""), 
       device = cairo_pdf, height = 6, width = 6, units = "in")


## CT incidence across years with foi decline (& forecasting) ##

# Create df
dat <- data.frame(
  time = years,                    #sampling year
  ct_rel = ct_rel[timepoints],          #modelled ct incidence
  ct_rel_low = ct_rel_low[timepoints],  #modelled ct incidence lower interval
  ct_rel_up = ct_rel_up[timepoints]     #modelled ct incidence upper interval
)

# Plot
ct_allyears <- ggplot(data=dat, aes(x=time, y=ct_rel)) + 
  geom_line() + 
  geom_ribbon(aes(ymin = ct_rel_low, ymax = ct_rel_up), alpha=0.2) + 
  xlab("Year") + 
  ylab("Incidence per 10 000 live births") + 
  scale_x_continuous(expand = c(0,0)) + 
  scale_y_continuous(expand = c(0,0)) + 
  theme(plot.margin=unit(c(rep(.5,4)),"cm")) + 
  theme_light(base_size = 12, base_line_size = 0, base_family = "Times") + 
  theme(legend.position = "none")

# Save
ggsave(filename = paste("plots/", pars$country, "/", pars$temporal_foi, "_CT_allyears", ".pdf", sep=""), 
       device = cairo_pdf, height = 6, width = 6, units = "in")


######################################################
## Multipanel plot for all results from one country ##
######################################################
ggarrange(prev_fittingyears, prev_allyears, ct_fittingyears, ct_allyears, ncol=2, nrow=2, 
          labels=c("(a)", "(b)", "(c)", "(d)"), font.label=list(size=12, family="Times"), 
          hjust = -.2)

# Save
ggsave(filename = paste("plots/", pars$country, "/", pars$temporal_foi, "_multipanel", ".pdf", sep=""), 
       device = cairo_pdf, height = 6, width = 6, units = "in")
