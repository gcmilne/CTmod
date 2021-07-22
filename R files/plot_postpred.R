############################################################
## Script to plot posterior predictions from each country ##
############################################################

# READ IN the posterior distribution
rm(post.lambda, post.shape, post.tdecline)
post <- readRDS(file = paste("posteriors/", pars$country, "/", "posteriors_", pars$country, "_t", pars$temporal_foi, "_", "a", pars$age_foi, ".RDS", sep=""))
attach(post)

###################################
## Load in posterior predictions ##
###################################
n_samples  <- 6    #no. parameter sets run on each cluster job
niter      <- 100  #overall no. cluster jobs

prev_list <- ct_list <- vector("list", length=niter)

for(i in 1:niter){
  
  prev_list[[i]] <- readRDS(file = paste("posteriors/", pars$country, "/", "prev_predictions/", "prev_predictions_", pars$country, "_t", pars$temporal_foi, "_a", 
                                         pars$age_foi, "_", i, ".Rdata", sep = ""))
  
  ct_list[[i]] <- readRDS(file = paste("posteriors/", pars$country, "/", "ct_predictions/", "ct_predictions_", pars$country, "_t", pars$temporal_foi, "_a", 
                                       pars$age_foi, "_", i, ".Rdata", sep = ""))
}


## put data lists into dataframes
prev_mat <- do.call(cbind.data.frame, prev_list)
ct_mat <- do.call(cbind.data.frame, ct_list)

## Calculate median, upper & lower percentiles of model estimates
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
# plot(true_prev$prev, tp)

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
# make list to store all plots from different countries
if (which(countries == pars$country) == 1) {
  prev_fittingyears <- vector("list", length=length(countries))
}

prev_fittingyears[[which(countries == pars$country)]] <- ggplot(data=dat, aes(x=time, y=mod_prev)) + 
  geom_line(size=0.3) + 
  geom_ribbon(aes(ymin = mod_prev_low, ymax = mod_prev_up), alpha=0.2) +
  geom_point(aes(y=dat_prev), col="grey", size=.3) + 
  geom_errorbar(aes(width=.3, ymin = dat_low, ymax = dat_up),  col="grey", size=.25) + 
  xlab("Year") + 
  ylab("Seroprevalence") + 
  scale_x_continuous(expand = c(0,0)) + 
  # scale_y_continuous(limits=c(0.45, 0.85), expand = c(0,0)) +
  theme(plot.margin=unit(c(rep(.5,4)),"cm")) + 
  theme_light(base_size = 12, base_line_size = 0, base_family = "Times") + 
  theme(legend.position = "none")

# Save
# ggsave(filename = paste("plots/", pars$country, "/", pars$temporal_foi, "_prev_fittingyears", ".pdf", sep=""),
# device = cairo_pdf, height = 6, width = 6, units = "in")


## Prevalence across years of foi decrease ##
if (pars$forecast == 0) { #without forecasting
  
  # time points across which foi is declining
  timepoints <- (pars$burnin - round(exp(pars$log.tdecline), 0)) : max(time)
  
  # equivalent years
  years <- (min(fitting_data$year) - round(exp(pars$log.tdecline), 0)) : (max(fitting_data$year))
  
  
  # Create df for model estimates
  mod <- data.frame(
    time = years,  #years of foi decrease
    mod_prev     = prev_med[timepoints],   #modelled median prevalence
    mod_prev_low = prev_lower[timepoints], #modelled prevalence lower interval
    mod_prev_up = prev_upper[timepoints]   #modelled prevalence upper interval
  )
  
  # Create df for data
  dat <- data.frame(
    years = fitting_data$year, 
    dat_prev = true_prev$prev,
    dat_low = true_prev$ci_lo,
    dat_up = true_prev$ci_up
  )
  
  
} else if (pars$forecast == 1){  #with forecasting
  
  # time points across which foi is declining
  timepoints <- (pars$burnin - (round(max(exp(post$post.tdecline)), 0))-10) : max(time)
  
  # equivalent years
  years <- (min(fitting_data$year) - (round(max(exp(post$post.tdecline)), 0))-10) : (max(fitting_data$year)+pars$years_forecast)
  
  # period from beginning of foi decrease to final data timepoint
  sampling_period <- (min(fitting_data$year) - round(exp(pars$log.tdecline), 0)) : max(fitting_data$year)
  
  # Create df for model estimates
  mod <- data.frame(
    time = years,  #years of foi decrease
    mod_prev     = prev_med[timepoints],   #modelled median prevalence
    mod_prev_low = prev_lower[timepoints], #modelled prevalence lower interval
    mod_prev_up = prev_upper[timepoints]   #modelled prevalence upper interval
  )
  
  # Create df for data
  dat <- data.frame(
    years = fitting_data$year, 
    dat_prev = true_prev$prev,
    dat_low = true_prev$ci_lo,
    dat_up = true_prev$ci_up
  )
  
  
}
# plot the two datasets [note the empty call to ggplot()]
# make list to store all plots from different countries
if (which(countries == pars$country) == 1) {
  prev_allyears <- vector("list", length=length(countries))
}

prev_allyears[[which(countries == pars$country)]] <- ggplot() + 
  geom_line(data = mod, aes(x=years, y=mod_prev), size=0.3) + 
  geom_ribbon(data = mod, aes(x=years, ymin = mod_prev_low, ymax = mod_prev_up), alpha=0.2) +
  # geom_point(data = dat, aes(x=years, y=dat_prev), col="grey", size=.6) +
  # geom_errorbar(data = dat, aes(width=.1, x=years, ymin = dat_low, ymax = dat_up),  col="grey", size=.6) +
  xlab("Year") + 
  ylab("Seroprevalence") + 
  scale_x_continuous(expand = c(0,0)) + 
  # scale_y_continuous(limits=c(0.45, 0.85), expand = c(0,0)) +
  theme(plot.margin=unit(c(rep(.5,4)),"cm")) + 
  theme_light(base_size = 12, base_line_size = 0, base_family = "Times") + 
  theme(legend.position = "none")

# Save
# ggsave(filename = paste("plots/", pars$country, "/", pars$temporal_foi, "_prev_allyears", ".pdf", sep=""),
# device = cairo_pdf, height = 6, width = 6, units = "in")

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
# make list to store all plots from different countries
if (which(countries == pars$country) == 1) {
  ct_fittingyears <- vector("list", length=length(countries))
}

ct_fittingyears[[which(countries == pars$country)]] <- ggplot(data=dat, aes(x=time, y=ct_rel)) + 
  geom_line(size=.3) + 
  geom_ribbon(aes(ymin = ct_rel_low, ymax = ct_rel_up), alpha=0.2) + 
  xlab("Year") + 
  ylab("Incidence per 10 000 live births") + 
  scale_x_continuous(expand = c(0,0)) + 
  scale_y_continuous(expand = c(0,0)) + 
  theme(plot.margin=unit(c(rep(.5,4)),"cm")) + 
  theme_light(base_size = 12, base_line_size = 0, base_family = "Times") + 
  theme(legend.position = "none")

# Save
# ggsave(filename = paste("plots/", pars$country, "/", pars$temporal_foi, "_CT_fittingyears", ".pdf", sep=""),
# device = cairo_pdf, height = 6, width = 6, units = "in")


## CT incidence across years with foi decline (& forecasting) ##

if (pars$forecast == 0){
  
  # Create df
  dat <- data.frame(
    time = years,                    #sampling year
    ct_rel = ct_rel[timepoints],          #modelled ct incidence
    ct_rel_low = ct_rel_low[timepoints],  #modelled ct incidence lower interval
    ct_rel_up = ct_rel_up[timepoints]     #modelled ct incidence upper interval
  )
  
  
} else if (pars$forecast == 1){
  
  # time points across which foi is declining
  timepoints <- (pars$burnin - (round(max(exp(post$post.tdecline)), 0))-10) : max(time)
  
  # equivalent years
  years <- (min(fitting_data$year) - (round(max(exp(post$post.tdecline)), 0))-10) : (max(fitting_data$year)+pars$years_forecast)
  
  # Create df for model estimates
  dat <- data.frame(
    time = years,  #years of foi decrease
    ct_rel = ct_rel[timepoints],          #modelled ct incidence
    ct_rel_low = ct_rel_low[timepoints],  #modelled ct incidence lower interval
    ct_rel_up = ct_rel_up[timepoints]     #modelled ct incidence upper interval
  )
  
  
}


# Plot
# make list to store all plots from different countries
if (which(countries == pars$country) == 1) {
  ct_allyears <- vector("list", length=length(countries))
}

ct_allyears[[which(countries == pars$country)]] <- ggplot(data=dat, aes(x=time, y=ct_rel)) + 
  geom_line(size=.3) + 
  geom_ribbon(aes(ymin = ct_rel_low, ymax = ct_rel_up), alpha=0.2) + 
  xlab("Year") + 
  ylab("Incidence per 10 000 live births") + 
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) + 
  theme(plot.margin=unit(c(rep(.5,4)),"cm")) + 
  theme_light(base_size = 12, base_line_size = 0, base_family = "Times") + 
  theme(legend.position = "none")

# Save
# ggsave(filename = paste("plots/", pars$country, "/", pars$temporal_foi, "_CT_allyears", ".pdf", sep=""),
# device = cairo_pdf, height = 6, width = 6, units = "in")


######################################################
## Multipanel plot for all results from one country ##
######################################################
country_index <- which(countries == pars$country)
ggarrange(prev_fittingyears[[country_index]], prev_allyears[[country_index]], 
          ct_fittingyears[[country_index]], ct_allyears[[country_index]], ncol=2, nrow=2,
          labels=c("(a)", "(b)", "(c)", "(d)"), font.label=list(size=12, family="Times"),
          hjust = -.2)


# Save
# ggsave(filename = paste("plots/", pars$country, "/", pars$temporal_foi, "_multipanel", ".pdf", sep=""),
# device = cairo_pdf, height = 6, width = 6, units = "in")


### Multipanel with each modelled result for all countries
## (1) Prevalence (fitting years)
wrap_plots(prev_fittingyears[[1]], prev_fittingyears[[2]])
ggarrange(prev_fittingyears, prev_fittingyears[[2]], prev_fittingyears[[3]], ncol=3, nrow=4, 
          labels=c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)", "(h)", "(i)", "(j)", "(k)"),
          font.label=list(size=12, family="Times"), hjust = -.2)

## (2) Prevalence (all years)
wrap_plots(prev_allyears[[1]], prev_allyears[[2]])

## (3) CT incidence (fitting years)
wrap_plots(ct_fittingyears[[1]], ct_fittingyears[[2]])

## (4) CT incidence (all years)
wrap_plots(ct_allyears[[1]], ct_allyears[[2]])

