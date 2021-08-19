############################################################
## Script to plot posterior predictions from each country ##
############################################################

#################
# Load packages #
#################
library(ggplot2)
library(binom)
library(patchwork)
library(RColorBrewer)

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

births <- sum(pars$propfert*pars$Na)       #annual no. births
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


################
## Make plots ##
################

## Prevalence ##

# make list to store all plots from different countries
if (!exists("prev_allyears")) {  #only create if list not in existence
  prev_allyears <- vector("list", length=length(countries))
}

#Colour for geom_ribbon
# display.brewer.pal(8, "Pastel1")
ribbonColour <- brewer.pal(8, "Pastel2")[c(3, 7)]

prev_allyears[[which(countries == pars$country)]] <- ggplot(
  data=prev_all[[which(countries == pars$country)]]) + 
  geom_line(aes(x=time, y=mod_prev), size=0.3) + 
  geom_ribbon(aes(x=time, ymin = mod_prev_low, ymax = mod_prev_up), alpha=0.5, fill=ribbonColour[2]) +
  annotate("rect", xmin=min(fitting_data$year), xmax=max(fitting_data$year), #years with data
           ymin=0, ymax=Inf, alpha=0.3, fill=ribbonColour[1]) +
  xlab("Year") + 
  ylab("Seroprevalence") + 
  scale_x_continuous(expand = c(0,0), limits = c(1900, 2030), n.breaks = 3) + 
  scale_y_continuous(expand = c(0,0)) +
  theme(plot.margin=unit(c(rep(.5,4)),"cm")) + 
  theme_light(base_size = 12, base_line_size = 0, base_family = "Times") + 
  theme(legend.position = "none") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

## Calculate where to place inset within main plot area
space_below_plot <- min(prev_all[[which(countries == pars$country)]]$mod_prev_low[which(years == 1900):which(years == 1975)])
space_above_plot <- 1 - (max(prev_all[[which(countries == pars$country)]]$mod_prev_up[which(years == 1900):which(years == 1975)]))

if (space_below_plot > space_above_plot) { 
  
  ymin_val <- 0.05
  ymax_val <- space_below_plot * 0.92
  xmin_val <- 1900 + 10
  xmax_val <- min(fitting_data$year) - 10
  
  # make arrow
  arrow_vec <- data.frame(
    x1 = min(fitting_data$year),
    x2 = xmax_val, 
    y1 = prev_all[[which(countries == pars$country)]]$mod_prev_low[which(prev_all[[which(countries == pars$country)]]$time == min(fitting_data$year))], #find intersection
    y2 = ymax_val * 0.95
  )
  
} else if (space_below_plot < space_above_plot) { 
  
  ymin_val <- 1 - space_above_plot + (0.02 * space_above_plot)
  ymax_val <- 0.98
  xmin_val <- 1900 + 10
  xmax_val <- min(fitting_data$year) - 10
  
  # make arrow
  arrow_vec <- data.frame(
    x1 = min(fitting_data$year),
    x2 = xmax_val, 
    y1 = prev_all[[which(countries == pars$country)]]$mod_prev_up[which(prev_all[[which(countries == pars$country)]]$time == min(fitting_data$year))], #find intersection
    y2 = ymin_val
  )
  
  prev_allyears[[which(countries == pars$country)]] <- prev_allyears[[which(countries == pars$country)]] +
    scale_y_continuous(expand = c(0,0), limits=c(0,1))  # change main plot limits to allow inset to fit
  
}

## Prevalence inset graph
if (!exists("prev_inset")) {  #only create if list not in existence
  prev_inset <- vector("list", length=length(countries))
}

# Calculate no. breaks for y axis
if (ymax_val - ymin_val > 0.4) {
  num_breaks <- 4
} else if (ymax_val - ymin_val < 0.4) { 
  num_breaks <- 3
} 

if (pars$country == "Italy" | pars$country == "Turkey"){
  num_breaks <- 2
}

prev_inset[[which(countries == pars$country)]] <- ggplot(
  data=prev_fit[[which(countries == pars$country)]], aes(x=time, y=mod_prev)) + 
  geom_line(size=0.3) + 
  geom_ribbon(aes(ymin = mod_prev_low, ymax = mod_prev_up), alpha=0.5, fill = ribbonColour[2]) +
  geom_point(aes(y=dat_prev), col="grey", size=0.2) +
  geom_errorbar(aes(width=.3, ymin = dat_low, ymax = dat_up),  col="grey", size=0.2) +
  scale_x_continuous(breaks=c(min(fitting_data$year), max(fitting_data$year))) +
  coord_cartesian(xlim=c(min(fitting_data$year), max(fitting_data$year, 
                                                     ylim=c(min(prev_fit[[which(countries == pars$country)]]$dat_low), 
                                                            max(prev_fit[[which(countries == pars$country)]]$dat_up))))) + #avoids error bars getting clipped
  scale_y_continuous(n.breaks=num_breaks) +
  theme_light(base_size = 12, base_line_size = 0, base_family = "Times") + 
  theme(legend.position = "none", axis.title.x=element_blank(), axis.title.y=element_blank()) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme(plot.margin=grid::unit(c(rep(0.1,4)),"cm"))

  
## Combined main prevalence plot with inset
if (!exists("prev_combo")) {  #only create if list not in existence
  prev_combo <- vector("list", length=length(countries))
}

prev_combo[[which(countries == pars$country)]] <- prev_allyears[[which(countries == pars$country)]] + 
  annotation_custom(
    ggplotGrob(prev_inset[[which(countries == pars$country)]]), 
    xmin = xmin_val, xmax = xmax_val, ymin = ymin_val, ymax = ymax_val
  ) + 
  geom_segment(data = arrow_vec, aes(x = x1, y = y1,       #add arrow
                                     xend = x2, yend = y2),
               lineend = "round", 
               linejoin = "round",
               size = 0.3, 
               arrow = arrow(length = unit(0.04, "inches")),
               colour = "grey") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())



## CT incidence ##

# make list to store all plots from different countries
if (!exists("ct_allyears")) {  #only create if list not in existence
  ct_allyears <- vector("list", length=length(countries))
}

# Main CT incidence plot
if(pars$country != "China"){
  ct_allyears[[which(countries == pars$country)]] <- ggplot(
    data=ct_all[[which(countries == pars$country)]], aes(x=time, y=ct_rel)) + 
    geom_line(aes(y=ct_rel), size=.3) +                        # Overall CT cases
    geom_ribbon(aes(ymin = ct_rel_low, ymax = ct_rel_up), alpha=0.5, fill = ribbonColour[2]) +
    annotate("rect", xmin=min(fitting_data$year), xmax=max(fitting_data$year), #years with data
             ymin=0, ymax=Inf, alpha=0.3, fill=ribbonColour[1]) +
    xlab("Year") + 
    ylab("Incidence per 10 000 live births") + 
    scale_x_continuous(expand = c(0,0), limits = c(1900, 2030), n.breaks = 3) + 
    scale_y_continuous(expand = c(0,0), n.breaks = 5) + 
    theme(plot.margin=unit(c(rep(.5,4)),"cm")) + 
    theme_light(base_size = 12, base_line_size = 0, base_family = "Times") + 
    theme(legend.position = "none") + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
} else if(pars$country == "China"){
  ct_allyears[[which(countries == pars$country)]] <- ggplot(
    data=ct_all[[which(countries == pars$country)]], aes(x=time, y=ct_rel)) + 
    geom_line(aes(y=ct_rel), size=.3) +                        # Overall CT cases
    geom_ribbon(aes(ymin = ct_rel_low, ymax = ct_rel_up), alpha=0.5, fill = ribbonColour[2]) +
    annotate("rect", xmin=min(fitting_data$year), xmax=max(fitting_data$year), #years with data
             ymin=0, ymax=Inf, alpha=0.3, fill=ribbonColour[1]) +
    xlab("Year") + 
    ylab("Incidence per 10 000 live births") + 
    scale_x_continuous(expand = c(0,0), limits = c(1900, 2030), n.breaks = 3) + 
    scale_y_continuous(limits=c(0, 300), expand = c(0,0), n.breaks = 5) + 
    theme(plot.margin=unit(c(rep(.5,4)),"cm")) + 
    theme_light(base_size = 12, base_line_size = 0, base_family = "Times") + 
    theme(legend.position = "none") + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
}

# Save
# ggsave(filename = paste("plots/", pars$country, "/", pars$temporal_foi, "_CT_allyears", ".pdf", sep=""),
# device = cairo_pdf, height = 6, width = 6, units = "in")


# CT sequelae incidence inset plot
# make list to store all plots from different countries
if (!exists("ct_inset")) {  #only create if list not in existence
  ct_inset <- vector("list", length=length(countries))
}

# colours
library(RColorBrewer)
# display.brewer.pal(8, "Pastel2") # choose colours

myColors <- brewer.pal(7,"Set2") #set colours to use

inset_line_size = 0.3

max_incidence <- max(map_dbl(ct_all[[which(countries == pars$country)]][,-(1:4)], max)) # find max sequelae incidence

ct_inset[[which(countries == pars$country)]] <- ggplot(
  data=ct_all[[which(countries == pars$country)]], aes(x=time, y=ct_rel)) + 
  geom_line(aes(y=chorioretinitis_first), size=inset_line_size, col = myColors[1]) +         # Chorioretinitis in first year
  geom_line(aes(y=chorioretinitis_later), size=inset_line_size, col = myColors[2]) +         # Chorioretinitis in later life
  geom_line(aes(y=intracranial_calcifications), size=inset_line_size, col = myColors[3]) +   # Intracranial calcifications
  geom_line(aes(y=cns_abnormalities), size=inset_line_size, col = myColors[4]) +             # CNS abnormalities
  geom_line(aes(y=foetal_loss), size=inset_line_size, col = myColors[5]) +                   # Foetal loss
  geom_line(aes(y=hydrocephalus), size=inset_line_size, col = myColors[6]) +                 # Hydrocephalus
  geom_line(aes(y=neonatal_death), size=inset_line_size, col = myColors[7]) +                # Neonatal death
  annotate("rect", xmin=min(fitting_data$year), xmax=max(fitting_data$year), #years with data
           # ymin=min(ct_all[[which(countries == pars$country)]]$ct_rel_low), 
           ymin=0, ymax=Inf, alpha=0.3, fill=ribbonColour[1]) +
  scale_x_continuous(expand = c(0,0), limits = c(1900, 2030), n.breaks=2) + 
  scale_y_continuous(limits=c(0, max_incidence+0.05*max_incidence), expand = c(0,0), n.breaks=3) + 
  theme(plot.margin=unit(c(rep(.5,4)),"cm")) + 
  theme_light(base_size = 12, base_line_size = 0, base_family = "Times") + 
  theme(legend.position = "none", axis.title.x=element_blank(), axis.title.y=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme(plot.margin=grid::unit(c(0,0,0,0),"cm"))


# ct_inset[[which(countries == pars$country)]]


# Combine main plot with inset plot
# make list to store all plots from different countries
if (!exists("ct_combo")) {  #only create if list not in existence
  ct_combo <- vector("list", length=length(countries))
}

if (pars$country == "Cameroon" | pars$country == "Ethiopia") {
  ymin_val <-  max(ct_all[[which(countries == pars$country)]]$ct_rel_up[which(years == 1900):which(years == 1975)]) * 1.1
  ymax_val <-  max(ct_all[[which(countries == pars$country)]]$ct_rel_up[which(years == 1900):which(years == 2025)]) * 0.90
  xmin_val <- 1910
  xmax_val <- 1985
  
} else if (pars$country == "Iran (Islamic Republic of)") { 
  ymin_val <- 0
  ymax_val <- min(ct_all[[which(countries == pars$country)]]$ct_rel_low[which(years == 1940):which(years == 2000)]) * 0.91
  xmin_val <- 1940
  xmax_val <- 2000
  
} else if (pars$country == "China") {
  ymin_val <- 150
  ymax_val <- 300 * 0.91
  xmin_val <- 1910
  xmax_val <- 1985
  
  } else if (pars$country != "Cameroon" | pars$country != "Ethiopia" | 
             pars$country != "Iran (Islamic Republic of)" | pars$country != "China") { 
  ymin_val <- 0
  ymax_val <- min(ct_all[[which(countries == pars$country)]]$ct_rel_low[which(years == 1900):which(years == 1975)]) * 0.90
  xmin_val <- 1910
  xmax_val <- 1985
  
}

## Combined CT incidence plot with sequelae inset plot
ct_combo[[which(countries == pars$country)]] <- ct_allyears[[which(countries == pars$country)]] + 
  annotation_custom(
    ggplotGrob(ct_inset[[which(countries == pars$country)]]), 
    xmin = xmin_val, xmax = xmax_val, ymin = ymin_val, ymax = ymax_val
  ) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

###################### 
## Multipanel plots ##
######################

## Prevalence, with inset
wrap_plots(prev_combo, nrow=4, ncol=3) + 
  plot_annotation(tag_levels = list(c('(a)', '(b)', '(c)', '(d)', '(e)', '(f)', 
                                      '(g)', '(h)', '(i)', '(j)', '(h)')))

# letters[which(pars$country == countries)]  #which letter corresponds to which country?

ggsave(filename = "plots/prev_inset_multipanel.pdf",
       device = cairo_pdf, height = 8, width = 8, units = "in")

## CT, with inset
wrap_plots(ct_combo, nrow=4, ncol=3) & 
  scale_x_continuous(expand = c(0,0), limits = c(1900, 2030), n.breaks = 3) & 
  ylab("Incidence") &  
  plot_annotation(tag_levels = list(c('(a)', '(b)', '(c)', '(d)', '(e)', '(f)', 
                                      '(g)', '(h)', '(i)', '(j)', '(h)')))

ggsave(filename = "plots/ct_sequelae_multipanel.pdf",
       device = cairo_pdf, height = 8, width = 8, units = "in")
