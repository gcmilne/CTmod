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

##  Read in parameter sets & likelihood values 
nsim  <- 10
niter <- 10

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


##########################
## Exploratory plotting ##
##########################
## Plot parameter values vs. likelihood
par(mfrow=c(2,2))
plot(exp(out_likpar$log.lambda0), out_likpar$likelihood)
abline(lm(out_likpar$likelihood ~ exp(out_likpar$log.lambda0)), col="red")

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
top_pars <- head(out_likpar[order(out_likpar$likelihood),], 10)

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
names(params) <- c("log.lambda0", "log.shape", "log.tdecline", "likelihood")
sorted_lik    <- head(out_likpar[order(out_likpar$likelihood),], 100) #sort by likelihood value
params[1,]    <- head(sorted_lik, 1)                #best fit
params[3,]    <- tail(sorted_lik, 3)[1,]            #worst fit
params[2,]    <- sorted_lik[nrow(out_likpar)/2,]    #middling fit

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

par(mfrow=c(1,1))
plot(fitting_data$year, fitting_data$prev, ylab = "Seroprevalence", xlab = "Year", ylim = c(0,1))
for(i in 1:nrow(params)){
  points(fitting_data$year, mean_matched_prev[[i]], col = myColors[i])
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

#calculate number of CT cases in each trimester & time point
## t1
ct1 <- ct2 <- ct3 <- vector("numeric", length=3)
for(i in 1:3){
  ct1[i] <- round(sum(df[[i]][,"ct1"]))
  ct2[i] <- round(sum(df[[i]][,"ct2"]))
  ct3[i] <- round(sum(df[[i]][,"ct3"]))
}
ct_tot <- ct1+ct2+ct3

## t2
ct1 <- ct2 <- ct3 <- vector("numeric", length=3)
for(i in 1:3){
  ct1[i] <- round(sum(df2[[i]][,"ct1"]))
  ct2[i] <- round(sum(df2[[i]][,"ct2"]))
  ct3[i] <- round(sum(df2[[i]][,"ct3"]))
}
ct_tot2 <- ct1+ct2+ct3

ct_tot2/ct_tot
########
### trying out a linear decrease in foi over time
########
top_par    <- head(out_likpar[order(out_likpar$likelihood),], 1) #sort by likelihood value

pars$log.lambda0 <- top_par$log.lambda0
pars$log.lambda1 <- top_par$log.lambda1
pars$log.gradient <- top_par$log.gradient
pars$log.shape <- top_par$log.shape
pars$log.tdecline <- top_par$log.tdecline

sol <- ode(y = y, times = time, parms = pars,  func = age_si)  #save model solution

df   <- getit(pars$burnin)
matched_prev  <- df[,"obs_pI"][matched_indices]
df2  <- getit(pars$burnin+11)
matched_prev2  <- df2[,"obs_pI"][matched_indices]

logliks     <- loglik(k1 = neth_95$k, n1 = neth_95$n, prev1 = matched_prev, 
                        k2 = neth_06$k, n2 = neth_06$n, prev2 = matched_prev2)
sum(-logliks)

par(mfrow=c(2,1))
plot(pars$age[matched_indices], prev_list1, type='l', ylim=c(0,1), ylab="Prevalence",
     xlab="Age (years)", main = "Netherlands 1995")
points(neth_95$age_mid, neth_95$prevalence)

# plot second time point vs. data
plot(pars$age[matched_indices], prev_list2, type='l', ylim=c(0,1), ylab="Prevalence",
     xlab="Age (years)", main = "Netherlands 2006")
points(neth_06$age_mid, neth_06$prevalence)

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

sol[c((pars$burnin),(pars$burnin+11)), "ctt"] ## CT cases over fitting time period

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

############################################
## Plot prior and posterior distributions ##
############################################
library(ggplot2)
## lambda0
#prior
prior <- out_likpar$log.lambda0

#posterior
sorted_lik <- head(out_likpar[order(out_likpar$likelihood),], 100)
post <- sorted_lik$log.lambda0

# Calculate 95% posterior credible interval
# ci_95 <- quantile(exp(post), probs = c(0.025, 0.5, 0.975))

#dataframe
dat <- data.frame(dens = c(prior, post), 
                  lines = c(rep("a", length(prior)), rep("b", length(post))))

#plot
lambda0 <- ggplot(dat, aes(x = exp(dens), fill = lines)) + geom_density(alpha = 0.5) + 
  # geom_vline(xintercept = ci_95, linetype='dotted') +
  theme_light(base_size = 12, base_line_size = 0, base_family = "Times") + 
  theme(legend.position = "none", panel.grid = element_blank(), axis.ticks = element_blank()) +
  labs(fill = "") + 
  scale_x_continuous(breaks = seq(0, 0.2, 0.1), limits = c(0, 0.2), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 150, 75), limits = c(0, 150), expand = c(0, 0)) +
  xlab("Lambda0") + ylab("Density") + 
  theme(plot.margin=unit(c(rep(1,4)),"cm"))

## lambda1
#prior
prior <- out_likpar$log.lambda1

#posterior
post <- sorted_lik$log.lambda1

# Calculate 95% posterior credible interval
# ci_95 <- quantile(exp(post), probs = c(0.025, 0.5, 0.975))

#dataframe
dat <- data.frame(dens = c(prior, post), 
                  lines = c(rep("a", length(prior)), rep("b", length(post))))

#plot
lambda1 <- ggplot(dat, aes(x = exp(dens), fill = lines)) + geom_density(alpha = 0.5) + 
  # geom_vline(xintercept = ci_95, linetype='dotted') +
  theme_light(base_size = 12, base_line_size = 0, base_family = "Times") + 
  theme(legend.position = "none", panel.grid = element_blank(), axis.ticks = element_blank()) +
  labs(fill = "") + 
  scale_x_continuous(breaks = seq(0, 0.2, 0.1), limits = c(0, 0.2), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 7, 3), limits = c(0, 7), expand = c(0, 0)) +
  xlab("Lambda1") + ylab("Density") + 
  theme(plot.margin=unit(c(rep(1,4)),"cm"))

## gradient
#prior
prior <- out_likpar$log.gradient

#posterior
post <- sorted_lik$log.gradient

# Calculate 95% posterior credible interval
# ci_95 <- quantile(exp(post), probs = c(0.025, 0.5, 0.975))

#dataframe
dat <- data.frame(dens = c(prior, post), 
                  lines = c(rep("a", length(prior)), rep("b", length(post))))

#plot
gradient <-
  ggplot(dat, aes(x = exp(dens), fill = lines)) + geom_density(alpha = 0.5) + 
  # geom_vline(xintercept = ci_95, linetype='dotted') +
  theme_light(base_size = 12, base_line_size = 0, base_family = "Times") + 
  theme(legend.position = "none", panel.grid = element_blank(), axis.ticks = element_blank()) +
  labs(fill = "") + 
  scale_x_continuous(breaks = seq(0, 2, 1), limits = c(0, 2), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 1, .5), limits = c(0, 1), expand = c(0, 0)) +
  xlab("Gradient") + ylab("Density") + 
  theme(plot.margin=unit(c(rep(1,4)),"cm"))

## shape
#prior
prior <- out_likpar$log.shape

#posterior
post <- sorted_lik$log.shape

# Calculate 95% posterior credible interval
# ci_95 <- quantile(exp(post), probs = c(0.025, 0.5, 0.975))

#dataframe
dat <- data.frame(dens = c(prior, post), 
                  lines = c(rep("a", length(prior)), rep("b", length(post))))

#plot
shape <- ggplot(dat, aes(x = exp(dens), fill = lines)) + geom_density(alpha = 0.5) + 
  # geom_vline(xintercept = ci_95, linetype='dotted') +
  theme_light(base_size = 12, base_line_size = 0, base_family = "Times") + 
  theme(legend.position = "none", panel.grid = element_blank(), axis.ticks = element_blank()) +
  labs(fill = "") + 
  scale_x_continuous(breaks = seq(0, 0.8, 0.4), limits = c(0, 0.8), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 1.6, 0.75), limits = c(0, 1.6), expand = c(0, 0)) +
  xlab("Shape") + ylab("Density") + 
  theme(plot.margin=unit(c(rep(1,4)),"cm"))

## tdecline
#prior
prior <- out_likpar$log.tdecline

#posterior
post <- sorted_lik$log.tdecline

# Calculate 95% posterior credible interval
# ci_95 <- quantile(exp(post), probs = c(0.025, 0.5, 0.975))

#dataframe
dat <- data.frame(dens = c(prior, post), 
                  lines = c(rep("a", length(prior)), rep("b", length(post))))

#plot
tdecline <-ggplot(dat, aes(x = exp(dens), fill = lines)) + geom_density(alpha = 0.5) + 
  # geom_vline(xintercept = ci_95, linetype='dotted') +
  theme_light(base_size = 12, base_line_size = 0, base_family = "Times") + 
  theme(legend.position = "none", panel.grid = element_blank(), axis.ticks = element_blank()) +
  labs(fill = "") + 
  scale_x_continuous(breaks = seq(0, 100, 50), limits = c(0, 100), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 0.02, 0.01), limits = c(0, 0.02), expand = c(0, 0)) +
  xlab("tdecline") + ylab("Density") + 
  theme(plot.margin=unit(c(rep(1,4)),"cm"))

## Multipanel plot of priors vs. posteriors
library(ggpubr)

png("plots/posteriors_linear_100.png", height=4, width = 8, units="in", res=300)
ggarrange(lambda0, lambda1, gradient, shape, tdecline, ncol=3, nrow=2, 
          labels=c("(a)", "(b)", "(c)", "(d)", "(e)"), 
          font.label=list(size=12, family="Times"), hjust = -2)
dev.off()

#########################################################
##  Calculate 95% CIs for data (Agresti-Coull method)  ##
#########################################################
library(binom)

# '95 timepoint
cis <- binom.confint(x=neth_95$k, n=neth_95$n, conf.level=0.95, methods="exact")
lower_95 <- cis$lower
upper_95 <- cis$upper

# '06 timepoint
cis <- binom.confint(x=neth_06$k, n=neth_06$n, conf.level=0.95, methods="exact")
lower_06 <- cis$lower
upper_06 <- cis$upper

rm(cis)
##########################
## GGPLOTS OF MODEL FIT ##
##########################
library("ggplot2")

sorted_lik <- head(out_likpar[order(out_likpar$likelihood),], 1)

pars$log.lambda0  <- sorted_lik$log.lambda0
pars$log.lambda1  <- sorted_lik$log.lambda1
pars$log.gradient <- sorted_lik$log.gradient
pars$log.shape    <- sorted_lik$log.shape
pars$log.tdecline <- sorted_lik$log.tdecline

sol <- ode(y = y, times = time, parms = pars,  func = age_si)  #save model solution
df  <- getit(pars$burnin)
prev_list1  <- df[,"obs_pI"][matched_indices]  #select observed prevalence from relevant age categories
df2 <- getit(pars$burnin+11)
prev_list2  <- df2[,"obs_pI"][matched_indices]  #select observed prevalence from relevant age categories

#dataset for plotting '95 fit
fit95 <- setNames(data.frame(matrix(ncol = 5, nrow = length(matched_indices))), c("age", "dat", "fit"))
fit95$age <- neth_95$age_mid
fit95$dat <- neth_95$prevalence
fit95$dat_ci_low  <- lower_95
fit95$dat_ci_high <- upper_95
fit95$fit <- prev_list1

#dataset for plotting '06 fit
fit06 <- setNames(data.frame(matrix(ncol = 3, nrow = length(matched_indices))), c("age", "dat", "fit"))
fit06$age <- neth_06$age_mid
fit06$dat <- neth_06$prevalence
fit06$dat_ci_low  <- lower_06
fit06$dat_ci_high <- upper_06
fit06$fit <- prev_list2

#plot '95 fit
t1 <- ggplot(data=fit95, aes(x=age, y=fit)) + 
  geom_line(aes(y=fit)) +
  geom_point(aes(y=dat)) + 
  geom_errorbar(aes(ymin = dat_ci_low, ymax = dat_ci_high, width=.3)) + 
  theme_light(base_size = 12, base_line_size = 0, base_family = "Times") + 
  theme(legend.position = "none", panel.grid = element_blank(), axis.ticks = element_blank()) +
  scale_x_continuous(breaks = seq(0, 80, 20), limits = c(0, 80), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1), expand = c(0, 0)) +
  xlab("Age (years)") + ylab("Seroprevalence") + 
  theme(plot.margin=unit(c(rep(1,4)),"cm"))
t1

#plot '06 fit
t2 <- ggplot(data=fit06, aes(x=age, y=fit)) + 
  geom_line(aes(y=fit)) +
  geom_point(aes(y=dat)) + 
  geom_errorbar(aes(ymin = dat_ci_low, ymax = dat_ci_high, width=.3)) + 
  theme_light(base_size = 12, base_line_size = 0, base_family = "Times") + 
  theme(legend.position = "none", panel.grid = element_blank(), axis.ticks = element_blank()) +
  scale_x_continuous(breaks = seq(0, 80, 20), limits = c(0, 80), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1), expand = c(0, 0)) +
  xlab("Age (years)") + ylab("Seroprevalence") + 
  theme(plot.margin=unit(c(rep(1,4)),"cm"))

## Multipanel plot of observed vs. modelled prevalence
library(ggpubr)

ggarrange(t1, t2, ncol=1, nrow=2, 
          labels=c("(a)", "(b)"), 
          font.label=list(size=12, family="Times"), hjust = -1)

ggsave(filename = "plots/bfit_linear_netherlands.eps", width = 6, height = 6, 
       units = "in", family = "Times")

foi <- exp(pars$log.lambda0) + exp(pars$log.lambda1) * (pars$age * exp(-exp(pars$log.gradient)*pars$age))
plot(pars$age, foi, 'l')

#############################################################
## Calculate credible intervals around parameter estimates ##
#############################################################
lambda0  <- exp(pars$log.lambda0)
lambda1  <- 0.1
gradient <- 0
shape    <- exp(pars$log.shape)
tdecline <- exp(pars$log.tdecline)
foi <- (lambda0 + lambda1 * (pars$age * exp(-gradient*pars$age)))*shape #decrease foi after burnin
par(mfrow=c(1,1))
plot(pars$age, foi, 'l')

#############################
### Supplementary figures ###
#############################
sol <- ode(y = y, times = time, parms = pars,  func = age_si)  #save model solution
df <- getit(pars$burnin)

# Plot N over age
data <- data.frame(mod_a=df[,"a"], mod_Na=df[,"Na"], dat_a=head(pars$age, -1), dat_Na=head(pars$Na, -1))

pop_size <- ggplot(data=data, aes(x=mod_a, y=Nt)) + 
  geom_line(aes(x=mod_a, y=mod_Na), size=1) +
  geom_point(aes(x=dat_a, y=dat_Na), colour = "grey", size=1) + 
  theme_light(base_size = 12, base_line_size = 0, base_family = "Times") + 
  theme(legend.position = "none", panel.grid = element_blank(), axis.ticks = element_blank()) +
  # scale_x_continuous(breaks = seq(0, 80, 20), limits = c(0, 80), expand = c(0, 0)) +
  # scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1), expand = c(0, 0)) +
  xlab("Age (years)") + ylab("Population size") + 
  theme(plot.margin=unit(c(rep(1,4)),"cm"))
pop_size
# ggsave(filename = "plots/pop_size.eps", width = 6, height = 6, 
       # units = "in", family = "Times")

# Plot burn-in


###### Plot entire posterior of foi over age
sorted_lik <- exp(head(out_likpar[order(out_likpar$likelihood),], 100))

foi <- sorted_lik$log.lambda0[1] + sorted_lik$log.lambda1[1] * (pars$age * exp(-sorted_lik$log.gradient[1]*pars$age))
plot(pars$age, foi, 'l', ylim=c(0,0.1))

foi <- rep(list(vector(length=pars$agrps)), 100)
for(i in 1:100){
  foi[[i]] <- sorted_lik$log.lambda0[i] + sorted_lik$log.lambda1[i] * (pars$age * exp(-sorted_lik$log.gradient[i]*pars$age))
  lines(pars$age, foi[[i]])
}


