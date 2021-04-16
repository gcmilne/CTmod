######################################
## POST-CLUSTER analysis & plotting ##
######################################

##  Read in parameter sets & likelihood values 
out_likpar <- data.frame(matrix(ncol=npars+1, nrow=niter*nsim))
names(out_likpar) <- c("log.lambda0", "log.lambda1", "log.gradient", "log.shape", "log.tdecline", "likelihood")
counter <- seq(1, (nsim*niter), by = 1/nsim)  # used in for loop to pick correct parliks_ file

#works for multiples of 10
i_seq <- seq(1, niter*nsim, by=nsim)
for(i in i_seq){
  if(i==1){
    out_likpar[i:(i+nsim-1),] <- readRDS(file = paste("mod_output/fit4_5pars/parliks5-r2_", i, ".RData", sep = ""))
  } else if (i > 1 & i < (niter*nsim)){
    out_likpar[i:(i+nsim-1),] <- readRDS(file = paste("mod_output/fit4_5pars/parliks5-r2_", i-(i-counter[i]), ".RData", sep = ""))
  }
}

lambda1 <- 0.15
lambda0<-0.01
gradient <- 0.6
foi <- lambda0 + lambda1 * (pars$age * exp(-gradient*pars$age))
plot(pars$age, foi)

##########################
## Exploratory plotting ##
##########################
## Plot parameter values vs. likelihood
par(mfrow=c(2,2))
plot(exp(out_likpar$log.lambda0), out_likpar$likelihood)
abline(lm(out_likpar$likelihood ~ exp(out_likpar$log.lambda0)), col="red")

plot(exp(out_likpar$log.lambda1), out_likpar$likelihood)
abline(lm(out_likpar$likelihood ~ exp(out_likpar$log.lambda1)), col="red")

plot(exp(out_likpar$log.gradient), out_likpar$likelihood)
abline(lm(out_likpar$likelihood ~ exp(out_likpar$log.gradient)), col="red")

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
top_pars <- head(out_likpar[order(out_likpar$likelihood),], 1000)

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
names(params) <- c("log.lambda0", "log.lambda1", "log.gradient", "log.shape", "log.tdecline", "likelihood")
sorted_lik    <- head(out_likpar[order(out_likpar$likelihood),], 1000) #sort by likelihood value
params[1:3,]  <- head(sorted_lik, 3)           #best fit
# params[3,]    <- out_likpar[tail(sorted_lik, 1),]               #worst fit
# params[4,]    <- out_likpar[sorted_lik[nrow(out_likpar)/2],]    #middling fit

prev_list1 <- rep(list(matrix(nrow=length(neth_95$age_mid), ncol=1)), nrow(params)) #for timepoint 1
prev_list2 <- rep(list(matrix(nrow=length(neth_95$age_mid), ncol=1)), nrow(params)) #for timepoint 2

df  <- rep(list(matrix(nrow=pars$agrps-1, ncol=16)), nrow(params)) #for timepoint 1
df2 <- rep(list(matrix(nrow=pars$agrps-1, ncol=16)), nrow(params)) #for timepoint 2

for(i in 1:nrow(params)){
  pars$log.lambda0  <- params[i, "log.lambda0"]
  pars$log.lambda1  <- params[i, "log.lambda1"]
  pars$log.gradient <- params[i, "log.gradient"]
  pars$log.shape    <- params[i, "log.shape"]
  pars$log.tdecline <- params[i, "log.tdecline"]
  sol <- ode(y = y, times = time, parms = pars,  func = age_si)  #save model solution
  df[[i]]  <- getit(pars$burnin)
  prev_list1[[i]]  <- df[[i]][,"obs_pI"][matched_indices]  #select observed prevalence from relevant age categories
  df2[[i]] <- getit(pars$burnin+11)
  prev_list2[[i]]  <- df2[[i]][,"obs_pI"][matched_indices]  #select observed prevalence from relevant age categories
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
pars$log.shape <- log(0.10)
pars$log.tdecline <- log(10)

sol <- ode(y = y, times = time, parms = pars,  func = age_si)  #save model solution
#### find out why time dependence of foi returning NA
plot(sol[,"time"], sol[,"pIt"], type='l')

df   <- getit(pars$burnin)
prev_list1  <- df[,"obs_pI"][matched_indices]
df2  <- getit(pars$burnin+11)
prev_list2  <- df2[,"obs_pI"][matched_indices]

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
sorted_lik <- head(out_likpar[order(out_likpar$likelihood),], 10)
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
  # scale_x_continuous(breaks = seq(0, 0.2, 0.1), limits = c(0, 0.2), expand = c(0, 0)) +
  # scale_y_continuous(breaks = seq(0, 54, 20), limits = c(0, 54), expand = c(0, 0)) +
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
  # scale_x_continuous(breaks = seq(0, 0.1, 0.05), limits = c(0, 0.1), expand = c(0, 0)) +
  # scale_y_continuous(breaks = seq(0, 12, 4), limits = c(0, 12), expand = c(0, 0)) +
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
  # scale_x_continuous(breaks = seq(0, 1, 0.5), limits = c(0, 1), expand = c(0, 0)) +
  # scale_y_continuous(breaks = seq(0, 1.8, .5), limits = c(0, 1.8), expand = c(0, 0)) +
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
  # scale_x_continuous(breaks = seq(0, 0.8, 0.4), limits = c(0, 0.8), expand = c(0, 0)) +
  # scale_y_continuous(breaks = seq(0, 3, 1), limits = c(0, 3), expand = c(0, 0)) +
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
  # scale_x_continuous(breaks = seq(0, 100, 50), limits = c(0, 100), expand = c(0, 0)) +
  # scale_y_continuous(breaks = seq(0, 0.03, 0.01), limits = c(0, 0.03), expand = c(0, 0)) +
  xlab("tdecline") + ylab("Density") + 
  theme(plot.margin=unit(c(rep(1,4)),"cm"))

## Multipanel plot of priors vs. posteriors
library(ggpubr)

png("plots/posteriors_100.png", height=4, width = 8, units="in", res=300)
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

ggsave(filename = "plots/bfit_stepwise_netherlands.eps", width = 6, height = 6, 
       units = "in", family = "Times")


#############################################################
## Calculate credible intervals around parameter estimates ##
#############################################################

