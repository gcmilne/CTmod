library("dplyr")
library("bayesplot")
library("ggplot2")
library("rstan")
###############
#summarise fit#
###############
fit$summary("comp_pI")
fit$summary("lambda0")

##############################
## View draws of parameters ##
##############################
fit$draws("lambda0")
#plot density
mcmc_dens(fit$draws("lambda0"))
#scatter plot one parameter against another
mcmc_scatter(fit$draws(c("pars1", "pars2")), alpha = 0.3)

##############################
#create object containing....#
##############################
##...fitted serprevalence
draws_array <- fit$draws("comp_pI")
str(draws_array)
draws_df <- as_draws_df(draws_array) # as_draws_matrix() for matrix
#plot model estimate vs. seroprevalence data
plot(1:80, draws_df[1,1:80], xlab="age (years)", ylab="seroprevalence", main="1 parameter best-fit model")
points(full_data$age_mid, full_data$k/full_data$n, col="red")

##...Na
draws_array <- fit$draws("comp_Na")
str(draws_array)
draws_df <- as_draws_df(draws_array) # as_draws_matrix() for matrix
#plot model estimate vs. seroprevalence data
plot(1:80, draws_df[1,1:80])
points(full_data$age_mid, full_data$k/full_data$n, col="red")

mcmc_hist(fit$draws("lambda0"), binwidth = 100)

#############################
## Diagnose fitting issues ##
#############################
fit$cmdstan_diagnose()
fit$cmdstan_summary()

#create stanfit object
stanfit <- rstan::read_stan_csv(fit$output_files())

#Plot histogram of two parameters on a grid
bayesplot_grid(
  mcmc_hist(fit$draws("lambda0"), binwidth = 0.025),
  mcmc_hist(fit$draws("lambda0"), binwidth = 0.025),
  titles = c("Posterior distribution from MCMC", ""),
  xlim = c(0, .1)
)

################
## Trace plot ##
################
trace <- stan_trace(fit)
trace +
  scale_color_manual(values = c("red", "blue", "green", "black"))
trace +
  scale_color_brewer(type = "div") +
  theme(legend.position = "none")

############################
## Calculate no. ct cases ##
############################
#define & unlist pI from model fitting
draws_array <- fit$draws("comp_pI")
draws_df <- as_draws_df(draws_array)
draws_df<-as.vector(draws_df[1:pars$agrps])
unlist_pI <- vector("numeric", length=pars$agrps)
for(i in 1:pars$agrps){
  unlist_pI[i] <- draws_df[[i]][1]
}

#define and unlist Na from model fitting
Na <- fit$draws("comp_Na")
Na <- as_draws_df(Na)
unlist_Na <- vector("numeric", length=pars$agrps)
for(i in 1:pars$agrps){
  unlist_Na[i] <- Na[[i]][1]
}

#define dprev, ct, etc.
dprev     <- vector("numeric", length=pars$agrps)
seroconv1 <- vector("numeric", length=pars$agrps)
matAb1     <- vector("numeric", length=pars$agrps)
c1       <- vector("numeric", length=pars$agrps)
ct1       <- vector("numeric", length=pars$agrps)

## total deaths, using model-fitted Na
deaths <- sum(pars$d*unlist_Na)

## births distributed among age groups according to fertility
births_age <-  deaths*pars$propfert
births <- sum(births_age)

## calculate conception distribution
for(i in 1:(pars$agrps-1)){
  c1[i] <- births_age[i+1] 
}

# Calculating seroconversions in pregnancy and cases of congenital disease
for(i in 1:(pars$agrps-1)){
  if(i==1){
    dprev[i] <- 0
    seroconv1[i] <- 0
    ct1[i] <- 0
    matAb1[i] <- 0
    
  } else {
    dprev[i] <- unlist_pI[i]-unlist_pI[i-1]                     # change in prevalence (must be positive)
    seroconv1[i] <- dprev[i]*c1[i]                # pregnant women seroconverting in trimester 1
    ct1[i+1] <- seroconv1[i]*mean(pars$mctr)         # likelihood of transmission trimester 1
    matAb1[i+1] <- seroconv1[i]*(1-mean(pars$mctr))  # maternal Ab trimester 1
  }
}

# total number of antibody positive and congenitally diseased births
matAbt <- sum(matAb1)  #17 per year
ctt <- sum(ct1)        #13 per year

#plot estimates
par(mfrow=c(2,2))
plot(dprev)
plot(seroconv1)
plot(ct1)
plot(matAb1)
