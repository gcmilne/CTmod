library("bayesplot")
library("rstanarm")
library("ggplot2")
library("rstan")
library("dplyr")

source("R files/demogdat.R")
source("R files/setparms.R")

data <- read.csv("data/netherlands_95.csv")

# data needed for parameter estimation
clean_dat <- matrix(nrow=length(data$age_mid), ncol=3)
clean_dat[,1:3] <- c(data$age_mid, data$k, data$n)
colnames(clean_dat) <- c("age_mid", "k", "n")

# create larger matrix 
age_groups <- seq(0,80,0.5)  #from 0 -> 80 in 0.5 year intervals
mod_dat <- matrix(nrow=length(age_groups), ncol=3)
colnames(mod_dat) = c("age_mid", "k", "n")
mod_dat[,"age_mid"] <- age_groups

# merge data into larger matrix for simpler fitting
x <- cbind(t=clean_dat[,"age_mid"], as.data.frame(unname(clean_dat[,2:3])))
y <- cbind(t=mod_dat[,"age_mid"], as.data.frame(unname(mod_dat[,2:3])))
xy <- merge(x, y, by='t', all=TRUE)
xy[is.na(xy)] <- 0
full_data <- data.frame("age_mid"=xy[,1], "k"=xy[,2], "n"=xy[,3])

# data needed for model simulation
agrps <- length(full_data$age_mid)
cases <- full_data$k
n <- full_data$n
#seroprev <- cases/n
#seroprev[is.na(seroprev)] <- 0
  
# other parameters for model
t <- seq(0,200,1)
t0 = 0 
ts <- t[-1]
t <-max(ts)

tot_pop <- sum(pars$Na)
age_prop <- pars$Na/tot_pop   ## needs to be same length as nrow(full_data)
  
# data list for Stan
data_si = list(
  agrps = agrps, 
  age_prop=age_prop,
  tot_pop=tot_pop, 
  K=3,  #no. state variables
  t0 = t0,
  ts = ts, 
  t=t,
  p_lambda0=c(0.5, 100),  #priors for lambda0
  p_seroprev=c(2, 5), #priors for seroprevalence
  n=n,          # n
  cases=cases,  # k
  inference=1, 
  doprint=0)

# number of MCMC steps
niter <- 2000

## Compile the model saved in 'R files/...' ##
model <- stan_model("R files/stan-basic.stan")

## Run MCMC ##
fit_si_negbin <- sampling(object = model, 
                          seed = 123,
                          data = data_si,
                          iter = niter,
                          chains = 1)

#################################################
##  check why MCMC rejects initial Stan value ##
#################################################
init<-vector("numeric",length=agrps*3)
for(i in 1:agrps){
  init[i] = (age_prop[i] * (1-seroprev[i]))
  init[agrps+i] = (age_prop[i] * seroprev[i])
  init[2*agrps+i] = (age_prop[i] * seroprev[i])
}
plot(1:(16*3), init)  ## all values close to 0 — why initial values are rejected?


## Check inference ##
# Specify parameters of interest
pars=c('lambda0')
print(fit_si_negbin, pars = pars)

#summary stats
par_summary <- summary(fit_si_negbin, pars = c("lambda0"), probs = c(0.1, 0.9))$summary
print(par_summary)

# Plot marginal posterior density & confirm Markov chains in agreement with each other
stan_dens(fit_si_negbin, pars = pars, separate_chains = TRUE)

#neat posterior plot
plot_title <- ggtitle("Posterior distributions",
                      "with medians and 95% intervals")
mcmc_areas(fit_si_negbin,
           pars = c("lambda0"),
           prob = 0.95) + plot_title
#trace plot
posterior2 <- extract(fit_si_negbin, inc_warmup = TRUE, permuted = FALSE)
color_scheme_set("mix-blue-pink")
p <- mcmc_trace(posterior2,  pars = "lambda0", n_warmup = 300,
                facet_args = list(nrow = 2, labeller = label_parsed))
p + facet_text(size = 15)
