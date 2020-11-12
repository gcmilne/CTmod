library("bayesplot")
library("rstanarm")
library("ggplot2")
library("rstan")
library("dplyr")

data <- read.csv("data/netherlands_95.csv")

# data needed for parameter estimation
clean_dat <- matrix(nrow=length(data$age_mid), ncol=3)
clean_dat[,1:3] <- c(data$age_mid, data$k, data$n)
colnames(clean_dat) <- c("age_mid", "k", "n")
data_agrps <- length(clean_dat[,"age_mid"])

# # create larger matrix 
age_groups <- seq(0,80,0.5)  #from 0 -> 80 in 0.5 year intervals
mod_dat <- matrix(nrow=length(age_groups), ncol=3)
colnames(mod_dat) = c("age_mid", "k", "n")
mod_dat[,"age_mid"] <- age_groups
# 
# # merge data into larger matrix for simpler fitting
x <- cbind(t=clean_dat[,"age_mid"], as.data.frame(unname(clean_dat[,2:3])))
y <- cbind(t=mod_dat[,"age_mid"], as.data.frame(unname(mod_dat[,2:3])))
xy <- merge(x, y, by='t', all=TRUE)
xy[is.na(xy)] <- 0
full_data <- data.frame("age_mid"=xy[,1], "k"=xy[,2], "n"=xy[,3])

## v2: creating larger dataset for full range of ages ##
# mod_dat <- matrix(nrow=length(pars$age), ncol=3)
# colnames(mod_dat) = c("age_mid", "k", "n")
# mod_dat[,"age_mid"] <- pars$age
# 
# x <- cbind(t=clean_dat[,"age_mid"], as.data.frame(unname(clean_dat[,2:3])))
# y <- cbind(t=mod_dat[,"age_mid"], as.data.frame(unname(mod_dat[,2:3])))
# xy <- merge(x, y, by='t', all=TRUE)
# xy[is.na(xy)] <- 0
# full_data <- data.frame("age_mid"=xy[,1], "k"=xy[,2], "n"=xy[,3])

##### problem === doing age groups by trimester creates age midpoints that can't match 
##### data from the Netherlands



# data needed for model simulation
agrps <- length(full_data$age_mid)
cases <- full_data$k
n <- full_data$n
  
# other parameters for model
t <- seq(0,10,1)
t0 = 0 
ts <- t[-1]
t <-max(ts)

source("R files/setparms.R")
source("R files/demogdat.R")
N <- sum(pars$Na)
age_prop <- pars$Na/N   ## needs to be same length as nrow(full_data)

# data list for Stan
### EACH OF THESE PARAMETERS NEED TO BE CORRECT LENGTH ##
data_si = list(
  y=y,
  agrps = pars$agrps, 
  data_agrps = data_agrps,
  age_prop=age_prop,
  tot_pop=N, 
  da=pars$da,
  d=pars$d,
  r=pars$r,
  propfert=pars$propfert,
  K=3,  #no. state variables
  t0 = t0,
  ts = ts, 
  t=t,
  n=n, #n
  cases=cases, #k
  inference=1, 
  doprint=0)

# number of MCMC steps
niter <- 100

## Compile the model saved in 'R files/...' ##
model <- stan_model("R files/stan-mod.stan")

## Run MCMC ##
fit_mod <- sampling(object = model, 
                    seed = 123,
                    data = data_si,
                    iter = niter,
                    chains = 2 
                    #control = list(adapt_delta = 0.99)
                    )
system("say -v Karen Surfs up, bro!")

# Warning messages:
#   1: There were 5 divergent transitions after warmup. See
# http://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
# to find out why this is a problem and how to eliminate them. 
# 2: Examine the pairs() plot to diagnose sampling problems
pairs(fit_mod, pars = c("lambda0"), las = 1) # below the diagonal


#################################################
##  check why MCMC rejects initial Stan value ##
#################################################
# init<-vector("numeric",length=agrps*3)
# for(i in 1:agrps){
#   init[i] = (age_prop[i] * (1-seroprev[i]))
#   init[agrps+i] = (age_prop[i] * seroprev[i])
#   init[2*agrps+i] = (age_prop[i] * seroprev[i])
# }
# plot(1:(16*3), init)  ## all values close to 0 — why initial values are rejected?


## Check inference ##
# Specify parameters of interest
pars=c('lambda0')
print(fit_mod, pars = pars)

#summary stats
par_summary <- summary(fit_mod, pars = c("lambda0"), probs = c(0.1, 0.9))$summary
print(par_summary)

# Plot marginal posterior density & confirm Markov chains in agreement with each other
stan_dens(fit_mod, pars = pars, separate_chains = TRUE)

#neat posterior plot
plot_title <- ggtitle("Posterior distributions",
                      "with medians and 95% intervals")
mcmc_areas(fit_mod,
           pars = c("lambda0"),
           prob = 0.95) + plot_title
#trace plot
posterior2 <- extract(fit_mod, inc_warmup = TRUE, permuted = FALSE)
color_scheme_set("mix-blue-pink")
p <- mcmc_trace(posterior2,  pars = "lambda0", n_warmup = 300,
                facet_args = list(nrow = 2, labeller = label_parsed))
p + facet_text(size = 15)
