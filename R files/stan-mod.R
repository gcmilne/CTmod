library("bayesplot")
library("rstanarm")
library("ggplot2")
library("rstan")
library("dplyr")

#read in data
data <- read.csv("data/netherlands_95.csv")

# data needed for parameter estimation
clean_dat <- data.frame("age_mid"=data$age_mid, "k"=data$k, "n"=data$n)
data_agrps <- length(clean_dat[,"age_mid"])

######################################################
## modify data age groups to match model age groups ##
######################################################
#create new dataset
matched_dat <- clean_dat
# x[,2] increases by one each time data age midpoint is closest match to modelled age midpoint
x <- cbind(pars$age, findInterval(pars$age, matched_dat$age_mid))
#head(x, n=12)
#returns FALSE if there's change between element i and element i+1
y <- diff(x[,2]) <= 0   #so save i+1 element of x[,1] to matched_ages[i]
#head(y, n=12)

# each time x[,2] increases by 1, save value of x[,1][i+1] to matched_dat$age_mid[i]
matched_ages <- vector("numeric", pars$agrps)
for(i in 1:length(y)){
  if(y[i]==T){
    matched_ages[i] <- NA
    
  }else if(y[i]==F){
    matched_ages[i] <- x[,1][i+1]
  }
}
#removes last element (which is 0 because of indexing)
matched_ages <- head(matched_ages, -1)
#remove NAs
matched_ages <-matched_ages[!is.na(matched_ages)]
#save age mid points to dataset
matched_dat$age_mid <- matched_ages

###################################################
## merge data into larger df for simpler fitting ##
###################################################
mod_dat <- data.frame("age_mid"=pars$age,"k"=0, "n"=0)
merged_dat <- merge(matched_dat, mod_dat, by='age_mid', all=TRUE)
merged_dat[is.na(merged_dat)] <- 0
full_data <- data.frame("age_mid"=merged_dat[,1], "k"=merged_dat[,2], "n"=merged_dat[,3])

############################################
# read in data needed for model simulation #
############################################
agrps <- length(full_data$age_mid)
cases <- full_data$k
n <- full_data$n

######################################
# read in other parameters for model #
######################################
t <- seq(0,10,1)
t0 = 0 
ts <- t[-1]
t <-max(ts)

source("R files/setparms.R")
source("R files/demogdat.R")
N <- sum(pars$Na)
age_prop <- pars$Na/N   ## needs to be same length as nrow(full_data)

###################### 
# data list for Stan #
######################
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


#####################
## check inference ##
#####################
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
