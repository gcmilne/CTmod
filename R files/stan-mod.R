#library("bayesplot")
#library("rstanarm")
library("ggplot2")
library("rstan")
library("dplyr")

#install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
library(cmdstanr)
library(posterior)
library(bayesplot)

#read in data
data <- read.csv("data/netherlands_95.csv")
source("R files/setparms.R")
source("R files/demogdat.R")

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
t0 = 0 
ts <- seq(1,10, 1)
t <-max(ts)
N <- sum(pars$Na)
age_prop <- pars$Na/N   ## needs to be same length as nrow(full_data)

#index of rows of expanded df in which data exist
data_rows <- c(which(full_data$k!=0), which(full_data$k!=0)+agrps, which(full_data$k!=0)+2*agrps)

d <- pars$d

###################### 
# data list for Stan #
######################
data_si = list(
  agrps = pars$agrps, 
  data_agrps = data_agrps,
  data_rows=data_rows,
  age_prop=age_prop,
  tot_pop=N, 
  # S_r = 5, #no real arrays
  # starts_r = c(1, 2, 3, 6, 6+pars$agrps, 6+(2*pars$agrps)),  #start indices of real arrays
  # r_array=c(pars$r, pars$da, pars$mctr, pars$d, pars$propfert), 
  age=pars$age,
  da=pars$da,
  d=pars$d,
  r=pars$r,
  mctr=pars$mctr,
  propfert=pars$propfert,
  K=3,  #no. state variables
  t0 = t0,
  ts = ts, 
  t=t,
  n=n, #n
  cases=cases, #k
  rel_tol = 1.0E-10, 
  abs_tol = 1.0E-10,
  max_num_steps = 1.0E3,
  inference=0, 
  doprint=1)

###################
# CmdStan running #
###################
file <- "R files/stan-mod-simple.stan"
mod <- cmdstan_model(file)

fit <- mod$sample(
  data = data_si,
  seed = 123,
  chains = 2,
  parallel_chains = 1, 
  iter_warmup = 10,
  iter_sampling = 40,
  refresh = 1
)
# system("say -v Karen Surfs up, bro!")

fit$summary()
draws_array <- fit$draws()
str(draws_array)
draws_df <- as_draws_df(draws_array) # as_draws_matrix() for matrix
print(draws_df)
mcmc_hist(fit$draws("theta"))

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


# For info on posterior predictive checks (including plotting modelled estimates vs. data), see:
# https://arxiv.org/pdf/2006.02985.pdf
