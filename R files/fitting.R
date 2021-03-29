########### Model fitting

# Load scripts #
source("R files/demogdat.R")
source("R files/setparms.R")
source("R files/model.R")

# Read in data #
data <- read.csv("data/netherlands_95.csv")
number_of_data_points = length(data$n)

# set prior on lambda0 #
set.seed(1001)
lambda0_sample <- runif(6, min=0, max=0.10)
# hist(lambda0_sample)

# set prior on lambda1 #
lambda1_sample <- runif(10000, min=1e-6, max=1e-5)
# hist(lambda1_sample)

# set prior on gradient #
lambda1_sample <- runif(10000, min=1e-6, max=1e-5)
# hist(lambda1_sample)


pars$lambda0 <- 0
pars$lambda1 <- 0.002
pars$gradient <- 0.008
# pars$shape <- 1e-5

# new model foi form
# foi <- (pars$lambda0 + pars$lambda1*(pars$age^2) * (pars$age * exp(-pars$gradient*pars$age)))*pars$shape

# tgerp model foi form
foi <- pars$lambda0 + pars$lambda1 * (pars$age * exp(-pars$gradient*pars$age))
plot(pars$age, foi)



### simulate a dataset from the parameters
nsim <- length(lambda0_sample) # no. of simulations
sol_age <- data.frame(matrix(NA, nrow=length(pars$age)-1, ncol = 16)) #dimensions of output from age profile function
names(sol_age) <- names(df) #set names
store_sim <- rep(list(sol_age, nsim)) # create list to store output from all simulations

### loop through samples from prior and store model output
for(i in 1:nsim){
  pars$lambda0 <- lambda0_sample[i] #replace lambda0 with sample from the prior
  sol <- ode(times = time, y = y,  parms = pars, func = age_si) #run ODE model
  store_sim[[i]] <- getit(max(time)) #store age profile after burnin period
}

### plot results for each of the simulations
# par(mfrow=c(2,3))
# for(i in 1:nsim){
#   plot(store_sim[[i]]$a, store_sim[[i]]$pI, type="l", xlab = "age (years)", ylab = "prevalence",
#        main = signif(lambda0_sample[i], 3))
# }

### select age groups from model output that match data age groups
clean_dat <- data.frame("age_mid"=data$age_mid, "k"=data$k, "n"=data$n, "prev"=data$prev)

# first create a list to store the model output from the correct age categories #
mod_matched <- data.frame(matrix(NA, nrow=length(data$age_mid), ncol = 3))
names(mod_matched) <- c("age_mid", "prev", "lambda0")
list_mod_matched <- rep(list(mod_matched), nsim)

# then find the closest age categories #
#create new dataset
matched_dat <- clean_dat
# x[,2] increases by one each time data age midpoint is closest match to modelled age midpoint
x <- cbind(pars$age, findInterval(pars$age, matched_dat$age_mid))
#returns FALSE if there's change between element i and element i+1
y1 <- diff(x[,2]) <= 0

# each time x[,2] increases by 1, save value of x[,1][i+1] to matched_dat$age_mid[i]
matched_ages <- vector("numeric", pars$agrps)
for(i in 1:length(y1)){
  if(y1[i]==T){
    matched_ages[i] <- NA
  }else if(y1[i]==F){
    matched_ages[i] <- x[,1][i+1]
  }
}

#removes last element (which is 0 because of indexing)
matched_ages <- head(matched_ages, -1)

# find the indices which match the age group most closely
matched_indices <- which(!is.na(matched_ages))

# store the relevant model output in a new list
for(i in 1:nsim){
  list_mod_matched[[i]]$age_mid <- store_sim[[i]][,"a"][matched_indices]
  list_mod_matched[[i]]$prev    <- store_sim[[i]][,"obs_pI"][matched_indices]
  list_mod_matched[[i]]$lambda0 <- lambda0_sample[i]
}

# plot modelled vs. observed
# par(mfrow=c(3,2))
# for(i in 1:nsim){
#   plot(list_mod_matched[[i]]$age_mid, list_mod_matched[[i]]$prev, ylim=c(0,1), 
#        type="l", xlab="age (years)", ylab="seroprevalence", main = signif(list_mod_matched[[i]]$lambda0[1], 3))
#   points(data$age_mid, data$prev)
# }

### Next: likelihood function

loglik <- function(kIgG, nIgG, prevIgG){ 
  dbinom(k, n, prev, log=T) 
}

## Likelihood wrapper function
loglik.wrap <- function(pars, data){ 
  
  pars$lambda0  <- par[1]
  pars$lambda1  <- par[2]
  pars$gradient <- par[3]

  sol          <- ode(y = state, times = time, parms = pars,  func = age_si)  #save model solution
  store_sim    <- getit(max(time))  #store age profile after burnin period
  matched_prev <- store_sim[,"obs_pI"][matched_indices]  #select observed prevalence from relevant age categories
  
  logliks <- loglik(k = data$k, n = data$n, prev = matched_prev)
  
  return(sum(-logliks))
}

## Check it works
par <- c(log(0.0002),log(0.03), log(0.05), log(0.0002))
#par <- c(0.02, 0.05, 0.1)

test<-loglik.wrap(par=par, data=v_data)
test


# Accept or reject based on the first method to compare a simulated sample to the observed data
accept_or_reject_with_squared_distance <- function (true, simulated, acceptance_threshold) {
  distance = compare_quantiles_with_squared_distance(compute_quantiles(true), compute_quantiles(simulated))
  if((distance < acceptance_threshold) ) return(T) else return(F)
}

