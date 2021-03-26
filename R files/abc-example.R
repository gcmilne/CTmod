#Example: basic ABC inference of the mean and the standard deviation of a Normal distribution
# Will use a rejection sampling algorithm, and then use a simple MCMC algorithm

#Generate data
number_of_data_points = 100
data =  rnorm(number_of_data_points, mean = 4.3, sd = 2.7)

#set priors on mu and sigma, parameters of the normal distribution
draw_mu <- function () {
  return (runif(1, min=0, max=10))
}
draw_sigma <- function () {
  return (runif(1, min=0, max=10))
}

#simulate a dataset from the parameters
simulate_data <- function (number_of_data_points, mu, sigma) { 
  return(rnorm(number_of_data_points, mean = mu, sd = sigma))
}

# Function to compute the quantiles.#
# We choose to use 3 quantiles.
compute_quantiles <- function(data) {
  return (quantile(data, probs=c(0.1, 0.5, 0.9)))
}
# First method to compare a simulated sample to the observed data
compare_quantiles_with_squared_distance <- function (true, simulated) {
  distance = sqrt(sum(mapply(function(x,y) (x-y)^2, true, simulated)))
  return(distance)
}

# Accept or reject based on the first method to compare a simulated sample to the observed data
accept_or_reject_with_squared_distance <- function (true, simulated, acceptance_threshold) {
  distance = compare_quantiles_with_squared_distance(compute_quantiles(true), compute_quantiles(simulated))
  if((distance < acceptance_threshold) ) return(T) else return(F)
}

## Full rejection sampler
sample_by_rejection <- function (true_data, n_iterations, acceptance_threshold, accept_or_reject_function) {
  number_of_data_points = length(true_data)
  accepted_or_rejected <- vector(length = n_iterations)
  sampled_mus <- vector(length = n_iterations, mode = "numeric")
  sampled_sigmas <- vector (length = n_iterations, mode = "numeric")
  for (i in 1:n_iterations){
    mu <- draw_mu()
    sigma <- draw_sigma()
    parameters = list("mu"=mu, "sigma"=sigma)
    simulated_data <- simulate_data(number_of_data_points, mu, sigma)
    accepted_or_rejected[i] = accept_or_reject_function(true_data, simulated_data, acceptance_threshold)
    sampled_mus[i] = mu
    sampled_sigmas[i] = sigma
  }
  return(data.frame(cbind("accepted_or_rejected" = accepted_or_rejected, "sampled_mus" = sampled_mus, "sampled_sigmas" = sampled_sigmas)))
}

#Now we perform inference using our distance functions, and analyze the results
#First, the squared distance:
sampled_parameter_values_squared_distances = sample_by_rejection(data, 20000, 0.5, accept_or_reject_with_squared_distance)

#How many samples have been accepted?
sum(sampled_parameter_values_squared_distances$accepted_or_rejected)

#Using Coda to summarize the samples and do some plots
# Useful library used here for plotting mostly.
# install.packages("coda")
library(coda)
rej_samples_squared_distances_as_mcmc = mcmc(sampled_parameter_values_squared_distances[which(sampled_parameter_values_squared_distances$accepted_or_rejected==1),c(2,3)])

summary(rej_samples_squared_distances_as_mcmc)

plot(rej_samples_squared_distances_as_mcmc)


##############################################
########### attempting with my model #########
##############################################

# Load scripts #
source("R files/demogdat.R")
source("R files/setparms.R")
source("R files/model.R")

# Read in data #
data <- read.csv("data/netherlands_95.csv")
number_of_data_points = length(data$n)

# set prior on lambda0 #
set.seed(1001)
lambda0_sample <- rnorm(3, mean=0.05, sd=0.01)
# hist(lambda0_sample)

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
# par(mfrow=c(2,2))
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
  list_mod_matched[[i]]$prev    <- store_sim[[i]][,"pI"][matched_indices]
  list_mod_matched[[i]]$lambda0 <- lambda0_sample[i]
}

# plot modelled vs. observed
par(mfrow=c(2,2))
for(i in 1:nsim){
  plot(list_mod_matched[[i]]$age_mid, list_mod_matched[[i]]$prev, ylim=c(0,1), 
       type="l", xlab="age (years)", ylab="seroprevalence", main = signif(list_mod_matched[[i]]$lambda0[1], 3))
  points(data$age_mid, data$prev)
}

