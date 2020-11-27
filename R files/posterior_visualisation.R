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
