#######################################################
#     MCMC : 3 SC VILLAGES - LONG TERM - V2           #
#######################################################

###Packages###
library(tidyverse) #graph
library(rstan)     #implements MCMC
library(shinystan) #diagnostics
library(readxl)    #excel
library(loo)       # WAIC

###Data###
#Read data
Slash_trial_LT2 <- read_excel("S&CT_DS1.xlsx",sheet = "Long term 2_UG") 
head(Slash_trial_LT2)

#Extract data
data_patille_DBR <- Slash_trial_LT2$`PATILLE (I)` #extraction of DBR
data_pwomunu_DBR <- Slash_trial_LT2$`PWOMUNU  (I)`
data_elugub_DBR <- Slash_trial_LT2$`ELUGU B (I)`

data_patille_time <- Slash_trial_LT2$Day            
data_patille_time_eq <- data_patille_time + 1500   
data_patille_time2 <- data_patille_time[!is.na(data_patille_DBR)]
data_patille_time_eq2 <- data_patille_time_eq[!is.na(data_patille_DBR)]

data_patille_DBR2 <- data_patille_DBR[!is.na(data_patille_DBR)]
data_pwomunu_DBR2 <- data_pwomunu_DBR[!is.na(data_patille_DBR)]
data_elugub_DBR2 <-data_elugub_DBR[!is.na(data_patille_DBR)]

#Initial conditions
E_0 = 20000
L1_0 = 5000
L2_0 = 5000
L3_0 = 5000
L4_0 = 5000
L5_0 = 5000
L6_0 = 5000
L7_0 = 5000
P_0 = 50000
N_0 = 10000
Pa_0 = 7500
initialpop = c(E = E_0, L1 = L1_0, L2 = L2_0, L3 = L3_0, L4 = L4_0, L5 = L5_0, L6 = L6_0, L7 = L7_0, P = P_0, N = N_0, Pa = Pa_0)
# data for Stan
data_stan <- list(n_days = length(data_patille_DBR2), initialpop = initialpop, t0 = 0,
                  ts = data_patille_time_eq2, DBRs_patille = data_patille_DBR2, DBRs_pwomunu = data_pwomunu_DBR2, DBRs_eluguB = data_elugub_DBR2,
                  muEo = 0.05, EpsN = 431.73, EpsP = 141.89,  #from Routledge et al.
                  Em = 0.17, muL = 0.24, muP = 0.1, #from Routledge et al.
                  HdivHBI = 541.69, temp = 26.34, g = 3.48, #from Palulu
                  #rainfallmortality1_bajere = 0.09, rainfalldecay1_bajere = 0.00, rainfallmortality2_bajere = 0.16, rainfalldecay2_bajere = 0.01, #from control villages (Sympa)
                  #rainfallmortality1_okidi = 0.09, rainfalldecay1_okidi = 0.00,rainfallmortality2_okidi = 0.15, rainfalldecay2_okidi = 0.00, #from control villages
                  #rainfallmortality1_elugua1 = 0.09, rainfalldecay1_elugua1 = 0.00, rainfallmortality2_elugua1 = 0.17, rainfalldecay2_elugua1 = 0.01) #from control villages
                  #DBRstar_bajere = 78.99, DBRstar_okidi = 101.24, DBRstar_elugua1 = 108.89) #from control villages
                 rainfallmortality1_bajere = 0.10, rainfallmortality2_bajere = 0.17, rainfalldecay2_bajere = 0.01, #from control villages
                 rainfallmortality1_okidi = 0.11, rainfallmortality2_okidi = 0.17, rainfalldecay2_okidi = 0.01, #from control villages
                 rainfallmortality1_elugua1 = 0.11, rainfallmortality2_elugua1 = 0.19, rainfalldecay2_elugua1 = 0.01)
                  # DBRstar_bajere = 69.36, DBRstar_okidi = 86.78, DBRstar_elugua1 = 93.04) #from control villages
### Execution of the model ###
# run in parallel on multiple processors
options(mc.cores=parallel::detectCores())
# number of MCMC steps
niter <- 200
#compilation
rstan_options(auto_write=TRUE)

###################################################### Want to run this for my model and get useful output that can fit to data
model_patille <- stan_model("LT_SCVillagesV2.stan") ##
######################################################

#run MCMC (the NUTS algorithm is used, a improved version of HMC)
fit_patille <- sampling(model_patille,data = data_stan,
                        iter = niter,chains = 1, #diagnostic_file = "Diag_patille",
                        cores = 4, algorithm = "NUTS")
print(fit_patille,probs = c(0.25,0.5,0.75),digits = 5)

###Diagnostics###

#png(filename = "Essai_patille")
#Dataframe 1 : Patille
DF_DBR_pred_patille <- cbind(as.data.frame(summary(
  fit_patille, pars = "pred_DBR_patille", probs = c(0.025, 0.5, 0.975))$summary),
  data_patille_time2, data_patille_DBR2)
colnames(DF_DBR_pred_patille) <- make.names(colnames(DF_DBR_pred_patille)) # to remove % in the col names 
#Posterior predictive check
ggplot(DF_DBR_pred_patille, mapping = aes(x = data_patille_time2)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), fill = "orange", alpha = 0.6) +
  geom_line(mapping = aes(x = data_patille_time2, y = X50.)) + 
  geom_point(mapping = aes(y = data_patille_DBR2)) +
  labs(x = "Day", y = "Daily Biting Rate", title="Posterior predictive check (Patille)")
#dev.off()

#png(filename = "Essai_pwomunu")
#Dataframe 2 : Pwomunu
DF_DBR_pred_pwomunu <- cbind(as.data.frame(summary(
  fit_patille, pars = "pred_DBR_pwomunu", probs = c(0.025, 0.5, 0.975))$summary),
  data_patille_time2, data_pwomunu_DBR2)
colnames(DF_DBR_pred_pwomunu) <- make.names(colnames(DF_DBR_pred_pwomunu)) # to remove % in the col names 
#Posterior predictive check
ggplot(DF_DBR_pred_pwomunu, mapping = aes(x = data_patille_time2)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), fill = "orange", alpha = 0.6) +
  geom_line(mapping = aes(x = data_patille_time2, y = X50.)) + 
  geom_point(mapping = aes(y = data_pwomunu_DBR2)) +
  labs(x = "Day", y = "Daily Biting Rate", title="Posterior predictive check (Pwomunu)")
#dev.off()

#png(filename = "Essai_elugub")
#Dataframe 3 : Elugu B
DF_DBR_pred_elugub <- cbind(as.data.frame(summary(
  fit_patille, pars = "pred_DBR_eluguB", probs = c(0.025, 0.5, 0.975))$summary),
  data_patille_time2, data_elugub_DBR2)
colnames(DF_DBR_pred_elugub) <- make.names(colnames(DF_DBR_pred_elugub)) # to remove % in the col names 
#Posterior predictive check
ggplot(DF_DBR_pred_elugub, mapping = aes(x = data_patille_time2)) +
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), fill = "orange", alpha = 0.6) +
  geom_line(mapping = aes(x = data_patille_time2, y = X50.)) + 
  geom_point(mapping = aes(y = data_elugub_DBR2)) +
  labs(x = "Day", y = "Daily Biting Rate", title="Posterior predictive check (Elugu B)")
#dev.off()
vecx <- 1511:1852
plot(x= vecx, y = -(log(1-0.48)/2)*exp(-0.09*(vecx-1511)))


#Diagnostics
#sFit <- as.shinystan(fit_patille)
#launch_shinystan(sFit)


###Warning messages###
