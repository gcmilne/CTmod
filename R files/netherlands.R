# Set demography to match years of data collection in 
# Hofhuis et al., 2011. Epidemiology & Infection 139(4)

## Netherlands '95/96 data
neth <- read.csv("data/netherlands_demog.csv")
fract <- vector("numeric", length=nrow(neth))
for(i in 1:nrow(neth)){
  fract[i] <- neth$N_95[i]/sum(neth$N_95)
}

data <- read.csv("data/netherlands_temporal.csv")

n <- round(fract*data$n_95[1], 0)  # presume study sampling representative of country's age-demography pattern 
k <- round(n*data$prev_95, 0)
prev <- k/n
age <- data$age_mid
neth_95 <- data.frame(age, k, n, prev)

## Netherlands '06/07 data
fract <- vector("numeric", length=nrow(neth))
for(i in 1:nrow(neth)){
  fract[i] <- neth$N_05[i]/sum(neth$N_05)
}

n <- round(fract*data$n_05[1], 0)  # presume study sampling representative of country's age-demography pattern 
k <- round(n*data$prev_05, 0)
prev <- k/n
neth_05 <- data.frame(age, k, n, prev)

# compare 95 to 05 data
plot(neth_95$age, neth_95$prev, type='l', ylim=c(0,1))
lines(neth_05$age, neth_05$prev, lty=2)