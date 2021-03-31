########### Model fitting

# Load scripts #
source("R files/demogdat.R")
source("R files/setparms.R")
source("R files/model.R")

# Load packages
library(deSolve)
require(lhs) #for latin hypercube sampling

# Read in data #
data <- read.csv("data/netherlands_95.csv")
number_of_data_points = length(data$n)

### select age groups from model output that match data age groups
clean_dat <- data.frame("age_mid"=data$age_mid, "k"=data$k, "n"=data$n, "prev"=data$prev)

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

## Function to get age profiles at given time
getit <- function(time) {
  row <- which(abs(sol[,"time"]-time)==min(abs(sol[,"time"]-time)))
  df <- sol[row,-1]
  S  <- df[1:pars$agrps]  
  I  <- df[(pars$agrps+1):(2*pars$agrps)]  
  Im <- df[(2*pars$agrps+1):(3*pars$agrps)]
  pI <- df[(3*pars$agrps+1):(4*pars$agrps)]
  obs_pI <- df[(4*pars$agrps+1):(5*pars$agrps)]
  dprev <- df[(5*pars$agrps+1):(6*pars$agrps)]
  seroconv1 <- df[(6*pars$agrps+1):(7*pars$agrps)] 
  seroconv2 <- df[(7*pars$agrps+1):(8*pars$agrps)] 
  seroconv3 <- df[(8*pars$agrps+1):(9*pars$agrps)] 
  matAb1 <- df[(9*pars$agrps+1):(10*pars$agrps)]
  matAb2 <- df[(10*pars$agrps+1):(11*pars$agrps)]
  matAb3 <- df[(11*pars$agrps+1):(12*pars$agrps)]
  ct1 <- df[(12*pars$agrps+1):(13*pars$agrps)]
  ct2 <- df[(13*pars$agrps+1):(14*pars$agrps)]
  ct3 <- df[(14*pars$agrps+1):(15*pars$agrps)]
  Na <- df[(15*pars$agrps+1):(16*pars$agrps)]
  age <- pars$age  
  out <- data.frame(a=age, I=I,Im=Im, pI=pI, obs_pI=obs_pI, dprev=dprev, 
                    seroconv1=seroconv1, seroconv2=seroconv2, seroconv3=seroconv3, 
                    matAb1=matAb1, matAb2=matAb2, matAb3=matAb3, 
                    ct1=ct1, ct2=ct2, ct3=ct3, Na=Na)
  ## remove last age category for aesthetics
  out <- out[-nrow(out),]
}

### Likelihood function
loglik <- function(k, n, prev){ 
  dbinom(k, n, prev, log=T) 
}

## Latin hypercube sampling
set.seed(1001)
nsim <- 3

par_arr <- randomLHS(nsim, 3) #create parameter array
par_arr[,1] <- log(qunif(par_arr[,1], min=0, max=0.2))        #log lambda0
par_arr[,2] <- log(rbeta(par_arr[,2], shape1 = 2, shape2=80)) #log lambda1
par_arr[,3] <- log(runif(par_arr[,3], min=0, max=0.2))        #log gradient

#par(mfrow=c(2,2))
#dummy <- apply(exp(par_arr), 2, hist, main = "") #plot prior distributions

lik_arr <- vector(mode="numeric", length=nsim) #create likelihood array
matched_prev <- rep(list(matrix(nrow=length(data$age_mid), ncol=2)),nsim) #create list to store model output matched to correct age bins
names(matched_prev) <- c("prev", "ct")

# Run model and store likelihoods
start.time <- Sys.time()
for(i in 1:nrow(par_arr)){
  pars$log.lambda0  <- par_arr[i,1]
  pars$log.lambda1  <- par_arr[i,2]
  pars$log.gradient <- par_arr[i,3]
  sol          <- ode(y = y, times = time, parms = pars,  func = age_si)  #save model solution
  store_sim    <- getit(max(time))  #store age profile after burnin period
  matched_prev[[i]]["prev"] <- store_sim[,"obs_pI"][matched_indices]  #select observed prevalence from relevant age categories
  matched_prev[[i]]["ct"] <- (store_sim[,"ct1"][matched_indices]+store_sim[,"ct2"][matched_indices]+store_sim[,"ct3"][matched_indices]) #select ct cases from relevant age categories
  logliks <- loglik(k = data$k, n = data$n, prev = matched_prev[[i]]) #run likelihood function
  lik_arr[i] <- sum(-logliks)
}

end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

# plot
bfit <- matched_prev[[which.min(lik_arr)]]["prev"] #prevalence with min likelihood
plot(data$age_mid, bfit, type='l')
points(data$age_mid, data$prev)


# time is in days, so need to sum ct cases over 365 days to get total for the year







## Plot priors ##
nsim <- 1000
# set prior on lambda0 #
set.seed(1001)
plot(density(par_arr[,1]), main = "lambda0 prior")
polygon(density(par_arr[,1]), col = "lightblue")

# set prior on lambda1 #
plot(density(par_arr[,2]), main = "lambda1 prior")
polygon(density(par_arr[,2]), col = "lightblue")

# set prior on gradient #
plot(density(par_arr[,3]), main = "gradient prior")
polygon(density(par_arr[,3]), col = "lightblue")

