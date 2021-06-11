# Clear environment
rm(list = ls())

# load required packages
library(MortalityLaws)
library(deSolve)
library(MASS)
library(ggplot2)
library(binom)
library(usethis)
library(wpp2019)

# library.dynam.unload("deSolve", libpath=paste(.libPaths()[1], "//deSolve", sep=""))
# library.dynam("deSolve", package="deSolve", lib.loc=.libPaths()[1])


# Add github environment
#edit_r_environ()  # GITHUB_PAT = 'dac92469d116f55dac8162b930c76aedde600b32'
#use_github(protocol = 'https', auth_token = Sys.getenv('GITHUB_PAT'))

source("R files/demogdat.R")
source("R files/setparms.R")
source("R files/model.R")

## for loop to calculate changing no. of CT cases with endemicity
# Set parms using best fit parameters for Netherlands '95/'96 data
# endemicity <- seq(0.0001, 0.01, by=0.0005)
# ct_cases<- vector("numeric", length= length(endemicity))
# for(i in 1:length(endemicity)){
#   pars$shape <- endemicity[i]
#   sol= ode(times = time, y = y,  parms = pars, func = age_si)
#   tp <- getit(249)
#   ct_cases[i] <- sum(tp[,"ct1"] + tp[,"ct2"] + tp[,"ct3"])
# }
# foi <- rep(0, length(endemicity))
# df <- setNames(data.frame(endemicity, ct_cases, foi), c("shapevalue", "ct_cases", "foi"))
# write.csv(df,"data/ct-changing-foi.csv", row.names = FALSE)

# plot changing CT cases with endemicity level
par(mfrow=c(1,1))
plot(endemicity, ct_cases, xlab = "shape value", ylab = "congenital cases", type="l")
points(5.647195e-04, 645.737, col = "red")  # best-fit values for Netherlands '95
points((5.647195e-04)*0.55, 444.8015, col = "red")  # best-fit values for Netherlands '06/07
text(x=(5.647195e-04)+0.0005, y=645.737, labels = "'95")
text(x=(5.647195e-04*.55)+0.0007, y=444.8015, labels = "'06/'07")

plot((pars$lambda0 + pars$lambda1*(pars$age^2) * (pars$age * exp(-pars$gradient*pars$age)))*0.0011,
     ylab="force of infection", xlab="age", ylim=c(0,0.07), type="l")

## solve model
start.time <- Sys.time()
#= ode(times = time, y = y,  parms = pars, func = rcpp_si)
sol= ode(times = time, y = y,  parms = pars, func = age_si)
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken
system("say -v Karen Surfs up, bro!")

## Time dynamics
# pdf(file="burnin-period.pdf")
par(mfrow=c(2,2))
# par(mar=rep(4, 4))

# plot(sol[,"Nt"]~sol[,"time"], type="l", ylab="pop size", xlab="time")
plot(sol[,"St"]~sol[,"time"], type="l", ylab="suceptible", xlab="time")
plot(sol[,"It"]~sol[,"time"], type="l", ylab="infected", xlab="time")
#plot(sol[,"Imt"]~sol[,"time"], type="l", ylab="maternal Ab +", xlab="time")
plot(sol[,"ctt"]~sol[,"time"], type="l", ylab="ct", xlab="time")
plot(sol[,"pIt"]~sol[,"time"], type="l", ylab="prevalence", xlab="time")

# dev.off()

## Age profiles at given time
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
  # out <- out[-nrow(out),]
}

par(mfrow=c(1,1))
df <- getit(150)
plot(df[,"Na"]~df[,"a"], type="l", ylab="pop size", xlab="age", ylim = c(20000, 70000))
points(pars$Na~pars$age)

plot(df[,"pI"]~df[,"a"], type="l", ylab="seroprevalence", xlab="age", ylim = c(0, 1))


df <- getit(1)
plot(df[,"Na"]~df[,"a"], type="l", ylab="pop size", xlab="age")
points(pars$Na~pars$age)
df <- getit(45)
plot(df[,"Na"]~df[,"a"], type="l", ylab="pop size", xlab="age", ylim=c(0,250000))
points(pars$Na~pars$age)
df <- getit(50)
plot(df[,"Na"]~df[,"a"], type="l", ylab="pop size", xlab="age", ylim=c(0,250000))
points(pars$Na~pars$age)
df <- getit(60)
plot(df[,"Na"]~df[,"a"], type="l", ylab="pop size", xlab="age", ylim=c(0,250000))
points(pars$Na~pars$age)
df <- getit(249)
plot(df[,"Na"]~df[,"a"], type="l", ylab="pop size", xlab="age", ylim=c(0,80000))
points(pars$Na~pars$age)

data <- data.frame(df[,"a"], df[,"Na"])

df <- getit(100)
par(mfrow=c(1,1))
#plot(df[,"pI"]~df[,"a"], type="l", ylab="prevalence", xlab="age", main="shape= .015")
#plot(df[,"pI"]~df[,"a"], type="l", ylab="prevalence", xlab="age", main="shape= .001")
par(mfrow=c(2,1))
plot(df[,"pI"]~df[,"a"], type="l", ylab="prevalence", xlab="age")
points(neth95$age, neth95$prev)
df <- getit(500)
plot(df[,"pI"]~df[,"a"], type="l", ylab="prevalence", xlab="age")
points(neth05$age, neth05$prev)

sum(df[,"ct1"])
sum(df[,"ct2"])
sum(df[,"ct3"])

df <- getit(249)
par(mfrow=c(3,2))
#plot(df[,"Na"]~df[,"a"], type="l", ylab="pop size", xlab="age")
#points(pars$Na~pars$age)
#plot((df[,"S"])~df[,"a"], type="l", ylab="susceptible", xlab="age")
plot(df[,"I"]~df[,"a"], type="l", ylab="infected", xlab="age")
plot(df[,"Im"]~df[,"a"], type="l", ylab="maternal Ab +", xlab="age", xlim=c(0,1))
plot(df[,"pI"]~df[,"a"], type="l", ylab="prevalence", xlab="age")
lines(df[,"obs_pI"]~df[,"a"], type="l", ylab="observed prevalence", xlab="age", lty=2, ylim=c(0,1))
#legend("topleft", c("Prevalence", "Observed prevalence"), lty=1:2)
plot(df[-nrow(df),"Na"]~df[-nrow(df),"a"], type="l", ylab="pop size", xlab="age", ylim=c(0,250000))
points(pars$Na~pars$age)
plot(df[,"seroconv3"]~df[,"a"], type="l", ylab="no. seroconversions", xlab="age")
lines(df[,"matAb3"]~df[,"a"], type="l", ylab="no. cases ct", xlab="age", lty=2)

# Total maternal Ab 
par(mfrow=c(1,1))
lines((df[,"matAb3"]+df[,"matAb2"]+df[,"matAb1"])~df[,"a"], type="l", ylab="no. cases", xlab="age", lty=1)
# seroconversions
plot((df[,"seroconv3"]+df[,"seroconv2"]+df[,"seroconv1"])~df[,"a"], type="l", ylab="no. cases", xlab="age", lty=2)
# Total CT cases
lines((df[,"ct3"]+df[,"ct2"]+df[,"ct1"])~df[,"a"], type="l", ylab="no. cases ct", xlab="age", lty=3)
legend("topright", c("Maternal IgG",  "seroconversions", "CT cases"), lty=1:3)

plot(((df[,"ct3"]+df[,"ct2"]+df[,"ct1"])+(df[,"matAb3"]+df[,"matAb2"]+df[,"matAb1"]))~df[,"a"], type="l", ylab="no. cases ct", xlab="age", lty=1)
lines((df[,"seroconv3"]+df[,"seroconv2"]+df[,"seroconv1"])~df[,"a"], type="l", ylab="no. cases", xlab="age", lty=2)


# CT cases by trimester of infection
par(mfrow=c(1,1))
plot(df[,"ct3"]~df[,"a"], type="l", ylab="no. cases ct", xlab="age", lty=1)
lines(df[,"ct2"]~df[,"a"], type="l", ylab="no. cases ct", xlab="age", lty=2)
lines(df[,"ct1"]~df[,"a"], type="l", ylab="no. cases ct", xlab="age", lty=3)
legend("topright", c("3rd trimester", "2nd trimester", "1st trimester"), lty=1:3)
