###############################################################
## Script to generate, save and plot posterior distributions ##
###############################################################

##################
## Load scripts ##
##################
source("R files/setparms.R")
# source("R files/seroprev_dat.R")
source("R files/model.R")   
source("R files/diagnostics.R")   

###################
## Load packages ##
###################
library(deSolve)
library(ggplot2)
library(bayestestR)
library(ggpubr)
library(Cairo)

###############################
## Load age profile function ##
###############################
getit <- function(time) {
  row <- which(abs(sol[,"time"]-time)==min(abs(sol[,"time"]-time)))
  df <- sol[row,-1]
  S  <- df[1:pars$agrps]
  I  <- df[(pars$agrps+1):(2*pars$agrps)]
  Im <- df[(2*pars$agrps+1):(3*pars$agrps)]
  pI <- df[(3*pars$agrps+1):(4*pars$agrps)]
  dprev <- df[(4*pars$agrps+1):(5*pars$agrps)]
  seroconv1 <- df[(5*pars$agrps+1):(6*pars$agrps)]
  seroconv2 <- df[(6*pars$agrps+1):(7*pars$agrps)] 
  seroconv3 <- df[(7*pars$agrps+1):(8*pars$agrps)] 
  matAb1 <- df[(8*pars$agrps+1):(9*pars$agrps)] 
  matAb2 <- df[(9*pars$agrps+1):(10*pars$agrps)]
  matAb3 <- df[(10*pars$agrps+1):(11*pars$agrps)]
  ct1 <- df[(11*pars$agrps+1):(12*pars$agrps)]
  ct2 <- df[(12*pars$agrps+1):(13*pars$agrps)]
  ct3 <- df[(13*pars$agrps+1):(14*pars$agrps)]
  Na <- df[(14*pars$agrps+1):(15*pars$agrps)]
  age <- pars$age  
  out <- data.frame(a=age, S=S, I=I,Im=Im, pI=pI, dprev=dprev, 
                    seroconv1=seroconv1, seroconv2=seroconv2, seroconv3=seroconv3, 
                    matAb1=matAb1, matAb2=matAb2, matAb3=matAb3, 
                    ct1=ct1, ct2=ct2, ct3=ct3, Na=Na)
  ## remove last age category for aesthetics
  # out <- out[-nrow(out),]
}


################################################
## Read in parameter sets & likelihood values ##
################################################
# niter <- 120                      # no. jobs to submit
# tot_iter <- 60000                 # total number of samples to run
# nsim <- ceiling(tot_iter / niter) # no. iterations on each for loop


if (pars$country == "Brazil" | pars$country == "Burkina Faso" | pars$country == "Cameroon") {
  niter <- 100
  
} else if (pars$country == "Iran (Islamic Republic of)") {
  niter <- 999
  
} else if (pars$country == "India") {
  niter <- 998
  
} else if (pars$country == "China" | pars$country == "Ethiopia" | pars$country == "Italy" | 
           pars$country == "Saudi Arabia" | pars$country == "Turkey" | pars$country == "United Kingdom") {
  niter <- 1000  # no. jobs to submit
}


if (pars$country == "Brazil" | pars$country == "Burkina Faso" | pars$country == "Cameroon") {
  nsim  <- 600

} else if (pars$country == "China" | pars$country == "Ethiopia" | pars$country == "Italy" | 
           pars$country == "Iran (Islamic Republic of)" | pars$country == "India" |
           pars$country == "Saudi Arabia" | pars$country == "Turkey" | pars$country == "United Kingdom") {
  nsim  <- 60
  
}

npars <- 3

out_likpar <- data.frame(matrix(ncol=npars+1, nrow=niter*nsim))
names(out_likpar) <- c("log.lambda0", "log.shape", "log.tdecline", "likelihood")
counter <- seq(1, (nsim*niter), by = 1/nsim)  # used in for loop to pick correct parliks_ file

#works for multiples of 10

## To do: read in differently for countries with missing iterations in the middle (e.g. India) ##
i_seq <- seq(1, niter*nsim, by=nsim)
for(i in i_seq){
  if(i==1){
    out_likpar[i:(i+nsim-1),] <- readRDS(file = paste("mod_output/", pars$country, "/parliks_", pars$country, "_t", pars$temporal_foi, "_a", pars$age_foi, "_", i, ".Rdata", sep = ""))
  } else if (i > 1 & i < (niter*nsim)){
    out_likpar[i:(i+nsim-1),] <- readRDS(file = paste("mod_output/", pars$country, "/parliks_", pars$country, "_t", pars$temporal_foi, "_a", pars$age_foi, "_", i-(i-counter[i]), ".Rdata", sep = ""))
  }
}

colnames(out_likpar) <- c("log.lambda0", "log.shape",  "log.tdecline", "likelihood")


sum(is.na(out_likpar$likelihood))


###########################################
## Summarise the posterior distributions ##
###########################################

#priors
prior.lambda   <- out_likpar$log.lambda0
prior.shape    <- out_likpar$log.shape
prior.tdecline <- out_likpar$log.tdecline

#posteriors
sorted_lik    <- head(out_likpar[order(out_likpar$likelihood),], nrow(out_likpar)*0.01) #take top 1%
post.lambda   <- sorted_lik$log.lambda0
post.shape    <- sorted_lik$log.shape
post.tdecline <- sorted_lik$log.tdecline

# SAVE the posterior distribution
post <- data.frame(post.lambda, post.shape, post.tdecline)
saveRDS(post, file = paste("posteriors/", pars$country, "/", "posteriors_", pars$country, "_t", pars$temporal_foi, "_", "a", pars$age_foi, ".RDS", sep=""))

## Medians & 95% CIs
# lambda0
x<-describe_posterior(post.lambda, centrality = "median")
post.median.lambda <- exp(x$Median)
ci.low.lambda    <- exp(x$CI_low)
ci.high.lambda   <- exp(x$CI_high)

# shape
x<-describe_posterior(post.shape, centrality = "median")
post.median.shape <- exp(x$Median)
ci.low.shape    <- exp(x$CI_low)
ci.high.shape   <- exp(x$CI_high)

# tdecline
x<-describe_posterior(post.tdecline, centrality = "median")
post.median.tdecline <- exp(x$Median)
ci.low.tdecline    <- exp(x$CI_low)
ci.high.tdecline   <- exp(x$CI_high)


############################################
## Plot prior and posterior distributions ##
############################################

## Lambda0

#dataframe
dat <- data.frame(dens = c(prior.lambda, post.lambda), 
                  ci_low = ci.low.lambda,
                  ci_high = ci.high.lambda,
                  lines = c(rep("a", length(prior.lambda)), rep("b", length(post.lambda))))

#plot
lambda0 <- ggplot(dat, aes(x = exp(dens), fill = lines)) + 
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = post.median.lambda) +                   # mode
  geom_vline(xintercept = ci.low.lambda, linetype='dotted') +   # lower ci
  geom_vline(xintercept = ci.high.lambda, linetype='dotted') +  # upper ci
  theme_light(base_size = 12, base_line_size = 0, base_family = "Times") + 
  theme(legend.position = "none", panel.grid = element_blank(), axis.ticks = element_blank()) +
  labs(fill = "") + 
  xlab(expression(lambda)) + 
  ylab("Density") + 
  theme(plot.margin=unit(c(rep(1,4)),"cm"))


## shape

#dataframe
dat <- data.frame(dens = c(prior.shape, post.shape), 
                  ci_low = ci.low.shape,
                  ci_high = ci.high.shape,
                  lines = c(rep("a", length(prior.shape)), rep("b", length(post.shape))))


#plot
shape <- ggplot(dat, aes(x = exp(dens), fill = lines)) + geom_density(alpha = 0.5) + 
  geom_vline(xintercept = post.median.shape) +                   # mode
  geom_vline(xintercept = ci.low.shape, linetype='dotted') +   # lower ci
  geom_vline(xintercept = ci.high.shape, linetype='dotted') +  # upper ci
  theme_light(base_size = 12, base_line_size = 0, base_family = "Times") + 
  theme(legend.position = "none", panel.grid = element_blank(), axis.ticks = element_blank()) +
  labs(fill = "") + 
  # scale_x_continuous(breaks = seq(0, 0.8, 0.4), limits = c(0, 0.8), expand = c(0, 0)) +
  # scale_y_continuous(breaks = seq(0, 1.6, 0.75), limits = c(0, 1.6), expand = c(0, 0)) +
  # xlab("Shape") + 
  xlab(expression(mu)) + 
  ylab("Density") + 
  theme(plot.margin=unit(c(rep(1,4)),"cm"))


## tdecline

#dataframe
dat <- data.frame(dens = c(prior.tdecline, post.tdecline), 
                  ci_low = ci.low.tdecline,
                  ci_high = ci.high.tdecline,
                  lines = c(rep("a", length(prior.tdecline)), rep("b", length(post.tdecline))))

#plot
tdecline <-ggplot(dat, aes(x = exp(dens), fill = lines)) + geom_density(alpha = 0.5) + 
  geom_vline(xintercept = post.median.tdecline) +                   # mode
  geom_vline(xintercept = ci.low.tdecline, linetype='dotted') +   # lower ci
  geom_vline(xintercept = ci.high.tdecline, linetype='dotted') +  # upper ci 
  theme_light(base_size = 12, base_line_size = 0, base_family = "Times") + 
  theme(legend.position = "none", panel.grid = element_blank(), axis.ticks = element_blank()) +
  labs(fill = "") + 
  xlab(expression(tau)) + 
  ylab("Density") + 
  theme(plot.margin=unit(c(rep(1,4)),"cm"))

## Multipanel plot of priors vs. posteriors
ggarrange(lambda0, shape, tdecline, ncol=2, nrow=2, 
          labels=c("(a)", "(b)", "(c)"), font.label=list(size=12, family="Times"), hjust = -2)

ggsave(filename = paste("plots/", pars$country, "/", pars$country, "_t", pars$temporal_foi, "_a", pars$age_foi, "_posteriors.pdf", sep=""), width = 6, height = 6,
units = "in", family = "Times")


####################################
#### Summary of the posteriors #####
####################################
posteriors <- data.frame(par = c("lambda0", "shape", "tdecline"), 
                         med = c(round(post.median.lambda, 3), round(post.median.shape, 3), round(post.median.tdecline, 0)),
                         lower_ci = c(round(ci.low.lambda, 3), round(ci.low.shape, 3), round(ci.low.tdecline, 0)), 
                         upper_ci = c(round(ci.high.lambda, 3), round(ci.high.shape, 3), round(ci.high.tdecline, 0))
)


