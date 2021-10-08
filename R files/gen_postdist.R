#####################################################
## Generate, save and plot posterior distributions ##
#####################################################

##################
## Load scripts ##
##################
source("R files/setparms.R")
source("R files/model.R")   
source("R files/diagnostics.R")
source("R files/funcs.R")

###################
## Load packages ##
###################
library(ggplot2)
library(bayestestR)
library(ggpubr)
library(Cairo)
library(patchwork)
library(RColorBrewer)

################################################
## Read in parameter sets & likelihood values ##
################################################
niter <- 1000
nsim  <- 60
npars <- 3

# Save posterior distributions for later plotting
post.lambda <- post.beta <- post.tau <- vector("list", length=length(countries))
posteriors  <- vector("list", length=length(countries)) #to store summaries of posterior distributions

for (j in 1:length(countries)) {
  pars$country      <- countries[j]
  out_likpar        <- data.frame(matrix(ncol=npars+1, nrow=niter*nsim))
  names(out_likpar) <- c("log.lambda0", "log.beta", "log.tau", "likelihood")
  counter           <- seq(1, (nsim*niter), by = 1/nsim)  # used in for loop to pick correct parliks_ file
  
  #works for multiples of 10
  i_seq <- seq(1, niter*nsim, by=nsim)
  
  for (i in i_seq) {
    if(i==1){
      out_likpar[i:(i+nsim-1),] <- readRDS(file = paste("mod_output/", pars$country, "/new_fit", "/parliks_", pars$country, "_t", pars$temporal_foi, "_a", pars$age_foi, "_", i, ".Rdata", sep = ""))
    
      } else if (i > 1 & i < (niter*nsim)) {
      out_likpar[i:(i+nsim-1),] <- readRDS(file = paste("mod_output/", pars$country,  "/new_fit", "/parliks_", pars$country, "_t", pars$temporal_foi, "_a", pars$age_foi, "_", i-(i-counter[i]), ".Rdata", sep = ""))
    }
  }
  
  
  colnames(out_likpar) <- c("log.lambda0", "log.beta",  "log.tau", "likelihood")
  
  # sum(is.na(out_likpar$likelihood))  #check if NAs in likelihood
  
  
  ###########################################
  ## Summarise the posterior distributions ##
  ###########################################
  
  #priors
  prior.lambda <- out_likpar$log.lambda0
  prior.beta   <- out_likpar$log.beta
  prior.tau    <- out_likpar$log.tau
  
  #posteriors
  sorted_lik       <- head(out_likpar[order(out_likpar$likelihood),], nrow(out_likpar)*0.01) #take top 1%
  post.lambda[[j]] <- sorted_lik$log.lambda0
  post.beta[[j]]   <- sorted_lik$log.beta
  post.tau[[j]]    <- sorted_lik$log.tau
  
  ## Medians & 95% CIs
  # lambda0
  x<-describe_posterior(post.lambda[[j]], centrality = "median")
  post.median.lambda <- exp(x$Median)
  ci.low.lambda      <- exp(x$CI_low)
  ci.high.lambda     <- exp(x$CI_high)
  
  # beta
  x<-describe_posterior(post.beta[[j]], centrality = "median")
  post.median.beta <- exp(x$Median)
  ci.low.beta      <- exp(x$CI_low)
  ci.high.beta     <- exp(x$CI_high)
  
  # tau
  x <- describe_posterior(post.tau[[j]], centrality = "median")
  post.median.tau <- exp(x$Median)
  ci.low.tau      <- exp(x$CI_low)
  ci.high.tau     <- exp(x$CI_high)
  
  #### Summary of the posteriors #####
  posteriors[[j]] <- 
    data.frame(par = c("lambda0", "beta", "tau"), 
               med = c(round(post.median.lambda, 3), round(post.median.beta, 3), round(post.median.tau, 0)),
               lower_ci = c(round(ci.low.lambda, 3), round(ci.low.beta, 3), round(ci.low.tau, 0)), 
               upper_ci = c(round(ci.high.lambda, 3), round(ci.high.beta, 3), round(ci.high.tau, 0))
    )
  
  
  ############################################
  ## Plot prior and posterior distributions ##
  ############################################
  if(.Platform$OS.type == "windows") { # set Times New Roman font on Windows
    # library(extrafont)
    # font_import()
    # loadfonts(device = "win")
    windowsFonts(Times=windowsFont("TT Times New Roman")) 
  }
  
  ## Lambda0
  
  #dataframe
  dat <- data.frame(
    dens    = c(prior.lambda, post.lambda[[j]]), 
    ci_low  = ci.low.lambda,
    ci_high = ci.high.lambda,
    lines   = c(rep("a", length(prior.lambda)), rep("b", length(post.lambda[[j]])))
  )
  
  #plot
  if (exists("plot_lambda0") == F) {  #only create if list not in existence
    plot_lambda0 <- vector("list", length=length(countries))
  }
  
  plot_lambda0[[which(pars$country == countries)]] <- ggplot(
    dat, aes(x = exp(dens), fill = lines)) + 
    geom_density(alpha = 0.5) +
    geom_vline(xintercept = post.median.lambda) +                   # mode
    geom_vline(xintercept = ci.low.lambda, linetype='dotted') +   # lower ci
    geom_vline(xintercept = ci.high.lambda, linetype='dotted') +  # upper ci
    theme_light(base_size = 12, base_line_size = 0, base_family = "Times") + 
    theme(legend.position = "none", panel.grid = element_blank(), axis.ticks = element_blank()) +
    labs(fill = "") + 
    xlab(expression(lambda[0])) + 
    ylab("Density") + 
    theme(plot.margin=unit(c(rep(1,4)),"cm"))
  
  
  ## beta
  
  #dataframe
  dat <- data.frame(dens = c(prior.beta, post.beta[[j]]), 
                    ci_low = ci.low.beta,
                    ci_high = ci.high.beta,
                    lines = c(rep("a", length(prior.beta)), rep("b", length(post.beta[[j]]))))
  
  
  #plot
  if (exists("plot_beta") == F) {  #only create if list not in existence
    plot_beta <- vector("list", length=length(countries))
  }
  
  plot_beta[[which(pars$country == countries)]] <- ggplot(
    dat, aes(x = exp(dens), fill = lines)) + geom_density(alpha = 0.5) + 
    geom_vline(xintercept = post.median.beta) +                   # mode
    geom_vline(xintercept = ci.low.beta, linetype='dotted') +   # lower ci
    geom_vline(xintercept = ci.high.beta, linetype='dotted') +  # upper ci
    theme_light(base_size = 12, base_line_size = 0, base_family = "Times") + 
    theme(legend.position = "none", panel.grid = element_blank(), axis.ticks = element_blank()) +
    labs(fill = "") + 
    xlab(expression(beta)) + 
    ylab("Density") + 
    theme(plot.margin=unit(c(rep(1,4)),"cm"))
  
  
  ## tau
  
  #dataframe
  dat <- data.frame(dens = c(prior.tau, post.tau[[j]]), 
                    ci_low = ci.low.tau,
                    ci_high = ci.high.tau,
                    lines = c(rep("a", length(prior.tau)), rep("b", length(post.tau[[j]]))))
  
  #plot
  if (exists("plot_tau") == F) {  #only create if list not in existence
    plot_tau <- vector("list", length=length(countries))
  }
  
  plot_tau[[which(pars$country == countries)]] <- ggplot(
    dat, aes(x = exp(dens), fill = lines)) + geom_density(alpha = 0.5) + 
    geom_vline(xintercept = post.median.tau) +                   # mode
    geom_vline(xintercept = ci.low.tau, linetype='dotted') +   # lower ci
    geom_vline(xintercept = ci.high.tau, linetype='dotted') +  # upper ci 
    theme_light(base_size = 12, base_line_size = 0, base_family = "Times") + 
    theme(legend.position = "none", panel.grid = element_blank(), axis.ticks = element_blank()) +
    labs(fill = "") + 
    xlab(expression(tau)) + 
    ylab("Density") + 
    theme(plot.margin=unit(c(rep(1,4)),"cm"))
  
  
}

#####################################
## SAVE the posterior distribution ##
#####################################
post <- vector("list", length=length(countries))
for(j in 1:length(countries)){
  
  pars$country <- countries[j]
  
  post[[j]] <- data.frame(lambda = exp(post.lambda[[j]]), beta = exp(post.beta[[j]]), 
                          tau = exp(post.tau[[j]]))
  
  saveRDS(post[[j]], file = paste("posteriors/", pars$country, "/", 
                                  "posteriors_", pars$country, "_t", pars$temporal_foi, 
                                  "_", "a", pars$age_foi, ".RDS", sep=""))
}


#############################################################################
## Multipanel plot of priors vs. posteriors (one parameter, all countries) ##
#############################################################################

#lambda
multipanel_lambda <- wrap_plots(plot_lambda0, nrow=4, ncol=3) +
  plot_annotation(tag_levels = list(c('(a)', '(b)', '(c)', '(d)', '(e)', '(f)', '(g)', '(h)', '(i)', '(j)', '(k)'))) &
  theme(plot.margin=unit(c(rep(0.2,4)),"cm"))& 
  scale_y_continuous(n.breaks = 4)

#beta
multipanel_beta <- wrap_plots(plot_beta, nrow=4, ncol=3) +
  plot_annotation(tag_levels = list(c('(a)', '(b)', '(c)', '(d)', '(e)', '(f)', '(g)', '(h)', '(i)', '(j)', '(k)'))) &
  theme(plot.margin=unit(c(rep(0.2,4)),"cm"))& 
  scale_y_continuous(n.breaks = 4)

#tau
multipanel_tau <- wrap_plots(plot_tau, nrow=4, ncol=3) +
  plot_annotation(tag_levels = list(c('(a)', '(b)', '(c)', '(d)', '(e)', '(f)', '(g)', '(h)', '(i)', '(j)', '(k)'))) &
  theme(plot.margin=unit(c(rep(0.2,4)),"cm")) & 
  scale_y_continuous(n.breaks = 4)

#Save
# PDF
ggsave(multipanel_lambda, filename = "plots/lambda_multipanel.pdf",
       device = cairo_pdf, height = 8, width = 8, units = "in")

ggsave(multipanel_beta, filename = "plots/beta_multipanel.pdf",
       device = cairo_pdf, height = 8, width = 8, units = "in")

ggsave(multipanel_tau, filename = "plots/tau_multipanel.pdf",
       device = cairo_pdf, height = 8, width = 8, units = "in")

# PNG
ggsave(multipanel_lambda, filename = "plots/lambda_multipanel.png",
       dpi=600, height = 8, width = 8, units = "in")

ggsave(multipanel_beta, filename = "plots/beta_multipanel.png",
       dpi=600, height = 8, width = 8, units = "in")

ggsave(multipanel_tau, filename = "plots/tau_multipanel.png",
       dpi=600, height = 8, width = 8, units = "in")


########################################################################
## Dot and whisker plot of posterior medians & 95% credible intervals ##
########################################################################
## Label country
for(i in 1:length(countries)){
  posteriors[[i]]$country <- countries[i] 
}

## Put data list into one dataframe
df <- do.call(rbind.data.frame, posteriors)
levels(df$country)[7] <- "Iran"

## Create list to store plots for each parameter
plot_dw <- vector("list", length=3)

## Create colourblind friendly custom color scale
myColors    <- brewer.pal(11,"PuOr")
myColors[6] <- brewer.pal(11,"PRGn")[5]  #make pale colour bolder
names(myColors) <- levels(df$country)    #link colours to countries
colScale <- scale_colour_manual(name = "plotcols",values = myColors)

## Lambda0 posteriors plot ##
p <- ggplot(data = subset(df, par == "lambda0"), aes(x=par, y=med, color = country)) +
  geom_point(size = 2, position=position_dodge(width=1)) +
  geom_errorbar(
    aes(ymin = lower_ci, ymax = upper_ci),
    width = 0.1,
    position=position_dodge(width=1)) +
  theme_light(base_size = 12, base_line_size = 0, base_family = "Times") + 
  coord_flip()

plot_dw[[1]] <- p + 
  theme(legend.position = "bottom", legend.direction= "vertical", legend.title = element_blank(), 
        panel.grid = element_blank(), axis.text.y=element_blank(),
        axis.ticks.y = element_blank()) + 
  labs(x = "", y = expression(lambda[0])) + 
  guides(colour = guide_legend(reverse=T)) +
  colScale


## beta posteriors plot ##
p <- ggplot(data = subset(df, par == "beta"), aes(x=par, y=med, color = country)) +
  geom_hline(yintercept=1, linetype="dotted") +
  geom_point(size = 2, position=position_dodge(width=1)) +
  geom_errorbar(
    aes(ymin = lower_ci, ymax = upper_ci),
    width = 0.1,
    position=position_dodge(width=1)) +
  theme_light(base_size = 12, base_line_size = 0, base_family = "Times") + 
  coord_flip()

plot_dw[[2]] <- p + 
  theme(legend.position = "none",
        panel.grid = element_blank(), axis.text.y=element_blank(),
        axis.ticks.y = element_blank()) + 
  labs(x = "", y = expression(beta)) + 
  colScale


## tau posteriors plot ##
p <- ggplot(data = subset(df, par == "tau"), aes(x=par, y=med, color = country)) +
  geom_point(size = 2, position=position_dodge(width=1)) +
  geom_errorbar(
    aes(ymin = lower_ci, ymax = upper_ci),
    width = 0.1,
    position=position_dodge(width=1)) +
  theme_light(base_size = 12, base_line_size = 0, base_family = "Times") + 
  coord_flip()

plot_dw[[3]] <- p + 
  theme(legend.position = "none",
        panel.grid = element_blank(), axis.text.y=element_blank(),
        axis.ticks.y = element_blank()) + 
  labs(x = "", y = expression(tau)) + 
  colScale


## Calculate annual change in FoI for each country ##

#read in seroprevalence data
df <- readRDS("data/global_data.rds")

#set burn-in & threshold
pars$burnin <- 2000

#definte yearly rate vector (length of posterior distribution)
yr_rate <- vector("list", length=length(countries))

for(i in 1:length(countries)){
  
  #set country
  pars$country <- countries[i]
  
  #subset seroprevalence data
  fitting_data <- subset(df, country == pars$country)
  
  #no. of years between 1st & last data time points
  pars$tdiff   <- max(fitting_data$year) - min(fitting_data$year)
  
  #read in posterior distribution
  post_dist    <- readRDS(file = paste("posteriors/", pars$country, "/", 
                                       "posteriors_", pars$country, "_t", pars$temporal_foi, 
                                       "_", "a", pars$age_foi, ".RDS", sep=""))
  
  #save single parameter set
  lambda0 <- post_dist$lambda
  beta    <- post_dist$beta 
  tau     <- post_dist$tau
  
  #calculate annual % change in FoI
  threshold       <- pars$burnin - tau
  yr_rate[[i]] <- ((1 - beta) / ((pars$burnin + pars$tdiff) - threshold)*100)
  
  #correct sign, so positive numbers equate to increases, whereas negatives equate to declines
  yr_rate[[i]] <- yr_rate[[i]] * -1 
}

## Calculate medians & 95% credible intervals
med.gradient <- ci.low.gradient <- ci.high.gradient <- vector("numeric", length=length(countries))
for(i in 1:length(countries)){
  x <- describe_posterior(yr_rate[[i]], centrality = "median")
  med.gradient    [i] <- x$Median
  ci.low.gradient [i] <- x$CI_low
  ci.high.gradient[i] <- x$CI_high
}

## gradient posteriors plot ##
df <- data.frame(par = rep("gradient", 11), "gradient" = med.gradient, "ci_low_gradient" = ci.low.gradient, 
                 "ci_high_gradient" = ci.high.gradient, country=countries)
levels(df$country)[7] <- "Iran"

p <- ggplot(data = df, aes(x=par, y=gradient, color = country)) +
  geom_hline(yintercept=0, linetype="dotted") +
  geom_point(size = 2, position=position_dodge(width=1)) +
  geom_errorbar(
    aes(ymin = ci_low_gradient, ymax = ci_high_gradient),
    width = 0.1,
    position=position_dodge(width=1)) +
  theme_light(base_size = 12, base_line_size = 0, base_family = "Times") + 
  coord_flip()

gradient_plot <- p + 
  theme(legend.position = "none",
        panel.grid = element_blank(), axis.text.y=element_blank(),
        axis.ticks.y = element_blank()) + 
  labs(x = "", y = "Annual FoI change (%)") + 
  colScale


############################
## Plot & save multipanel ##
############################

##  Make legend a ggplot object
my_legend <- get_legend(plot_dw[[1]])
my_legend <- as_ggplot(my_legend)

## Remove legend from lambda plot
plot_dw[[1]] <- plot_dw[[1]] + theme(legend.position = "none")

## Arrange multipanel plot of lambda0 + beta + tau
ggarrange(plot_dw[[1]], plot_dw[[2]], plot_dw[[3]], 
          my_legend, font.label = list(face="plain", family="Times"),
          labels = c("(a)", "(b)", "(c)"), hjust=.05
          )

## Save
# PDF
ggsave(filename = "plots/posteriors_dwplot.pdf",
       device = cairo_pdf, height = 8, width = 8, units = "in")

# PNG
ggsave(filename = "plots/posteriors_dwplot.png",
       height = 8, width = 8, units = "in", dpi=600)


## Arrange multipanel plot of lambda0 + annual FoI change
ggarrange(plot_dw[[1]], my_legend, gradient_plot, 
          nrow=1, font.label = list(face="plain", family="Times"),
          labels = c("(a)", "", "(b)"), hjust=.05
)

## Save
# PNG
ggsave(filename = "plots/posteriors_FoIchange_dwplot.png",
       height = 6, width = 6, units = "in", dpi=600)

