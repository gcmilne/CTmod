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
nsim <- 60
npars <- 3

# Save posterior distributions for later plotting
post.lambda <- post.shape <- post.tdecline <- vector("list", length=length(countries))
posteriors <- vector("list", length=length(countries)) #to store summaries of posterior distributions

for(j in 1:length(countries)){
  pars$country <- countries[j]
  out_likpar <- data.frame(matrix(ncol=npars+1, nrow=niter*nsim))
  names(out_likpar) <- c("log.lambda0", "log.shape", "log.tdecline", "likelihood")
  counter <- seq(1, (nsim*niter), by = 1/nsim)  # used in for loop to pick correct parliks_ file
  
  #works for multiples of 10
  i_seq <- seq(1, niter*nsim, by=nsim)
  
  for(i in i_seq){
    if(i==1){
      out_likpar[i:(i+nsim-1),] <- readRDS(file = paste("mod_output/", pars$country, "/new_fit", "/parliks_", pars$country, "_t", pars$temporal_foi, "_a", pars$age_foi, "_", i, ".Rdata", sep = ""))
    } else if (i > 1 & i < (niter*nsim)){
      out_likpar[i:(i+nsim-1),] <- readRDS(file = paste("mod_output/", pars$country,  "/new_fit", "/parliks_", pars$country, "_t", pars$temporal_foi, "_a", pars$age_foi, "_", i-(i-counter[i]), ".Rdata", sep = ""))
    }
  }
  
  
  colnames(out_likpar) <- c("log.lambda0", "log.shape",  "log.tdecline", "likelihood")
  
  # sum(is.na(out_likpar$likelihood))  #check if NAs in likelihood
  
  
  ###########################################
  ## Summarise the posterior distributions ##
  ###########################################
  
  #priors
  prior.lambda   <- out_likpar$log.lambda0
  prior.shape    <- out_likpar$log.shape
  prior.tdecline <- out_likpar$log.tdecline
  
  #posteriors
  sorted_lik         <- head(out_likpar[order(out_likpar$likelihood),], nrow(out_likpar)*0.01) #take top 1%
  post.lambda[[j]]   <- sorted_lik$log.lambda0
  post.shape[[j]]    <- sorted_lik$log.shape
  post.tdecline[[j]] <- sorted_lik$log.tdecline
  
  ## Medians & 95% CIs
  # lambda0
  x<-describe_posterior(post.lambda[[j]], centrality = "median")
  post.median.lambda <- exp(x$Median)
  ci.low.lambda    <- exp(x$CI_low)
  ci.high.lambda   <- exp(x$CI_high)
  
  # shape
  x<-describe_posterior(post.shape[[j]], centrality = "median")
  post.median.shape <- exp(x$Median)
  ci.low.shape    <- exp(x$CI_low)
  ci.high.shape   <- exp(x$CI_high)
  
  # tdecline
  x<-describe_posterior(post.tdecline[[j]], centrality = "median")
  post.median.tdecline <- exp(x$Median)
  ci.low.tdecline    <- exp(x$CI_low)
  ci.high.tdecline   <- exp(x$CI_high)
  
  #### Summary of the posteriors #####
  posteriors[[j]] <- 
    data.frame(par = c("lambda0", "shape", "tdecline"), 
               med = c(round(post.median.lambda, 3), round(post.median.shape, 3), round(post.median.tdecline, 0)),
               lower_ci = c(round(ci.low.lambda, 3), round(ci.low.shape, 3), round(ci.low.tdecline, 0)), 
               upper_ci = c(round(ci.high.lambda, 3), round(ci.high.shape, 3), round(ci.high.tdecline, 0))
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
    dens = c(prior.lambda, post.lambda[[j]]), 
    ci_low = ci.low.lambda,
    ci_high = ci.high.lambda,
    lines = c(rep("a", length(prior.lambda)), rep("b", length(post.lambda[[j]])))
  )
  
  #plot
  if (exists("lambda0") == F) {  #only create if list not in existence
    lambda0 <- vector("list", length=length(countries))
  }
  
  lambda0[[which(pars$country == countries)]] <- ggplot(
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
  
  
  ## shape
  
  #dataframe
  dat <- data.frame(dens = c(prior.shape, post.shape[[j]]), 
                    ci_low = ci.low.shape,
                    ci_high = ci.high.shape,
                    lines = c(rep("a", length(prior.shape)), rep("b", length(post.shape[[j]]))))
  
  
  #plot
  if (exists("shape") == F) {  #only create if list not in existence
    shape <- vector("list", length=length(countries))
  }
  
  shape[[which(pars$country == countries)]] <- ggplot(
    dat, aes(x = exp(dens), fill = lines)) + geom_density(alpha = 0.5) + 
    geom_vline(xintercept = post.median.shape) +                   # mode
    geom_vline(xintercept = ci.low.shape, linetype='dotted') +   # lower ci
    geom_vline(xintercept = ci.high.shape, linetype='dotted') +  # upper ci
    theme_light(base_size = 12, base_line_size = 0, base_family = "Times") + 
    theme(legend.position = "none", panel.grid = element_blank(), axis.ticks = element_blank()) +
    labs(fill = "") + 
    # scale_x_continuous(breaks = seq(0, 0.8, 0.4), limits = c(0, 0.8), expand = c(0, 0)) +
    # scale_y_continuous(breaks = seq(0, 1.6, 0.75), limits = c(0, 1.6), expand = c(0, 0)) +
    # xlab("Shape") + 
    xlab(expression(beta)) + 
    ylab("Density") + 
    theme(plot.margin=unit(c(rep(1,4)),"cm"))
  
  
  ## tdecline
  
  #dataframe
  dat <- data.frame(dens = c(prior.tdecline, post.tdecline[[j]]), 
                    ci_low = ci.low.tdecline,
                    ci_high = ci.high.tdecline,
                    lines = c(rep("a", length(prior.tdecline)), rep("b", length(post.tdecline[[j]]))))
  
  #plot
  if (exists("tdecline") == F) {  #only create if list not in existence
    tdecline <- vector("list", length=length(countries))
  }
  
  tdecline[[which(pars$country == countries)]] <- ggplot(
    dat, aes(x = exp(dens), fill = lines)) + geom_density(alpha = 0.5) + 
    geom_vline(xintercept = post.median.tdecline) +                   # mode
    geom_vline(xintercept = ci.low.tdecline, linetype='dotted') +   # lower ci
    geom_vline(xintercept = ci.high.tdecline, linetype='dotted') +  # upper ci 
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
  
  post[[j]] <- data.frame(lambda = exp(post.lambda[[j]]), shape = exp(post.shape[[j]]), 
                          tdecline = exp(post.tdecline[[j]]))
  
  saveRDS(post[[j]], file = paste("posteriors/", pars$country, "/", "new_fit/", 
                                  "posteriors_", pars$country, "_t", pars$temporal_foi, 
                                  "_", "a", pars$age_foi, ".RDS", sep=""))
}


#############################################################################
## Multipanel plot of priors vs. posteriors (one parameter, all countries) ##
#############################################################################

#lambda
multipanel_lambda   <- wrap_plots(lambda0, nrow=4, ncol=3) +
  plot_annotation(tag_levels = list(c('(a)', '(b)', '(c)', '(d)', '(e)', '(f)', '(g)', '(h)', '(i)', '(j)', '(k)'))) &
  theme(plot.margin=unit(c(rep(0.2,4)),"cm"))& 
  scale_y_continuous(n.breaks = 4)

#shape
multipanel_shape    <- wrap_plots(shape, nrow=4, ncol=3) +
  plot_annotation(tag_levels = list(c('(a)', '(b)', '(c)', '(d)', '(e)', '(f)', '(g)', '(h)', '(i)', '(j)', '(k)'))) &
  theme(plot.margin=unit(c(rep(0.2,4)),"cm"))& 
  scale_y_continuous(n.breaks = 4)

#tdecline
multipanel_tdecline    <- wrap_plots(tdecline, nrow=4, ncol=3) +
  plot_annotation(tag_levels = list(c('(a)', '(b)', '(c)', '(d)', '(e)', '(f)', '(g)', '(h)', '(i)', '(j)', '(k)'))) &
  theme(plot.margin=unit(c(rep(0.2,4)),"cm")) & 
  scale_y_continuous(n.breaks = 4)

#Save
# PDF
ggsave(multipanel_lambda, filename = "plots/lambda_multipanel.pdf",
       device = cairo_pdf, height = 8, width = 8, units = "in")

ggsave(multipanel_shape, filename = "plots/shape_multipanel.pdf",
       device = cairo_pdf, height = 8, width = 8, units = "in")

ggsave(multipanel_tdecline, filename = "plots/tdecline_multipanel.pdf",
       device = cairo_pdf, height = 8, width = 8, units = "in")

# PNG
ggsave(multipanel_lambda, filename = "plots/lambda_multipanel.png",
       dpi=600, height = 8, width = 8, units = "in")

ggsave(multipanel_shape, filename = "plots/shape_multipanel.png",
       dpi=600, height = 8, width = 8, units = "in")

ggsave(multipanel_tdecline, filename = "plots/tdecline_multipanel.png",
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

## Create custom color scale
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


## Shape posteriors plot ##
p <- ggplot(data = subset(df, par == "shape"), aes(x=par, y=med, color = country)) +
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


## Tdecline posteriors plot ##
p <- ggplot(data = subset(df, par == "tdecline"), aes(x=par, y=med, color = country)) +
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


############################
## Plot & save multipanel ##
############################

##  Make legend a ggplot object
my_legend <- get_legend(plot_dw[[1]])
my_legend <- as_ggplot(my_legend)

## Remove legend from lambda plot
plot_dw[[1]] <- plot_dw[[1]] + theme(legend.position = "none")

## Arrange multipanel plot
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
