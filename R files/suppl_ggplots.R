## ggplots
library(ggplot2)
library(ggpubr)
library(RColorBrewer)

#read in data
df <- readRDS("data/global_data.rds")

#### Plot of publishing vs. sampling year
  
#create custom color scale
library(RColorBrewer)
myColors <- brewer.pal(11,"PuOr")
myColors[6] <- brewer.pal(11,"PRGn")[5]   #make pale colours bolder
names(myColors) <- levels(df$country)
colScale <- scale_colour_manual(name = "plotcols",values = myColors)

#plot
p1 <- ggplot(data=df, aes(x=year_published, y=year, colour = country)) +
  geom_point() +
  scale_x_continuous(limits = c(1982, 2020), expand = c(0,0)) + 
  scale_y_continuous(limits = c(1982, 2020), expand = c(0,0)) + 
  ylab("Median sampling year") +
  xlab("Publishing year") +
  theme_light(base_size = 12, base_line_size = .5, base_family = "Times") +
  theme(legend.position = "none", axis.ticks = element_blank()) + 
  theme(plot.margin=unit(c(rep(.5,4)),"cm"))

#add colours
p1 + colScale

# save plot
# ggsave(filename = "plots/publishing_vs_sampling_year.png", width = 6, height = 6,
# units = "in", family = "Times")


#### Plots of demographic data

## change to get plots for mortality, fertility or population size

## make ggplots for each of the countries
p <- vector(mode = "list", length = length(countries))


for(i in 1:length(countries)){
  pars$country <- countries[i]
  
  # pars$country <- "Brazil"
  fitting_data <- subset(df, df$country == pars$country)
  
  # Load demographic data
  source("R files/demogdat.R")
  
  ### interpolated parameters from demographic data
  # Population size #
  f <- spline(age_pop, tot_pop, xout=pars$age)
  
  pars$Na <- f$y*(sum(tot_pop[1:which(age_pop==pars$amax-2.5)])/sum(f$y)) # select correct age specific population size depending on choice of pars$amax
  pars$N  <- sum(pars$Na) # total population size
  
  # Mortality rate #
  pars$d <- spline(x=age_mort, y=mort, xout=pars$age)$y # interpolated age specific per capita mortality rates
  pars$d[which(pars$d<0)] <- 0  # zero minus elements
  
  # Birth rate #
  pars$propfert <- approx(age_fert, prop_fert_age, xout=pars$age)$y
  pars$propfert[is.na(pars$propfert)] <- 0
  pars$propfert <- pars$propfert/sum(pars$propfert)
  
  # set data
  data <- data.frame(age = pars$age, Na = pars$Na)
  
  #save ggplots
  if(countries[i] != "Iran (Islamic Republic of)"){
    
    p[[i]] <- ggplot(data=data, aes(x=age, y=Na)) + 
      geom_line() + 
      # scale_y_continuous(limits = c(0, 0.15), breaks = seq(0, 0.15, 0.05)) +
      scale_x_continuous(expand = c(0,0), limits = c(0, 55), breaks = seq(0, 55, 25)) + 
      ylab("Population size") +
      xlab("Age (years)") +
      labs(title = paste( countries[i])) +
      theme_light(base_size = 9, base_line_size = 0, base_family = "Times") +
      theme(legend.position = "none", axis.ticks = element_blank()) + 
      theme(plot.margin=unit(c(rep(.5,4)),"cm"))
    
  } else if(countries[i] == "Iran (Islamic Republic of)"){
    
    p[[i]] <- ggplot(data=data, aes(x=age, y=Na)) + 
      geom_line() + 
      # scale_y_continuous(limits = c(0, 0.15), breaks = seq(0, 0.15, 0.05)) +
      scale_x_continuous(expand = c(0,0), limits = c(0, 55), breaks = seq(0, 55, 25)) + 
      ylab("Population size") +
      xlab("Age (years)") +
      labs(title = c(paste("Iran"))) +
      theme(legend.position = "none", axis.ticks = element_blank()) + 
      theme(plot.margin=unit(c(rep(.5,4)),"cm")) + 
      theme_light(base_size = 9, base_line_size = 0, base_family = "Times")
  
  }
}

wrap_plots(p)


# save plot
# ggsave(filename = "plots/pop_size.png", width = 6, height = 6,
#        units = "in", family = "Times")



####
## Temporal foi plot
####

plot(sol[,"time"], sol[,"pIt"], 'l')

res <- data.frame(t = sol[,"time"], prev =  sol[,"pIt"])

# create data for arrows showing key time points
arrow_x <- c(meta[1:15], rep(0, 187))
arrow_y <- rep(c(1.2, 0), times=c(15, 187))
arrow_yend <- rep(c(0.2, 0), times=c(15, 187))

ggplot(data=res, aes(x=t, y=prev)) + 
  geom_line() + 
  scale_x_continuous(expand = c(0,0), limits = c(min(time)-25, max(time)+25)) + 
  scale_y_continuous(expand = c(0,0), limits=c(0,0.05)) + 
  ylab("CT incidence") +
  xlab("Time") +
  
  ## to exemplify temporal parameters
  geom_vline(xintercept = pars$burnin, linetype = "dotted") + 
  geom_vline(xintercept = pars$burnin - exp(pars$log.tdecline), linetype = "dotted") + 
  geom_vline(xintercept = pars$burnin + pars$tdiff, linetype = "dotted") + 
  
  theme(plot.margin=unit(c(rep(.5,4)),"cm")) + 
  theme_light(base_size = 12, base_line_size = 0, base_family = "Times") + 
  theme(legend.position = "none", axis.ticks = element_blank(), 
        axis.text.x=element_blank(), axis.text.y=element_blank())


# stepwise foi decline
foi <-  c(rep(0.02, 100), rep(0.02*0.5, 100))
time <- seq(1, 200)

dat <- data.frame(t=time, foi=foi)

p_step <- ggplot(data=dat, aes(x=t, y=foi)) + 
  geom_line() + 
  # scale_x_continuous(expand = c(0,0), limits = c(min(time)-25, max(time)+25)) + 
  # scale_y_continuous(expand = c(0,0), limits=c(0,0.05)) + 
  ylab("Force of infection") +
  xlab("Time") +
  
  ## to exemplify temporal parameters
  # geom_vline(xintercept = pars$burnin, linetype = "dotted") + 
  # geom_vline(xintercept = pars$burnin - exp(pars$log.tdecline), linetype = "dotted") + 
  # geom_vline(xintercept = pars$burnin + pars$tdiff, linetype = "dotted") + 
  
  theme(plot.margin=unit(c(rep(.5,4)),"cm")) + 
  theme_light(base_size = 12, base_line_size = 0, base_family = "Times") + 
  theme(legend.position = "none", axis.ticks = element_blank(), 
        axis.text.x=element_blank(), axis.text.y=element_blank())


# linear foi decline
foi <-  c(rep(0.02, 100), seq(0.02 -  (0.02-0.01)/100 , 0.02*0.5, by = - (0.02-0.01)/100 ))
time <- seq(1, 200)

dat <- data.frame(t=time, foi=foi)

p_lin <- ggplot(data=dat, aes(x=t, y=foi)) + 
  geom_line() + 
  scale_x_continuous(expand = c(0,0)) +
  # scale_y_continuous(expand = c(0,0)) +
  ylab("Force of infection") +
  xlab("Time") +
  
  ## to exemplify temporal parameters
  # geom_vline(xintercept = pars$burnin, linetype = "dotted") + 
  # geom_vline(xintercept = pars$burnin - exp(pars$log.tdecline), linetype = "dotted") + 
  # geom_vline(xintercept = pars$burnin + pars$tdiff, linetype = "dotted") + 
  
  theme(plot.margin=unit(c(rep(.5,4)),"cm")) + 
  theme_light(base_size = 12, base_line_size = 0, base_family = "Times") + 
  theme(legend.position = "none", axis.ticks = element_blank(), 
        axis.text.x=element_blank(), axis.text.y=element_blank())


wrap_plots(p_step, p_lin, ncol=1)


setEPS()
postscript("plots/temporal_foi.eps", family="Times", width = 6, height = 6)
ggarrange(p_step, p_lin, ncol=1, nrow=2, labels=c("(a)", "(b)"), font.label=list(size=12, family="Times"), hjust = .02)



dev.off()


### Burn-in time plots
sol <- ode(y=y, times=time, func=age_si, parms=pars)

data <- data.frame(time = sol[,"time"], prev = sol[,"pIt"], ct = sol[,"ctt"])

# Seroprevalence over time
p1 <- ggplot(data=data, aes(x=time, y=prev)) + 
  geom_line() + 
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), limits=c(0,0.5)) +
  ylab("Seroprevalence") +
  xlab("Time") +
  
  geom_vline(xintercept = pars$burnin, linetype = "dotted") +
  
  theme_light(base_size = 12, base_line_size = 0, base_family = "Times") + 
  theme(legend.position = "none", axis.ticks = element_blank(), plot.margin=unit(c(rep(.5,4)),"cm"))

# CT incidence over time
p2 <- ggplot(data=data, aes(x=time, y=ct)) + 
  geom_line() + 
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), limits=c(0,15000)) +
  ylab("Incidence of CT") +
  xlab("Time") +
  
  geom_vline(xintercept = pars$burnin, linetype = "dotted") +
  
  theme_light(base_size = 12, base_line_size = 0, base_family = "Times") + 
  theme(legend.position = "none", axis.ticks = element_blank(), plot.margin=unit(c(rep(.5,4)),"cm"))



setEPS()
postscript("plots/burnin.eps", family="Times", width = 6, height = 6)
ggarrange(p1, p2, ncol=1, nrow=2, labels=c("(a)", "(b)"), font.label=list(size=12, family="Times"), hjust = .02)
dev.off()