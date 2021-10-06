## ggplots
library(ggplot2)
library(ggpubr)
library(RColorBrewer)

#############
# Set fonts #
#############
if(.Platform$OS.type == "windows") { # set Times New Roman font on Windows
  # library(extrafont)
  # font_import()
  # loadfonts(device = "win")
  windowsFonts(Times=windowsFont("TT Times New Roman")) 
}

#read in data
df <- readRDS("data/global_data.rds")

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


### Burn-in time plots
data <- data.frame(time = sol[,"time"], prev = sol[,"pIt"], ct = sol[,"ctt"])

# Seroprevalence over time
p1 <- ggplot(data=data, aes(x=time, y=prev)) + 
  geom_line() + 
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), limits=c(0,0.5)) +
  ylab("Seroprevalence") +
  xlab("Time") +
  
  geom_vline(xintercept = 2000, linetype = "dotted") +
  
  theme_light(base_size = 12, base_line_size = 0, base_family = "Times") + 
  theme(legend.position = "none", axis.ticks = element_blank(), plot.margin=unit(c(rep(.5,4)),"cm"))

# CT incidence over time
p2 <- ggplot(data=data, aes(x=time, y=ct)) + 
  geom_line() + 
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), limits=c(0,15000)) +
  ylab("CT cases") +
  xlab("Time") +
  
  geom_vline(xintercept = 2000, linetype = "dotted") +
  
  theme_light(base_size = 12, base_line_size = 0, base_family = "Times") + 
  theme(legend.position = "none", axis.ticks = element_blank(), plot.margin=unit(c(rep(.5,4)),"cm"))

# Save (eps)
setEPS()
postscript("plots/burnin.eps", family="Times", width = 6, height = 6)
ggarrange(p1, p2, ncol=1, nrow=2, labels=c("(a)", "(b)"), font.label=list(size=12, family="Times"), hjust = .02)
dev.off()

# Save (png)
ggsave("plots/burnin.png", family="Times", width = 6, height = 6, dpi=600)

