## ggplots
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(patchwork)

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
ggsave(filename = "plots/pop_size.png", width = 6, height = 6,
       dpi=600, units = "in", family = "Times")

#######################
## FoI profile plots ##
#######################
FoI_plots  <- vector("list", length=length(countries))
foi_recent <- vector("list", length=length(foi)) #list to store recent foi changes (1980-2030)
country_names    <- levels(countries)
country_names[7] <- "Iran"
min_year <- 1980
max_year <- 2030

for(k in 1:length(countries)){
  
  #set country
  k <- 2
  pars$country <- countries[k]
  
  #subset seroprevalence data
  fitting_data <- subset(df, country == pars$country)
  
  #no. of years between 1st & last data time points
  pars$tdiff   <- max(fitting_data$year) - min(fitting_data$year)
  
  #read in posterior distribution
  post_dist    <- readRDS(file = paste("posteriors/", pars$country, "/", 
                                       "posteriors_", pars$country, "_t", pars$temporal_foi, 
                                       "_", "a", pars$age_foi, ".RDS", sep=""))
  
  #store parameter values
  lambda0_vec <- post_dist$lambda
  beta_vec    <- post_dist$beta 
  tau_vec    <- post_dist$tau
  
  #set burn-in & time
  pars$burnin <- 2000
  tspan <- seq(1, pars$burnin + pars$tdiff + pars$years_forecast, 1)
  foi   <- vector("list", length=length(lambda0_vec)) #list to store FoI profile
  
  #store FoI for each parameter set
  for(i in 1:length(lambda0_vec)){
    
    lambda0 <- lambda0_vec[i]
    beta    <- beta_vec[i]
    tau     <- tau_vec[i]
    
    threshold <- pars$burnin - tau #set time of FoI change
    
    for(j in 1:length(tspan)){
      
      time <- tspan[j]
      
      if (time < threshold) {
        foi[[i]][j] <- lambda0 
        
      } else if (time >= threshold) {
        yr_rate <- (1 - beta) / ((pars$burnin + pars$tdiff) - threshold)  #define yearly rate of decline
        t_current <- time - threshold  #time difference
        foi[[i]][j] <- lambda0 * (1-(yr_rate * t_current)) 
        
        if(foi[[i]][j] <= 0){ # to avoid foi < 0 when forecasting
          foi[[i]][j] <- 0
        }
      }
      
    }
  }
  
  ## Plot all the FoI lines from the posterior
  
  ## First, calculate years from model time
  # min sampling year of the data
  data_min_year <- min(fitting_data$year)
  
  # equivalent model timepoint is burn-in, so calculate equivalent years by...
  time_diff <- data_min_year - pars$burnin
  years     <- tspan + time_diff
  
  # create dataset for 1980 to 2030
  year_indices <- which(years == 1980) : which(years == 2030)
  
  for(i in 1:length(foi_recent)){
    foi_recent[[i]] <- foi[[i]][year_indices]
  }
  
  foi_df[[k]] <- data.frame(
    "year"  = rep(years [ which(years == 1980) : which(years == 2030) ], length(lambda0_vec)), 
    "foi"   = unlist(foi_recent), 
    "group" = rep(1:600, each=length(1980:2030))
  )
  
  FoI_plots[[k]] <- ggplot(data=foi_df[[k]], aes(year, foi, group=group)) + 
    ggtitle(country_names[k]) + 
    geom_line(alpha=0.05) +
    # annotate("rect", xmin=min(fitting_data$year), xmax=max(fitting_data$year), #years with data
    #          ymin=0, ymax=Inf, alpha=0.3, fill=ribbonColour[1]) +
    xlab("Year") + 
    ylab(expression(lambda(t))) + 
    scale_x_continuous(expand = c(0,0), breaks = c(min_year, 2005, max_year)) + 
    scale_y_continuous(expand = c(0,0), n.breaks = 3) +
    theme(plot.margin=unit(c(rep(.5,4)),"cm")) + 
    theme_light(base_size = 12, base_line_size = 0, base_family = "Times") + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    theme(plot.margin=grid::unit(c(0, 0.5, 0, 0),"cm"))
  
}

## Plot & save multipanel FoI plot
wrap_plots(FoI_plots, nrow=4, ncol=3) + 
  plot_annotation(
    tag_levels = list(c('(a)', '(b)', '(c)', '(d)', '(e)', '(f)', 
                        '(g)', '(h)', '(i)', '(j)', '(k)')))

#PDF
ggsave(filename = "plots/foi_multipanel.pdf",
       device = cairo_pdf, height = 8, width = 8, units = "in")

# PNG
ggsave(filename = "plots/foi_multipanel.png",
       dpi=600, height = 8, width = 8, units = "in")
