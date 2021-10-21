###############################################
## ggplots supplementary to the main results ##
###############################################

###################
## Load packages ##
###################
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(patchwork)
library(deSolve)

##################
## Load scripts ##
##################
source("R files/setparms.R")
source("R files/model.R")
source("R files/funcs.R")

###############
## Set fonts ##
###############
if(.Platform$OS.type == "windows") { # set Times New Roman font on Windows
  # library(extrafont)
  # font_import()
  # loadfonts(device = "win")
  windowsFonts(Times=windowsFont("TT Times New Roman")) 
}

##################
## Read in data ##
##################
df <- readRDS("data/global_data.rds")

###############################
## Plots of demographic data ##
###############################
## change to get plots for mortality, fertility or population size

## Make list to store ggplots for each of the countries
# (1) Fertility rate; (2) Mortality rate; (3) Population size
p_fert <- p_mort <-  p_pop <- vector(mode = "list", length = length(countries))
data   <- vector(mode = "list", length = length(countries))  #data object for plotting

## Set global parameters
pars$temporal_foi  <- "none"    #TEMPORAL FoI decline (options: "none", "stepwise" or "linear" )
time <- seq(1,2000)

## Run the model for each country & store output (only needs to be run once)
# for(i in 1:length(countries)){
# 
#   pars$country <- countries[i]
#   fitting_data <- subset(df, df$country == pars$country)
#   
#   # Load demographic data
#   source("R files/demogdat.R")
#   
#   ### interpolated parameters from demographic data
#   # Population size #
#   f <- spline(age_pop, tot_pop, xout=pars$age)
#   
#   pars$Na <- f$y*(sum(tot_pop[1:which(age_pop==pars$amax-2.5)])/sum(f$y)) # select correct age specific population size depending on choice of pars$amax
#   pars$N  <- sum(pars$Na) # total population size
#   
#   # Mortality rate #
#   pars$d <- spline(x=age_mort, y=mort, xout=pars$age)$y # interpolated age specific per capita mortality rates
#   pars$d[which(pars$d<0)] <- 0  # zero minus elements
#   
#   # Birth rate #
#   pars$propfert <- approx(age_fert, prop_fert_age, xout=pars$age)$y
#   pars$propfert[is.na(pars$propfert)] <- 0
#   pars$propfert <- pars$propfert/sum(pars$propfert)
#   
#   # Initial values (state variables) #
#   S0 <- I0 <- Im0 <- vector("numeric", length=pars$agrps)
#   
#   # Entire population initially susceptible
#   S0[1:length(S0)] <- pars$Na
#   y                <- c(S0, I0, Im0)  
#   
#   # Run model to calculate population size after burn-in
#   sol <- ode(y=y, times=time, parms=pars, func = age_si)  #get model solution
#   age_dist  <- getit(2000)                                #get age distribution @ burn-in
# 
#   # set data
#   data[[i]] <- data.frame(age = pars$age, 
#                      fert     = pars$propfert*100,  #fertility (as %)
#                      mort     = pars$d,             #mortality
#                      dat_na   = pars$Na,            #pop size (data)
#                      mod_na   = age_dist[,"Na"])    #pop size (model)
#   
# }
# 
# # save data
# saveRDS(data, file = "data/demographic_data_allcountries.RDS")

## Make ggplots
data <- readRDS(data, file = "data/demographic_data_allcountries.RDS")

for(i in 1:length(countries)){
  
  #fertility plot
  p_fert[[i]] <- ggplot(data=data[[i]], aes(x=age, y=fert)) + 
    ggtitle(levels(countries)[i]) + 
    geom_line() + 
    scale_x_continuous(expand = c(0,0), limits = c(0, 55), breaks = seq(0, 55, 25)) + 
    scale_y_continuous(n.breaks = 3) + 
    ylab("Pregnant (%)") +
    xlab("Age (years)") +
    theme(plot.margin=unit(c(rep(.5,4)),"cm")) + 
    theme_light(base_size = 12, base_line_size = 0, base_family = "Times") + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    theme(plot.margin=grid::unit(c(0, 0.5, 0, 0),"cm"))
  
  #mortality plot
  p_mort[[i]] <- ggplot(data=data[[i]], aes(x=age, y=mort)) + 
    ggtitle(levels(countries)[i]) + 
    geom_line() + 
    scale_x_continuous(expand = c(0,0), limits = c(0, 55), breaks = seq(0, 55, 25)) + 
    scale_y_continuous(n.breaks = 3) + 
    ylab("Mortality rate") +
    xlab("Age (years)") +
    theme(plot.margin=unit(c(rep(.5,4)),"cm")) + 
    theme_light(base_size = 12, base_line_size = 0, base_family = "Times") + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    theme(plot.margin=grid::unit(c(0, 0.5, 0, 0),"cm"))
  
  #population size plot
  #set minimum y axis value
  if ( min(data[[i]]$dat_na) < min(data[[i]]$mod_na) ) { 
    y_min <- min(data[[i]]$dat_na) - .05*min(data[[i]]$dat_na)
  } else {
    y_min <- min(data[[i]]$mod_na) - .05*min(data[[i]]$mod_na)
  }
  
  #set maximum y axis value
  if ( max(data[[i]]$dat_na) < max(data[[i]]$mod_na) ) { 
    y_max <- max(data[[i]]$mod_na) + .05*max(data[[i]]$mod_na)
  } else {
    y_max <- max(data[[i]]$dat_na) + .05*max(data[[i]]$dat_na)
  }
  
  p_pop[[i]] <- ggplot(data[[i]], aes(x=age, y=dat_na)) + 
    ggtitle(levels(countries)[i]) + 
    geom_bar(stat="identity", width=0.1, col="grey")+
    geom_line(aes(x=age, y=mod_na, color=fert)) +
    labs(x = "Age (years)", y = "Population size", color = "Pregnant (%)") +
    theme(plot.margin=unit(c(rep(.5,4)),"cm")) + 
    theme_light(base_size = 12, base_line_size = 0, base_family = "Times") + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    theme(plot.margin=grid::unit(c(0, 0.5, 0, 0),"cm")) +
    scale_color_viridis_c(option = "cividis")
  
}

## Shorten Iran's name
p_fert[[7]] <- p_fert[[7]] + ggtitle("Iran")
p_mort[[7]] <- p_mort[[7]] + ggtitle("Iran")
p_pop[[7]]  <- p_pop[[7]]  + ggtitle("Iran")

## Multipanel plots
## (1) Fertility rate ##
x <- wrap_plots(p_fert, nrow=4, ncol=3) &
  scale_y_continuous(n.breaks=5)

x + plot_annotation(
  tag_levels = list(c('(a)', '(b)', '(c)', '(d)', '(e)', '(f)', 
                      '(g)', '(h)', '(i)', '(j)', '(k)')))

#save
#PNG
ggsave(filename = "plots/fertility_rates.png", width = 8, height = 8,
       dpi=600, units = "in", family = "Times")

#PDF
ggsave(filename = "plots/fertility_rates.pdf",
       device = cairo_pdf, height = 8, width = 8, units = "in")

## (2) Mortality rate ##
x <- wrap_plots(p_mort, nrow=4, ncol=3) &
  scale_y_continuous(n.breaks=5)

x + plot_annotation(
  tag_levels = list(c('(a)', '(b)', '(c)', '(d)', '(e)', '(f)', 
                      '(g)', '(h)', '(i)', '(j)', '(k)')))

#save
#PNG
ggsave(filename = "plots/mortality_rates.png", width = 8, height = 8,
       dpi=600, units = "in", family = "Times")

#PDF
ggsave(filename = "plots/mortality_rates.pdf",
       device = cairo_pdf, height = 8, width = 8, units = "in")

## (3) Population size ##
x <- wrap_plots(p_pop, nrow=4, ncol=3) &
scale_y_continuous(n.breaks=5)

x + plot_annotation(
    tag_levels = list(c('(a)', '(b)', '(c)', '(d)', '(e)', '(f)', 
                        '(g)', '(h)', '(i)', '(j)', '(k)')))

#save
#PNG
ggsave(filename = "plots/pop_size.png", width = 10, height = 8,
       dpi=600, units = "in", family = "Times")

#PDF
ggsave(filename = "plots/pop_size.pdf",
       device = cairo_pdf, height = 10, width = 8, units = "in")


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
  tau_vec     <- post_dist$tau
  
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

#PNG
ggsave(filename = "plots/foi_multipanel.png",
       dpi=600, height = 8, width = 8, units = "in")
