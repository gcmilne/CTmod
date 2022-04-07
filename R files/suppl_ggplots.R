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

###########################
## Set working directory ##
###########################
cluster <- "none"

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

## ggplots
data <- readRDS(data, file = "data/demographic_data_allcountries.RDS")

for(i in 1:length(countries)){
  
  # Add year of demographic data
  pars$country <- countries[i]                            #set country
  fitting_data <- subset(df, df$country == pars$country)  #subset data
  data[[i]]$year <- round(min(fitting_data$year)/5)*5
  
  #fertility plot
  p_fert[[i]] <- ggplot(data=data[[i]], aes(x=age, y=fert)) + 
    ggtitle(paste(c(levels(countries)[i]), " (", data[[i]]$year[1], ")", sep="")) + 
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
    ggtitle(paste(c(levels(countries)[i]), " (", data[[i]]$year[1], ")", sep="")) + 
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
p_fert[[7]] <- p_fert[[7]] + ggtitle(paste(c("Iran"), " (", data[[7]]$year[1], ")", sep=""))
p_mort[[7]] <- p_mort[[7]] + ggtitle(paste(c("Iran"), " (", data[[7]]$year[1], ")", sep=""))
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
# change back to linear FoI change
pars$temporal_foi  <- "linear"

foi_plots  <- foi_df <- vector("list", length=length(countries))
foi_recent <- vector("list", length=600) #list to store recent foi changes (1980-2030)
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
  
  foi_plots[[k]] <- ggplot(data=foi_df[[k]], aes(year, foi, group=group)) + 
    ggtitle(country_names[k]) + 
    geom_line(alpha=0.05) +
    xlab("Year") + 
    ylab(expression(lambda(t))) + 
    scale_x_continuous(expand = c(0,0), breaks = c(min_year, 2005, max_year)) + 
    scale_y_continuous(expand = c(0,0), n.breaks = 3) +
    theme(plot.margin=unit(c(rep(.5,4)),"cm")) + 
    theme_light(base_size = 12, base_line_size = 0, base_family = "Times") + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    theme(plot.margin=grid::unit(c(0, 0.5, 0, 0),"cm"))
  
}


## Calculate mean FoI profile for each country & plot on one graph
med_foi <- rep(list(setNames(data.frame(matrix(ncol=2, nrow=length(1980:2030))), c("year", "foi"))), length(countries))

for(i in 1:length(countries)){  #each country
  for(j in 1:length(1980:2030)){
    med_foi[[i]]$foi[j] <- median(foi_df[[i]]$foi[which(foi_df[[i]]$year == 1979+j)])
  }
  med_foi[[i]]$year <- 1980:2030
}

## transform list into dataframe
med_foi <- do.call(rbind.data.frame, med_foi)
med_foi$country <- rep(countries, each=length(1980:2030))

p <- ggplot(data=med_foi, aes(year, foi, group=country)) + 
  ggtitle("All countries") + 
  geom_line(alpha=0.05) +
  xlab("Year") + 
  ylab(expression(lambda(t))) + 
  scale_x_continuous(expand = c(0,0), breaks = c(min_year, 2005, max_year)) + 
  scale_y_continuous(expand = c(0,0), n.breaks = 3) +
  theme(plot.margin=unit(c(rep(.5,4)),"cm")) + 
  theme_light(base_size = 12, base_line_size = 0, base_family = "Times") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme(plot.margin=grid::unit(c(0, 0.5, 0, 0),"cm"))

## Plot & save multipanel FoI plot
foi_plots[[12]] <- p 

wrap_plots(foi_plots, nrow=4, ncol=3) + 
  plot_annotation(
    tag_levels = list(c('(a)', '(b)', '(c)', '(d)', '(e)', '(f)', 
                        '(g)', '(h)', '(i)', '(j)', '(k)', '(l)')))

#PDF
ggsave(filename = "plots/foi_multipanel.pdf",
       device = cairo_pdf, height = 8, width = 8, units = "in")

#PNG
ggsave(filename = "plots/foi_multipanel.png",
       dpi=600, height = 8, width = 8, units = "in")


############################################################
## Linear regression of publication year on sampling year ##
############################################################

## Estimate accuracy of publication year as proxy for sampling year
t.test(df$year_published, df$year, paired=T)

# Make colours
# create custom color scale
myColors        <- brewer.pal(11,"PuOr")
myColors[6]     <- brewer.pal(11,"PRGn")[5]   #make pale colour bolder
names(myColors) <- levels(df$country)
colScale <- scale_colour_manual(name = "plotcols",values = myColors)

# Run linear model
mod <- lm(year~year_published, data=df)

# Make y=x line to show what perfect condordance looks like
y_x <- min(df$year_published) : max(df$year_published)
y_x <- c(y_x, rep(NA, 100 - length(y_x))) # make vector same length as other data

# Plot
p1 <- ggplot(data=df, aes(x=year_published, y=year, colour = country)) +
  geom_point() +
  geom_line(aes(x=year_published, y=mod$fitted.values), col="darkgrey") + #regression line
  geom_line(aes(x=y_x, y=y_x), col="darkgrey", linetype="dashed") + #y=x line
  scale_x_continuous(limits = c(1982, 2020), expand = c(0,0)) + 
  scale_y_continuous(limits = c(1982, 2020), expand = c(0,0)) + 
  ylab("Median sampling year") +
  xlab("Publishing year") +
  theme_light(base_size = 12, base_line_size = .5, base_family = "Times") +
  theme(legend.position = "none", axis.ticks = element_blank()) + 
  theme(plot.margin=unit(c(rep(.5,4)),"cm"))

# Add colours & zoom plot
p1 + colScale + 
  coord_cartesian(
  xlim = c(1990, 2020),
  ylim = c(1982, 2020))

# Save plot
ggsave(filename = "plots/publishing_vs_sampling_year.png", width = 6, height = 6, 
       units = "in", dpi=600, family = "Times")
