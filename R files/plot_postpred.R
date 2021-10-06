############################################################
## Script to plot posterior predictions from each country ##
############################################################

#################
# Load packages #
#################
library(ggplot2)
library(binom)
library(patchwork)
library(RColorBrewer)
library(purrr)
library(RColorBrewer)

#############
# Load data #
#############
ct_all   <- readRDS("data/ct_predictions.rds")
prev_all <- readRDS("data/prev_predictions.rds")
prev_fit <- readRDS("data/prev_data.rds")

#############
# Set fonts #
#############
if(.Platform$OS.type == "windows") { # set Times New Roman font on Windows
  # library(extrafont)
  # font_import()
  # loadfonts(device = "win")
  windowsFonts(Times=windowsFont("TT Times New Roman")) 
}


######################################
## Set country socioeconomic status ##
######################################
# Below from World Bank (https://data.worldbank.org/country)
income_status <- list("high" = c("Italy", "Saudi Arabia", "United Kingdom"), 
                      "upper_middle" = c("Brazil", "China", "Turkey"),
                      "lower_middle" = c("Cameroon", "India", "Iran"),
                      "low" = c("Burkina Faso", "Ethiopia"))


################
## Make plots ##
################
min_year <- 1980 #Minimum year to display modelled estimates from
max_year <- 2030 #Maximum year to display modelled estimates until

## Prevalence ##

# make list to store all plots from different countries
if (!exists("prev_allyears")) {  #only create if list not in existence
  prev_allyears <- vector("list", length=length(countries))
}

#Colour for geom_ribbon
ribbonColour <- brewer.pal(8, "Pastel2")[c(3, 7)]

#Colours for income status
income_cols <- brewer.pal(4, "Blues")

levels(countries)[7] <- "Iran" ## Make Iran's name shorter for plotting

for (i in 1:length(countries)) {
  
  pars$country <- levels(countries)[i]  #set country
  
  ## store colour according to income status
  if (pars$country %in% income_status$high) {
    ses_colour <- income_cols[1]
  } else if (pars$country %in% income_status$upper_middle) {
    ses_colour <- income_cols[2]
  } else if (pars$country %in% income_status$lower_middle) {
    ses_colour <- income_cols[3]
  } else if (pars$country %in% income_status$low) {
    ses_colour <- income_cols[4]
  }
  
  prev_allyears[[i]] <- ggplot(
    data=prev_all[[i]]) + 
    ggtitle(levels(countries)[i]) + 
    geom_line(aes(x=time, y=mod_prev), size=0.3) + 
    geom_ribbon(aes(x=time, ymin = mod_prev_low, ymax = mod_prev_up), alpha=0.5, fill=ribbonColour[2]) +
    geom_point(data=prev_fit[[i]], 
               aes(x=time, y=dat_prev), col="grey", size=0.2) +
    geom_errorbar(data=prev_fit[[i]], 
                  aes(width=.3, x=time, ymin = dat_low, ymax = dat_up),  col="grey", size=0.2) +
    xlab("Year") + 
    ylab("Seroprevalence (%)") + 
    scale_x_continuous(expand = c(0,0), limits = c(min_year, max_year), breaks = c(min_year, 2005, max_year)) + 
    scale_y_continuous(expand = c(0,0), limits = c(0,100)) +
    theme(plot.margin=unit(c(rep(.5,4)),"cm")) + 
    theme_light(base_size = 12, base_line_size = 0, base_family = "Times") + 
    theme(legend.position = "none") + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    theme(plot.margin=grid::unit(c(0, 0.5, 0, 0),"cm"))
  
  prev_allyears[[i]] <- prev_allyears[[i]] + 
    theme(panel.background = element_rect(colour = ses_colour, 
                                          size = 2, linetype = "solid"))
}

## Plot & save multipanel prevalence plot
wrap_plots(prev_allyears, nrow=4, ncol=3) + 
  plot_annotation(
    tag_levels = list(c('(a)', '(b)', '(c)', '(d)', '(e)', '(f)', 
                        '(g)', '(h)', '(i)', '(j)', '(k)')))

#PDF
ggsave(filename = "plots/prev_multipanel.pdf",
       device = cairo_pdf, height = 8, width = 8, units = "in")

# PNG
ggsave(filename = "plots/prev_multipanel.png",
       dpi=600, height = 8, width = 8, units = "in")

## CT incidence ##

# make list to store all plots from different countries
if (!exists("ct_allyears")) {  #only create if list not in existence
  ct_allyears <- vector("list", length=length(countries))
}

# Main CT incidence plot
for (i in 1:length(countries)) {
  
  pars$country <- levels(countries)[i]  #set country
  
  ## store colour according to income status
  if (pars$country %in% income_status$high) {
    ses_colour <- income_cols[1]
  } else if (pars$country %in% income_status$upper_middle) {
    ses_colour <- income_cols[2]
  } else if (pars$country %in% income_status$lower_middle) {
    ses_colour <- income_cols[3]
  } else if (pars$country %in% income_status$low) {
    ses_colour <- income_cols[4]
  }
  
  if(pars$country != "China" & pars$country != "United Kingdom"){
    ct_allyears[[i]] <- ggplot(
      data=ct_all[[i]], aes(x=time, y=ct_rel)) + 
      ggtitle(levels(countries)[i]) + 
      geom_line(aes(y=ct_rel), size=.3) +  #Overall CT cases
      geom_ribbon(aes(ymin = ct_rel_low, ymax = ct_rel_up), alpha=0.5, fill = ribbonColour[2]) +
      annotate("rect", xmin=min(prev_fit[[i]]$time), xmax=max(prev_fit[[i]]$time), #years with data
               ymin=0, ymax=Inf, alpha=0.3, fill=ribbonColour[1]) +
      xlab("Year") + 
      ylab("Incidence per 10 000 live births") + 
      scale_x_continuous(expand = c(0,0), limits = c(min_year, max_year), breaks = c(min_year, 2005, max_year)) + 
      scale_y_continuous(expand = c(0,0), n.breaks = 5) + 
      theme(plot.margin=unit(c(rep(.5,4)),"cm")) + 
      theme_light(base_size = 12, base_line_size = 0, base_family = "Times") + 
      theme(legend.position = "none") + 
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
      theme(plot.margin=grid::unit(c(0, 0.5, 0, 0),"cm"))
    
    ct_allyears[[i]] <- ct_allyears[[i]] + 
      theme(panel.background = element_rect(colour = ses_colour, 
                                            size = 2, linetype = "solid"))
    
  } else if(pars$country == "China" | pars$country == "United Kingdom"){ #create space at top of graph for inset plots
    
    if (pars$country == "China") y_limits <- c(0, 200) else y_limits <- c(0,40)
    
    ct_allyears[[i]] <- ggplot(
      data=ct_all[[i]], aes(x=time, y=ct_rel)) + 
      ggtitle(levels(countries)[i]) + 
      geom_line(aes(y=ct_rel), size=.3) +   #Overall CT cases
      geom_ribbon(aes(ymin = ct_rel_low, ymax = ct_rel_up), alpha=0.5, fill = ribbonColour[2]) +
      annotate("rect", xmin=min(prev_fit[[i]]$time), xmax=max(prev_fit[[i]]$time), #years with data
               ymin=0, ymax=Inf, alpha=0.3, fill=ribbonColour[1]) +
      xlab("Year") + 
      ylab("Incidence per 10 000 live births") + 
      scale_x_continuous(expand = c(0,0), limits = c(min_year, max_year), breaks = c(min_year, 2005, max_year)) + 
      scale_y_continuous(limits=y_limits, expand = c(0,0), n.breaks = 5) + 
      theme(plot.margin=unit(c(rep(.5,4)),"cm")) + 
      theme_light(base_size = 12, base_line_size = 0, base_family = "Times") + 
      theme(legend.position = "none") + 
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
      theme(plot.margin=grid::unit(c(0, 0.5, 0, 0),"cm"))
    
    ct_allyears[[i]] <- ct_allyears[[i]] + 
      theme(panel.background = element_rect(colour = ses_colour, 
                                            size = 2, linetype = "solid"))
  }
}

# CT sequelae incidence inset plot
# make list to store all plots from different countries
if (!exists("ct_inset")) {  #only create if list not in existence
  ct_inset <- vector("list", length=length(countries))
}

# colours
myColors <- brewer.pal(7,"Set2") #set colours to use

inset_line_size = 0.3

for(i in 1:length(countries)){
  max_incidence <- max(map_dbl(ct_all[[i]][,-(1:4)], max)) # find max sequelae incidence
  
  ct_inset[[i]] <- ggplot(
    data=ct_all[[i]], aes(x=time, y=ct_rel)) + 
    geom_line(aes(y=chorioretinitis_first), size=inset_line_size, col = myColors[1]) +         # Chorioretinitis in first year
    geom_line(aes(y=chorioretinitis_later), size=inset_line_size, col = myColors[2]) +         # Chorioretinitis in later life
    geom_line(aes(y=intracranial_calcifications), size=inset_line_size, col = myColors[3]) +   # Intracranial calcifications
    geom_line(aes(y=cns_abnormalities), size=inset_line_size, col = myColors[4]) +             # CNS abnormalities
    geom_line(aes(y=foetal_loss), size=inset_line_size, col = myColors[5]) +                   # Foetal loss
    geom_line(aes(y=hydrocephalus), size=inset_line_size, col = myColors[6]) +                 # Hydrocephalus
    geom_line(aes(y=neonatal_death), size=inset_line_size, col = myColors[7]) +                # Neonatal death
    annotate("rect", xmin=min(prev_fit[[i]]$time), xmax=max(prev_fit[[i]]$time), #years with data
             ymin=0, ymax=Inf, alpha=0.3, fill=ribbonColour[1]) +
    scale_x_continuous(expand = c(0,0), limits = c(min_year, max_year), breaks=c(min_year, max_year)) + 
    scale_y_continuous(limits=c(0, max_incidence+0.05*max_incidence), expand = c(0,0), n.breaks=3) + 
    theme(plot.margin=unit(c(rep(.5,4)),"cm")) + 
    theme_light(base_size = 12, base_line_size = 0, base_family = "Times") + 
    theme(legend.position = "none", axis.title.x=element_blank(), axis.title.y=element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    theme(plot.margin=grid::unit(c(0,0,0,0),"cm"))
}

# Combine main plot with inset plot
# make list to store all plots from different countries
if (!exists("ct_combo")) {  #only create if list not in existence
  ct_combo <- vector("list", length=length(countries))
}

ymin_val <- ymax_val <- xmin_val <- xmax_val <- vector("numeric", length=length(countries))

for(i in 1:length(countries)){
  if (countries[i] == "China") {
    ymin_val[i] <- 60
    ymax_val[i] <- 200 * 0.91
    xmin_val[i] <- min_year + 2
    xmax_val[i] <- min_year + (max_year - min_year)/2 + 2
    
  } else if (countries[i] ==  "Cameroon" | countries[i] ==  "Ethiopia"){
    ymin_val[i] <- 4
    ymax_val[i] <- 61
    xmin_val[i] <- max_year - (max_year - min_year)/2 - 6
    xmax_val[i] <- max_year - 6
    
  } else if (countries[i] ==  "United Kingdom") {
    ymin_val[i] <- 18
    ymax_val[i] <- 38
    xmin_val[i] <- max_year - (max_year - min_year)/2 - 6
    xmax_val[i] <- max_year - 6
    
  } else if (countries[i] != "China" | countries[i] != "Cameroon" | 
             countries[i] != "Ethiopia" | countries[i] != "United Kingdom") { 
    ymin_val[i] <- 5
    ymax_val[i] <- min(ct_all[[i]]$ct_rel_low[which(years == min_year):which(years == 2007)]) * 0.90
    xmin_val[i] <- min_year + 2
    xmax_val[i] <- min_year + (max_year - min_year)/2 + 2
    
  }
}

## Combined CT incidence plot with sequelae inset plot
for(i in 1:length(countries)){
  ct_combo[[i]] <- ct_allyears[[i]] + 
    annotation_custom(
      ggplotGrob(ct_inset[[i]]), 
      xmin = xmin_val[i], xmax = xmax_val[i], ymin = ymin_val[i], ymax = ymax_val[i]
    ) + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
}

## Plot & save CT incidence estimates, with inset graph
wrap_plots(ct_combo, nrow=4, ncol=3) & 
  ylab("Incidence") &  
  plot_annotation(tag_levels = list(c('(a)', '(b)', '(c)', '(d)', '(e)', '(f)', 
                                      '(g)', '(h)', '(i)', '(j)', '(k)')))

ggsave(filename = "plots/ct_sequelae_multipanel.pdf",
       device = cairo_pdf, height = 8, width = 8, units = "in")

ggsave(filename = "plots/ct_sequelae_multipanel.png",
       height = 8, width = 8, units = "in", dpi=600)
