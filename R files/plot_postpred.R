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
ct_all <- readRDS("data/ct_predictions.rds")
prev_all <- readRDS("data/prev_predictions.rds")


################
## Make plots ##
################

## Prevalence ##

# make list to store all plots from different countries
if (!exists("prev_allyears")) {  #only create if list not in existence
  prev_allyears <- vector("list", length=length(countries))
}

#Colour for geom_ribbon
# display.brewer.pal(8, "Pastel1")
ribbonColour <- brewer.pal(8, "Pastel2")[c(3, 7)]

prev_allyears[[which(countries == pars$country)]] <- ggplot(
  data=prev_all[[which(countries == pars$country)]]) + 
  geom_line(aes(x=time, y=mod_prev), size=0.3) + 
  geom_ribbon(aes(x=time, ymin = mod_prev_low, ymax = mod_prev_up), alpha=0.5, fill=ribbonColour[2]) +
  annotate("rect", xmin=min(fitting_data$year), xmax=max(fitting_data$year), #years with data
           ymin=0, ymax=Inf, alpha=0.3, fill=ribbonColour[1]) +
  xlab("Year") + 
  ylab("Seroprevalence") + 
  scale_x_continuous(expand = c(0,0), limits = c(1900, 2030), n.breaks = 3) + 
  scale_y_continuous(expand = c(0,0)) +
  theme(plot.margin=unit(c(rep(.5,4)),"cm")) + 
  theme_light(base_size = 12, base_line_size = 0, base_family = "Times") + 
  theme(legend.position = "none") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

## Calculate where to place inset within main plot area
space_below_plot <- min(prev_all[[which(countries == pars$country)]]$mod_prev_low[which(years == 1900):which(years == 1975)])
space_above_plot <- 1 - (max(prev_all[[which(countries == pars$country)]]$mod_prev_up[which(years == 1900):which(years == 1975)]))

if (space_below_plot > space_above_plot) { 
  
  ymin_val <- 0.05
  ymax_val <- space_below_plot * 0.92
  xmin_val <- 1900 + 10
  xmax_val <- min(fitting_data$year) - 10
  
  # make arrow
  arrow_vec <- data.frame(
    x1 = min(fitting_data$year),
    x2 = xmax_val, 
    y1 = prev_all[[which(countries == pars$country)]]$mod_prev_low[which(prev_all[[which(countries == pars$country)]]$time == min(fitting_data$year))], #find intersection
    y2 = ymax_val * 0.95
  )
  
} else if (space_below_plot < space_above_plot) { 
  
  ymin_val <- 1 - space_above_plot + (0.02 * space_above_plot)
  ymax_val <- 0.98
  xmin_val <- 1900 + 10
  xmax_val <- min(fitting_data$year) - 10
  
  # make arrow
  arrow_vec <- data.frame(
    x1 = min(fitting_data$year),
    x2 = xmax_val, 
    y1 = prev_all[[which(countries == pars$country)]]$mod_prev_up[which(prev_all[[which(countries == pars$country)]]$time == min(fitting_data$year))], #find intersection
    y2 = ymin_val
  )
  
  prev_allyears[[which(countries == pars$country)]] <- prev_allyears[[which(countries == pars$country)]] +
    scale_y_continuous(expand = c(0,0), limits=c(0,1))  # change main plot limits to allow inset to fit
  
}

## Prevalence inset graph
if (!exists("prev_inset")) {  #only create if list not in existence
  prev_inset <- vector("list", length=length(countries))
}

# Calculate no. breaks for y axis
if (ymax_val - ymin_val > 0.4) {
  num_breaks <- 4
} else if (ymax_val - ymin_val < 0.4) { 
  num_breaks <- 3
} 

if (pars$country == "Italy" | pars$country == "Turkey"){
  num_breaks <- 2
}

prev_inset[[which(countries == pars$country)]] <- ggplot(
  data=prev_fit[[which(countries == pars$country)]], aes(x=time, y=mod_prev)) + 
  geom_line(size=0.3) + 
  geom_ribbon(aes(ymin = mod_prev_low, ymax = mod_prev_up), alpha=0.5, fill = ribbonColour[2]) +
  geom_point(aes(y=dat_prev), col="grey", size=0.2) +
  geom_errorbar(aes(width=.3, ymin = dat_low, ymax = dat_up),  col="grey", size=0.2) +
  scale_x_continuous(breaks=c(min(fitting_data$year), max(fitting_data$year))) +
  coord_cartesian(xlim=c(min(fitting_data$year), max(fitting_data$year, 
                                                     ylim=c(min(prev_fit[[which(countries == pars$country)]]$dat_low), 
                                                            max(prev_fit[[which(countries == pars$country)]]$dat_up))))) + #avoids error bars getting clipped
  scale_y_continuous(n.breaks=num_breaks) +
  theme_light(base_size = 12, base_line_size = 0, base_family = "Times") + 
  theme(legend.position = "none", axis.title.x=element_blank(), axis.title.y=element_blank()) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme(plot.margin=grid::unit(c(rep(0.1,4)),"cm"))

  
## Combined main prevalence plot with inset
if (!exists("prev_combo")) {  #only create if list not in existence
  prev_combo <- vector("list", length=length(countries))
}

prev_combo[[which(countries == pars$country)]] <- prev_allyears[[which(countries == pars$country)]] + 
  annotation_custom(
    ggplotGrob(prev_inset[[which(countries == pars$country)]]), 
    xmin = xmin_val, xmax = xmax_val, ymin = ymin_val, ymax = ymax_val
  ) + 
  geom_segment(data = arrow_vec, aes(x = x1, y = y1,       #add arrow
                                     xend = x2, yend = y2),
               lineend = "round", 
               linejoin = "round",
               size = 0.3, 
               arrow = arrow(length = unit(0.04, "inches")),
               colour = "grey") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())



## CT incidence ##

# make list to store all plots from different countries
if (!exists("ct_allyears")) {  #only create if list not in existence
  ct_allyears <- vector("list", length=length(countries))
}

# Main CT incidence plot
if(pars$country != "China"){
  ct_allyears[[which(countries == pars$country)]] <- ggplot(
    data=ct_all[[which(countries == pars$country)]], aes(x=time, y=ct_rel)) + 
    geom_line(aes(y=ct_rel), size=.3) +                        # Overall CT cases
    geom_ribbon(aes(ymin = ct_rel_low, ymax = ct_rel_up), alpha=0.5, fill = ribbonColour[2]) +
    annotate("rect", xmin=min(fitting_data$year), xmax=max(fitting_data$year), #years with data
             ymin=0, ymax=Inf, alpha=0.3, fill=ribbonColour[1]) +
    xlab("Year") + 
    ylab("Incidence per 10 000 live births") + 
    scale_x_continuous(expand = c(0,0), limits = c(1900, 2030), n.breaks = 3) + 
    scale_y_continuous(expand = c(0,0), n.breaks = 5) + 
    theme(plot.margin=unit(c(rep(.5,4)),"cm")) + 
    theme_light(base_size = 12, base_line_size = 0, base_family = "Times") + 
    theme(legend.position = "none") + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
} else if(pars$country == "China"){
  ct_allyears[[which(countries == pars$country)]] <- ggplot(
    data=ct_all[[which(countries == pars$country)]], aes(x=time, y=ct_rel)) + 
    geom_line(aes(y=ct_rel), size=.3) +                        # Overall CT cases
    geom_ribbon(aes(ymin = ct_rel_low, ymax = ct_rel_up), alpha=0.5, fill = ribbonColour[2]) +
    annotate("rect", xmin=min(fitting_data$year), xmax=max(fitting_data$year), #years with data
             ymin=0, ymax=Inf, alpha=0.3, fill=ribbonColour[1]) +
    xlab("Year") + 
    ylab("Incidence per 10 000 live births") + 
    scale_x_continuous(expand = c(0,0), limits = c(1900, 2030), n.breaks = 3) + 
    scale_y_continuous(limits=c(0, 300), expand = c(0,0), n.breaks = 5) + 
    theme(plot.margin=unit(c(rep(.5,4)),"cm")) + 
    theme_light(base_size = 12, base_line_size = 0, base_family = "Times") + 
    theme(legend.position = "none") + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
}

# Save
# ggsave(filename = paste("plots/", pars$country, "/", pars$temporal_foi, "_CT_allyears", ".pdf", sep=""),
# device = cairo_pdf, height = 6, width = 6, units = "in")


# CT sequelae incidence inset plot
# make list to store all plots from different countries
if (!exists("ct_inset")) {  #only create if list not in existence
  ct_inset <- vector("list", length=length(countries))
}

# colours
# display.brewer.pal(8, "Pastel2") # choose colours

myColors <- brewer.pal(7,"Set2") #set colours to use

inset_line_size = 0.3

max_incidence <- max(map_dbl(ct_all[[which(countries == pars$country)]][,-(1:4)], max)) # find max sequelae incidence

ct_inset[[which(countries == pars$country)]] <- ggplot(
  data=ct_all[[which(countries == pars$country)]], aes(x=time, y=ct_rel)) + 
  geom_line(aes(y=chorioretinitis_first), size=inset_line_size, col = myColors[1]) +         # Chorioretinitis in first year
  geom_line(aes(y=chorioretinitis_later), size=inset_line_size, col = myColors[2]) +         # Chorioretinitis in later life
  geom_line(aes(y=intracranial_calcifications), size=inset_line_size, col = myColors[3]) +   # Intracranial calcifications
  geom_line(aes(y=cns_abnormalities), size=inset_line_size, col = myColors[4]) +             # CNS abnormalities
  geom_line(aes(y=foetal_loss), size=inset_line_size, col = myColors[5]) +                   # Foetal loss
  geom_line(aes(y=hydrocephalus), size=inset_line_size, col = myColors[6]) +                 # Hydrocephalus
  geom_line(aes(y=neonatal_death), size=inset_line_size, col = myColors[7]) +                # Neonatal death
  annotate("rect", xmin=min(fitting_data$year), xmax=max(fitting_data$year), #years with data
           # ymin=min(ct_all[[which(countries == pars$country)]]$ct_rel_low), 
           ymin=0, ymax=Inf, alpha=0.3, fill=ribbonColour[1]) +
  scale_x_continuous(expand = c(0,0), limits = c(1900, 2030), n.breaks=2) + 
  scale_y_continuous(limits=c(0, max_incidence+0.05*max_incidence), expand = c(0,0), n.breaks=3) + 
  theme(plot.margin=unit(c(rep(.5,4)),"cm")) + 
  theme_light(base_size = 12, base_line_size = 0, base_family = "Times") + 
  theme(legend.position = "none", axis.title.x=element_blank(), axis.title.y=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme(plot.margin=grid::unit(c(0,0,0,0),"cm"))


# ct_inset[[which(countries == pars$country)]]


# Combine main plot with inset plot
# make list to store all plots from different countries
if (!exists("ct_combo")) {  #only create if list not in existence
  ct_combo <- vector("list", length=length(countries))
}

if (pars$country == "Cameroon" | pars$country == "Ethiopia") {
  ymin_val <-  max(ct_all[[which(countries == pars$country)]]$ct_rel_up[which(years == 1900):which(years == 1975)]) * 1.1
  ymax_val <-  max(ct_all[[which(countries == pars$country)]]$ct_rel_up[which(years == 1900):which(years == 2025)]) * 0.90
  xmin_val <- 1910
  xmax_val <- 1985
  
} else if (pars$country == "Iran (Islamic Republic of)") { 
  ymin_val <- 0
  ymax_val <- min(ct_all[[which(countries == pars$country)]]$ct_rel_low[which(years == 1940):which(years == 2000)]) * 0.91
  xmin_val <- 1940
  xmax_val <- 2000
  
} else if (pars$country == "China") {
  ymin_val <- 150
  ymax_val <- 300 * 0.91
  xmin_val <- 1910
  xmax_val <- 1985
  
  } else if (pars$country != "Cameroon" | pars$country != "Ethiopia" | 
             pars$country != "Iran (Islamic Republic of)" | pars$country != "China") { 
  ymin_val <- 0
  ymax_val <- min(ct_all[[which(countries == pars$country)]]$ct_rel_low[which(years == 1900):which(years == 1975)]) * 0.90
  xmin_val <- 1910
  xmax_val <- 1985
  
}

## Combined CT incidence plot with sequelae inset plot
ct_combo[[which(countries == pars$country)]] <- ct_allyears[[which(countries == pars$country)]] + 
  annotation_custom(
    ggplotGrob(ct_inset[[which(countries == pars$country)]]), 
    xmin = xmin_val, xmax = xmax_val, ymin = ymin_val, ymax = ymax_val
  ) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


###################### 
## Multipanel plots ##
######################

## Prevalence, with inset
wrap_plots(prev_combo, nrow=4, ncol=3) + 
  plot_annotation(tag_levels = list(c('(a)', '(b)', '(c)', '(d)', '(e)', '(f)', 
                                      '(g)', '(h)', '(i)', '(j)', '(k)')))

# letters[which(pars$country == countries)]  #which letter corresponds to which country?

ggsave(filename = "plots/prev_inset_multipanel.pdf",
       device = cairo_pdf, height = 8, width = 8, units = "in")

## CT, with inset
wrap_plots(ct_combo, nrow=4, ncol=3) & 
  scale_x_continuous(expand = c(0,0), limits = c(1900, 2030), n.breaks = 3) & 
  ylab("Incidence") &  
  plot_annotation(tag_levels = list(c('(a)', '(b)', '(c)', '(d)', '(e)', '(f)', 
                                      '(g)', '(h)', '(i)', '(j)', '(k)')))

ggsave(filename = "plots/ct_sequelae_multipanel.pdf",
       device = cairo_pdf, height = 8, width = 8, units = "in")


#####################################
## Temporal change in CT incidence ##
#####################################
if (!exists("ct_change")) {  #only create if list not in existence
  ct_change <- vector("list", length=length(countries))
}

min_index <- max_index <- vector("list", length=length(countries))

# Create dataset with just past & present casting (excluding forecasts)
ct_past <- vector("list", length=length(countries))

for (i in 1:length(ct_past)) {
  ct_past[[i]]$time   <- ct_all[[i]]$time      [ which(ct_all[[i]]$time == 1900) : which(ct_all[[i]]$time == 2021)]
  ct_past[[i]]$ct     <- ct_all[[i]]$ct_rel    [ which(ct_all[[i]]$time == 1900) : which(ct_all[[i]]$time == 2021)]
  ct_past[[i]]$ct_low <- ct_all[[i]]$ct_rel_low[ which(ct_all[[i]]$time == 1900) : which(ct_all[[i]]$time == 2021)]
  ct_past[[i]]$ct_up  <- ct_all[[i]]$ct_rel_up [ which(ct_all[[i]]$time == 1900) : which(ct_all[[i]]$time == 2021)]
}

# Calculate min and max incidence for all countries
for(i in 1:length(countries)){
  ct_change[[i]]$baseline <- ct_past[[i]]$ct[1]  #incidence at year 1900
  ct_change[[i]]$max_ct   <- max(ct_past[[i]]$ct)#max incidence
  ct_change[[i]]$min_ct   <- min(ct_past[[i]]$ct)#min incidence
  
  min_index <- which(ct_past[[i]]$ct == ct_change[[i]]$min_ct) #find min index 
  ct_change[[i]]$min_ct_low <- ct_past[[i]]$ct_low[min_index]    #lower CI
  ct_change[[i]]$min_ct_up  <- ct_past[[i]]$ct_up[min_index]     #upper CI
  
  max_index <- which(ct_past[[i]]$ct == ct_change[[i]]$max_ct) #find max index 
  ct_change[[i]]$max_ct_low <- ct_past[[i]]$ct_low[max_index]    #lower CI
  ct_change[[i]]$max_ct_up  <- ct_past[[i]]$ct_up[max_index]     #upper CI
}

# Print ranges of incidence predictions
for(i in 1:length(countries)){
  print(paste(
    countries[i], ": ", round(ct_change[[i]]$min_ct, 1), 
    " (", round(ct_change[[i]]$min_ct_low, 1), ", ", round(ct_change[[i]]$min_ct_up, 1), ") ", 
    "– ", round(ct_change[[i]]$max_ct, 1),
    " (", round(ct_change[[i]]$max_ct_low, 1), ", ", round(ct_change[[i]]$max_ct_up, 1), ") ", 
    sep=""
    
  ))
}



# which countries have increasing CT?
for(i in 1:length(countries)){
  
  if (ct_change[[i]]$max_ct > ct_change[[i]]$baseline) {
    
    print(paste("Country: ", countries[i], "; ", 
                # "Baseline CT incidence: ", round(ct_change[[i]]$baseline, 2), 
                # "Max CT incidence: ", round(ct_change[[i]]$max_ct, 2), 
                "Fold-change=", round(ct_change[[i]]$max_ct / ct_change[[i]]$baseline, 1), 
                sep = ""))
    
  }
}


for (i in 1:length(countries)) {
  
  #if minimum CT incidence is at baseline
  if (ct_change[[i]]$baseline == ct_change[[i]]$min_ct) { 
    
    #calculate median percentage increase in incidence
    ct_change[[i]]$percent_change <- round ( 
      (ct_change[[i]]$max_ct / ct_change[[i]]$baseline) * 100, 2
    )
    
    # else if minimum CT incidence is post-baseline
  } else if (ct_change[[i]]$baseline != ct_change[[i]]$min_ct) {
    
    #calculate median percentage decrease in incidence
    ct_change[[i]]$percent_change <- round ( 
      (ct_change[[i]]$max_ct / ct_change[[i]]$baseline) * 100, 2
    )
    
  }
  
  
}
