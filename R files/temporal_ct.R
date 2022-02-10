###############################################
## Calculate temporal change in CT incidence ##
###############################################

#########################
# Set working directory #
#########################
cluster <- "none"

##################
## Load scripts ##
##################
source("R files/setparms.R")

#############
# Load data #
#############
ct_all   <- readRDS("data/ct_predictions.rds")    #CT incidence estimates
prev_all <- readRDS("data/prev_predictions.rds")  #Seroprevalence estimates
prev_fit <- readRDS("data/prev_data.rds")         #Serorevalence data

## Create dataset (from min country-specific year of data sampling to 2030)
ct_past <- vector("list", length=length(countries))
min_year <- 1980  #minimum year
max_year <- 2030  #maximum year

for (i in 1:length(ct_past)) {
  
  # min_year <- min(prev_fit[[i]]$time)  #country-specific minimum sampling year
  
  ct_past[[i]]$time   <- ct_all[[i]]$time      [ which(ct_all[[i]]$time == min_year) : which(ct_all[[i]]$time == max_year)]
  ct_past[[i]]$ct     <- ct_all[[i]]$ct_rel    [ which(ct_all[[i]]$time == min_year) : which(ct_all[[i]]$time == max_year)]
  ct_past[[i]]$ct_low <- ct_all[[i]]$ct_rel_low[ which(ct_all[[i]]$time == min_year) : which(ct_all[[i]]$time == max_year)]
  ct_past[[i]]$ct_up  <- ct_all[[i]]$ct_rel_up [ which(ct_all[[i]]$time == min_year) : which(ct_all[[i]]$time == max_year)]
}

## Calculate min and max incidence for all countries
min_index <- max_index <- vector("list", length=length(countries))

if (!exists("ct_change")) {  #only create if list not in existence
  ct_change <- vector("list", length=length(countries))
}

for(i in 1:length(countries)){
  ct_change[[i]]$baseline <- ct_past[[i]]$ct[1]   #incidence at first sampling year
  ct_change[[i]]$max_ct   <- max(ct_past[[i]]$ct) #max incidence
  ct_change[[i]]$min_ct   <- min(ct_past[[i]]$ct) #min incidence
  
  min_index <- which(ct_past[[i]]$ct == ct_change[[i]]$min_ct)   #find min index 
  ct_change[[i]]$min_ct_low <- ct_past[[i]]$ct_low[min_index]    #lower CI
  ct_change[[i]]$min_ct_up  <- ct_past[[i]]$ct_up[min_index]     #upper CI
  
  max_index <- which(ct_past[[i]]$ct == ct_change[[i]]$max_ct)   #find max index 
  ct_change[[i]]$max_ct_low <- ct_past[[i]]$ct_low[max_index]    #lower CI
  ct_change[[i]]$max_ct_up  <- ct_past[[i]]$ct_up[max_index]     #upper CI
}

## which countries have increasing CT?
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
  
  # Print change in CT incidence as %
  print(paste(countries[i], round(ct_change[[i]]$percent_change - 100, 2), sep=", "))
}

#############################################
## Extract CT estimates for specific years ##
#############################################

## (used to compare our estimates to those of a previous study)

# Min:max sampling years that systematic reivew extracted data (Rostami et al 2019 PloS NTD)
years_list <- list(1997:2014, 
                   2004:2014, 
                   2009, 
                   1993:2017, 
                   2011:2015, 
                   2002:2017, 
                   1999:2017, 
                   1987:2014, 
                   1999:2016, 
                   2002:2016,
                   1992)

# Create lists to store indices & matched CT estimates
matched_ct <- indices <- vector("list", length=length(years_list))

# Fill lists 
for (i in 1:length(years_list)) {
  
  # Find indices
  indices[[i]] <- which(ct_all[[i]]$time == years_list[[i]][1]): which(ct_all[[i]]$time == years_list[[i]][length(years_list[[i]])])
  
  
  # Find CT incidence & uncertainties
  matched_ct[[i]]$ct      <- ct_all[[i]]$ct_rel    [indices[[i]]]
  matched_ct[[i]]$ct_low  <- ct_all[[i]]$ct_rel_low[indices[[i]]]
  matched_ct[[i]]$ct_high <- ct_all[[i]]$ct_rel_up [indices[[i]]]
  
}

# Calculate median CT incidence over whole period
point_estimates <- vector("list", length=length(matched_ct))

for (i in 1:length(point_estimates)) {
  
  point_estimates[[i]]$central <- round(median(matched_ct[[i]]$ct), 1)
  point_estimates[[i]]$low     <- round(median(matched_ct[[i]]$ct_low), 1)
  point_estimates[[i]]$up      <- round(median(matched_ct[[i]]$ct_high), 1)
  
}

# Display in useful format for table
for (i in 1:length(point_estimates)) {
  print(paste(countries[i], ": ", point_estimates[[i]]$central, " (", 
              point_estimates[[i]]$low, ", ", point_estimates[[i]]$up, ")", sep=""))
}
