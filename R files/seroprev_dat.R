################################################################
################## Global seroprevalence data ##################
################################################################

###################
## Load packages ##
###################
library(dplyr)
library(binom)

###############
## Load data ##
###############
data <- read.csv("data/seroprev_global.csv", fileEncoding="UTF-8-BOM") #to avoid erroneous labelling of col 1

#####################################################
## Clean data & apply inclusion/exclusion criteria ##
#####################################################
# minimum number of samples per country
min_samples <- 4 

# clean data
data <- data %>%
  filter(exclude=="n") %>%              # exclude rows with epidemiological biases or other errors
  filter(!is.na(first_sample)) %>%      # remove rows with no info on sampling year
  filter(method != "" | !is.na(method)) # only include studies with specified immunoassay method

# calculate median sampling year
for (i in 1:nrow(data)) {
  data$year[i] <- floor(median(seq(data$first_sample[i],data$last_sample[i], 1)))
}

# include countries with minimum no. of observations (defined above)
data <- data %>%
  group_by(country) %>%
  filter(n() >= min_samples)

# create separate data frames for countries, group by year
countries <- sort(unique(data$country))
data_list <- vector(mode = "list", length = length(countries))

## Store data for all countries in a list
for (i in 1:length(countries)) {
  data_list[[i]] <- 
    data %>%
    filter(country==countries[i]) %>% 
    filter(!is.na(year)) %>%  #omit data with missing sampling year
    arrange(year) %>%
    summarise(
      prev=k/n, #calculate seroprevalence
      year=year,  #sampling year
      year_published=year_published, #publishing year
      k=k,  #no. positive tests
      n=n,  #total no. tests
      country=countries[i], #country from which samples came
      method=method,  #immunoassay method 
      method2=method2,  #second immunoassay method (if applicable)
      method3=method3,  #third immunoassay method (if applicable)
      prop_method1=prop_method1,  #proportion of study's time spent using first immunoassay method
      prop_method2=prop_method2,  #proportion of study's time spent using second immunoassay method
      prop_method3=prop_method3) %>%  #proportion of study's time spent using third immunoassay method
    as.data.frame()
}

## transform data list into dataframe
df <- do.call(rbind.data.frame, data_list)
df$country <- as.factor(df$country)

## Calculate 95% CIs
cis <- binom.confint(x=df$k, n=df$n, conf.level=0.95, methods="exact")

## Store lower and upper estimates
df$ci_lo <- cis$lower
df$ci_up <- cis$upper

## Save data as R data file
saveRDS(df, "data/global_data.rds")
