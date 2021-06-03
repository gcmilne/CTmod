################################################################
################## Global seroprevalence data ##################
################################################################
library(dplyr)
library(ggplot2)
library(patchwork)
library(lme4)

data <- read.csv("data/seroprev_global.csv")
min_samples <- 4 #minimum number of longitudinal samples

# clean data
data <- data %>%
  filter(exclude=="n") %>%              # exclude rows with epidemiological biases or other errors
  filter(!is.na(first_sample)) %>%      # remove rows with no info on sampling year
  filter(method != "" | !is.na(method)) # only include studies with specified immunoassay method

# calculate prevalence
data$prev <- data$k/data$n

# calculate median sampling year
for(i in 1:nrow(data)){
  data$year[i] <- floor(median(seq(data$first_sample[i],data$last_sample[i], 1)))
}

# include countries with minimum of x observations
data <- data %>%
  group_by(country) %>%
  filter(n() >= min_samples)

# create separate dfs for countries, group by year & calculate weighted prevalence
countries <- sort(unique(data$country))
data_list <- vector(mode = "list", length = length(countries))

for(i in 1:length(countries)){
  data_list[[i]] <- 
    data %>%
    filter(country==countries[i]) %>% 
    filter(!is.na(year)) %>%
    arrange(year) %>%
    group_by(year) %>%
    summarise(w.prev = weighted.mean(prev, n), k=sum(k), n = sum(n), country=countries[i], 
              state=state, method=method, method2=method2, method3=method3, prop_method1=prop_method1,
              prop_method2=prop_method2, prop_method3=prop_method3) %>%
    # distinct(k, n, .keep_all = T) %>% #remove duplicate rows (where median sampling years are equal)
    as.data.frame()
}

## combine methods into one row, for matching median sampling years
f <- function(x) {
  x <- na.omit(x)
  if (length(x) > 0) paste(x,collapse=', ') else NA
}

res <- vector(mode = "list", length = length(countries))

for(i in 1:length(data_list)){
  res[[i]] <- data_list[[i]] %>% group_by(year) %>% summarise_all(list(f))
}

## correct country names to include only one instance in each row
for(i in 1:length(data_list)){
  res[[i]]$country <- countries[i]
}

## correct k to include only one instance in each row
for(i in 1:length(res)){
  for(j in 1:nrow(res[[i]])){
    x <- strsplit(res[[i]]$k, ", ")
    if(length(x[[j]]) > 1) {
      res[[i]]$k[j] <- x[[j]][1] #select first mention
    }
  }
}

## correct n to include only one instance in each row
for(i in 1:length(res)){
  for(j in 1:nrow(res[[i]])){
    x <- strsplit(res[[i]]$n, ", ")
    if(length(x[[j]]) > 1) {
      res[[i]]$n[j] <- x[[j]][1] #select first mention
    }
  }
}

## correct weighted seroprevalence to include only one instance in each row
for(i in 1:length(res)){
  for(j in 1:nrow(res[[i]])){
    x <- strsplit(res[[i]]$w.prev, ", ")
    if(length(x[[j]]) > 1) {
      res[[i]]$w.prev[j] <- x[[j]][1] #select first mention
    }
  }
}
  
## make altered columns numeric
for(i in 1:length(res)){
  res[[i]]$k <- as.numeric(res[[i]]$k)
  res[[i]]$n <- as.numeric(res[[i]]$n)
  res[[i]]$w.prev <- as.numeric(res[[i]]$w.prev)
}
  

## make ggplots for each of the countries
p <- vector(mode = "list", length = length(countries))

# for(j in 1:length(who_regions)){  ## if you want to order plots by region too
for(i in 1:length(countries)){
  p[[i]] <- ggplot(res[[i]], aes(x=year, y=w.prev)) +
    geom_point() +
    geom_smooth(method = "lm", formula = 'y~x') +
    xlab("Year") + ylab("Seroprevalence") +
    scale_y_continuous(breaks = seq(0,1, 0.2)) +
    labs(title = paste(countries[i])) +
    theme_light()
}

wrap_plots(p)


## put data list into one dataframe
df <- do.call(rbind.data.frame, res)
df$country <- as.factor(df$country)

# saveRDS(df, "data/global_data.rds")

####################
## Regional plots ##
####################

# # omit rows with no regional data
# for(i in 1:length(countries)){
#   data_list[[i]] <- data_list[[i]] %>%
#     filter(!is.na(year)) %>%
#     filter(!is.na(state))
# }
# 
# # store names of unique states
# states <- vector(mode = "list", length = length(countries))
# for(i in 1:length(countries)){
#   states[[i]] <- unique(data_list[[i]]$state)
# }
# 
# # find n times each unique state is present in dataset
# n_states <- table(data$state)
# 
# # store names of states if >n number of samples
# n <- 3
# 
# freq_states <- vector("numeric", length=length(n_states))
# 
# for(i in 1:length(n_states)){
#   if(n_states[i] >= n){
#     print(n_states[i])
#     freq_states[i] <- names(n_states[i])
#   }else{
#     freq_states[i] <- NA
#   }
# }
# 
# freq_states <- freq_states[!is.na(freq_states)] #omit na
# 
# ## store data & calculate weighted prevalence for these states
# regional_data <- vector(mode = "list", length = length(freq_states))
# 
# for(i in 1:length(freq_states)){
#   regional_data[[i]] <- data %>%
#     filter(!is.na(year)) %>%
#     filter(!is.na(state)) %>%
#     filter(state==freq_states[[i]]) %>%
#     arrange(year) %>%
#     group_by(year) %>%
#     summarise(w.prev = weighted.mean(prev, n), k=sum(k), n = sum(n), country=country, state=freq_states[[i]]) %>%
#     as.data.frame()
# }
# 
# # omit those with fewer than n observations
# for(i in 1:length(freq_states)){
#   if(length(unique(regional_data[[i]]$year)) < n){
#     regional_data[[i]] <- NA
#   }
# }
# 
# regional_data <- regional_data[!is.na(regional_data)] #omit na
# 
# # ggplot for each state
# p <- vector(mode = "list", length = length(regional_data))
# for(i in 1:length(regional_data)){
#   p[[i]] <- vector(mode = "list", length = length(regional_data[[i]]))
# }
# 
# for(i in 1:length(regional_data)){
#   p[[i]] <- ggplot(regional_data[[i]], aes(x=year, y=w.prev)) + 
#     geom_point() + 
#     geom_smooth(method = "lm", formula = 'y~x') +
#     scale_y_continuous(breaks = seq(0,1, 0.2)) +
#     labs(title = paste(c(freq_states[i]), regional_data[[i]]$country, sep=", ")) + 
#          xlab("Year") +
#          ylab("Seroprevalence")
# }
# 
# # wrap_plots(p)
