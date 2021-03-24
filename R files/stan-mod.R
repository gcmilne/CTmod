#library("rstanarm")
library("ggplot2")
# library("rstan")
library("dplyr")
# install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
# install_cmdstan(overwrite=T)
library("cmdstanr")

#add option of using old Stan compiler in installation (avoid binary .exe error in compilation)
# cmdstan_make_local(dir = cmdstan_path(), cpp_options = "STANC2=true")

#rebuild the package
# rebuild_cmdstan(
#   dir = cmdstan_path(),
#   cores = getOption("mc.cores", 2),
#   quiet = FALSE,
#   timeout = 600
# )

#directory when using cluster
setwd("/storage/users/gmilne/test")

#directory when not using the cluster
# setwd("~/Desktop/R Projects/stan")

#read in data
# data <- read.csv("data/netherlands_95.csv")
# read in scripts
# source("R files/setparms.R")
# source("R files/demogdat.R")

#NB change of directory for using cluster
data <- read.csv("netherlands_95.csv")
source("setparms.R")
source("demogdat.R")

# data needed for parameter estimation
clean_dat <- data.frame("age_mid"=data$age_mid, "k"=data$k, "n"=data$n)
data_agrps <- length(clean_dat[,"age_mid"])

######################################################
## modify data age groups to match model age groups ##
######################################################
#create new dataset
matched_dat <- clean_dat
# x[,2] increases by one each time data age midpoint is closest match to modelled age midpoint
x <- cbind(pars$age, findInterval(pars$age, matched_dat$age_mid))
#head(x, n=12)
#returns FALSE if there's change between element i and element i+1
y1 <- diff(x[,2]) <= 0
#head(y, n=12)

# each time x[,2] increases by 1, save value of x[,1][i+1] to matched_dat$age_mid[i]
matched_ages <- vector("numeric", pars$agrps)
for(i in 1:length(y1)){
  if(y1[i]==T){
    matched_ages[i] <- NA
    
  }else if(y1[i]==F){
    matched_ages[i] <- x[,1][i+1]
  }
}
#removes last element (which is 0 because of indexing)
matched_ages <- head(matched_ages, -1)
#remove NAs
matched_ages <-matched_ages[!is.na(matched_ages)]
#save age mid points to dataset
matched_dat$age_mid <- matched_ages

###################################################
## merge data into larger df for simpler fitting ##
###################################################
mod_dat <- data.frame("age_mid"=pars$age,"k"=0, "n"=0)
merged_dat <- merge(matched_dat, mod_dat, by='age_mid', all=TRUE)
merged_dat[is.na(merged_dat)] <- 0
full_data <- data.frame("age_mid"=merged_dat[,1], "k"=merged_dat[,2], "n"=merged_dat[,3])

##check age combined as expected (differences should be the same)
# age_diff <- vector("numeric", length=length(full_data$age_mid)-1)
# for(i in 1:length(full_data$age_mid)-1){
#   if(i==1){
#     age_diff[i] <- 0
#   }else{
#     age_diff[i] <- full_data$age_mid[i+1] - full_data$age_mid[i]
#   }
# }

############################################
# read in data needed for model simulation #
############################################
agrps <- length(full_data$age_mid)
# cases <- full_data$k
# n <- full_data$n
cases <- data$k
n <- data$n
######################################
# read in other parameters for model #
######################################
t0 = 0 
ts <- seq(1,249, 1)
t <-max(ts)
N <- sum(pars$Na)
age_prop <- pars$Na/N

#index of rows of expanded df in which data exist
data_rows <- c(which(full_data$k!=0), which(full_data$k!=0)+agrps, which(full_data$k!=0)+2*agrps)

###################### 
# data list for Stan #
######################
data_si = list(
  agrps = pars$agrps, 
  data_agrps = data_agrps,
  data_rows=data_rows,
  age_prop=age_prop,
  tot_pop=N, 
  # age=pars$age,
  da=pars$da,
  d=pars$d,
  r=pars$r,
  mctr=pars$mctr,
  #mean_mctr = mean(pars$mctr),
  propfert=pars$propfert,
  K=3,  #no. state variables
  t0 = t0,
  ts = ts, 
  t=t,
  n=n, #n
  cases=cases, #k
  rel_tol = 1.0E-10, 
  abs_tol = 1.0E-10,
  max_num_steps = 1.0E3,
  inference=1, 
  doprint=0)

###################
# CmdStan running #
###################
# file <- "R files/stan_mod_simple.stan"
# file <- "R files/stan_mod_complex.stan"

#set cmdstan path (have moved cmdstan files to the FileZilla folder)
# set_cmdstan_path("/storage/users/gmilne/test/.cmdstanr/cmdstan-2.25.0")

#change of directory for cluster
# file <- "stan_mod_simple.stan"

# compile the model and save as an object
# mod <- cmdstan_model(file)

# save the compiled model as an .RData file
# save(mod, file = 'mod.RData')

# load the compiled model
load('mod.RData')

# fit the model to the data
fit <- mod$sample(
  data = data_si,
  seed = 123,
  chains = 1,
  parallel_chains = 3, 
  iter_warmup = 5,
  iter_sampling = 10,
  refresh = 1
  )

# system("say -v Karen Surfs up, bro!")

#save model output
save(fit, file = "modoutput.Rdata")

#############################
## Plotting model vs. data ##
#############################
# par(mfrow=c(2,3))
# #prevalence
# x<-c(0.140017,0.380873,0.554329,0.679229,0.769145,0.83386,0.880419,0.913903,0.937978,0.955293,0.967755,0.976734,0.983214,0.987894,0.991275,0.993714)
# plot(matched_ages, x, main="prevalence")
# lines(matched_ages, data$prev, col="red")
# #S
# x<-c(39796.5,28714.3,20794.3,15122.6,11032.2,8071.41,5913.92,4328.13,3152.86,2276.53,1622.81,1138.79,784.631,529.474,347.603,217.802)
# plot(matched_ages, x, main="S")
# #I
# x<-c(6479.4,17664.4,25864.1,32021.9,36756.3,40510.6,43541.5,45942,47681.9,48644.6,48704.4,47808.6,45957.8,43207.5,39492.9,34429.3)
# plot(matched_ages, x, main="I")
# #Im
# x<-c(5.06417e-05,3.39563e-07,2.5041e-07,-6.51204e-07,-1.44312e-06,-1.2249e-06,-9.27533e-07,-1.09858e-06,1.25964e-08,2.22018e-06,3.08956e-06,3.39268e-06,3.33301e-06,6.68076e-06,1.87276e-05,2.34077e-05)
# plot(matched_ages, x, main="Im")
# #Na
# x<-c(46275.9,46378.7,46658.4,47144.5,47788.6,48582,49455.4,50270.1,50834.8,50921.2,50327.2,48947.4,46742.4,43736.9,39840.5,34647.1)
# plot(matched_ages, x, main="Na", ylim=c(0,70000))
# lines(matched_ages, pars$Na[data_rows[1:16]], col="red")
# 
# #ggplot
# #at t=249, with births & deaths corrected in the model
# x<-c(46456,46418.2,46385.5,46357.5,46333.8,46314.2,46298.2,46285.5,46275.9,46269.1,46264.7,46262.5,46262.3,46263.8,46266.7,46270.9,46276.2,46282.5,46289.5,46297.1,46305.1,46313.6,46322.4,46331.4,46340.6,46350,46359.4,46369,46378.7,46388.5,46398.6,46408.9,46419.5,46430.5,46441.9,46453.8,46466.1,46478.9,46492.2,46506.1,46520.6,46535.7,46551.4,46567.7,46584.6,46602.1,46620.3,46639,46658.4,46678.3,46698.8,46719.8,46741.3,46763.3,46785.8,46808.8,46832.2,46856,46880.3,46905,46930.1,46955.5,46981.4,47007.7,47034.3,47061.3,47088.6,47116.4,47144.5,47173,47201.9,47231.1,47260.7,47290.8,47321.2,47352,47383.2,47414.8,47446.7,47479.1,47511.9,47545.1,47578.7,47612.7,47647.1,47681.9,47717,47752.6,47788.6,47824.9,47861.7,47898.8,47936.3,47974.2,48012.4,48051,48090,48129.3,48168.9,48208.9,48249.2,48289.8,48330.8,48372,48413.5,48455.2,48497.2,48539.5,48582,48624.7,48667.6,48710.7,48754,48797.4,48841,48884.6,48928.4,48972.3,49016.2,49060.2,49104.2,49148.3,49192.3,49236.3,49280.3,49324.2,49368,49411.8,49455.4,49498.9,49542.2,49585.4,49628.3,49671.1,49713.6,49755.8,49797.8,49839.4,49880.7,49921.7,49962.3,50002.4,50042.2,50081.5,50120.3,50158.6,50196.3,50233.5,50270.1,50306.1,50341.5,50376.2,50410.2,50443.4,50475.9,50507.6,50538.5,50568.5,50597.7,50626,50653.3,50679.7,50705.1,50729.4,50752.7,50774.9,50796,50816,50834.8,50852.4,50868.8,50883.9,50897.7,50910.2,50921.3,50931.1,50939.5,50946.4,50951.9,50956,50958.5,50959.4,50958.8,50956.7,50952.9,50947.5,50940.4,50931.6,50921.2,50909,50895.1,50879.4,50861.9,50842.6,50821.5,50798.6,50773.8,50747.2,50718.6,50688.2,50655.9,50621.6,50585.4,50547.3,50507.2,50465.1,50421.1,50375.1,50327.2,50277.2,50225.3,50171.4,50115.5,50057.6,49997.8,49935.9,49872,49806.1,49738.2,49668.3,49596.4,49522.4,49446.4,49368.4,49288.3,49206.1,49121.9,49035.7,48947.4,48857,48764.5,48669.9,48573.2,48474.5,48373.6,48270.6,48165.5,48058.3,47949.1,47837.7,47724.2,47608.7,47491.1,47371.4,47249.7,47125.9,47000.1,46872.3,46742.4,46610.7,46476.9,46341.2,46203.5,46064,45922.5,45779.1,45633.7,45486.5,45337.3,45186.3,45033.2,44878.2,44721.3,44562.3,44401.4,44238.4,44073.3,43906.2,43736.9,43565.5,43391.9,43216.1,43038,42857.6,42674.8,42489.6,42302,42111.9,41919.2,41723.9,41526,41325.4,41122,40915.8,40706.7,40494.7,40279.7,40061.6,39840.5,39616.1,39388.6,39157.7,38923.5,38685.9,38444.8,38200,37951.6,37699.4,37443.3,37183.3,36919.2,36650.9,36378.4,36101.4,35820,35533.9,35243.2,34947.6,34647.1,34341.6,34031,33715.2,33394,33067.5,32735.5,32397.9,32054.7,31705.9,31351.3,30990.9,30624.8,30252.8,29875,29491.4,29102,28706.7,28305.7,27899.1,27486.8,27068.9,26645.6,26216.9,25783.1,25344.1,24900.1,24451.3,23997.8,23539.8,23077.4,22610.9,22140.4,21666.1,21188.4,20707.3,20223.3,19736.6,19247.5,18756.3,18263.4,17769.2,17273.9,16778.2,16282.3,15786.7,15292.1,14798.8,14307.6,13818.8,13333.1,12851.1,12373.4,11900.5,11433.1,10971.7,10516.9,10069.3,9629.34,9197.6,8774.55,8360.66,7956.35,7562.02,7178.05,6804.74,6442.4,6091.27,5751.55,5423.43,5107.02,4802.42,4509.68,4228.8,3959.78,3702.54,3456.99,3223.01,3000.43,2789.07,2588.71,2399.12,2220.03,2051.16,1892.2,1742.86,1602.79,1471.66,1349.11,1234.8,1128.36,1029.43)
# length(x)
# df <- data.frame(pars$age, x, pars$Na)
# 
# ### or using data restricted to age groups found in the data
# x<-c(46275.9,46378.7,46658.4,47144.5,47788.6,48582,49455.4,50270.1,50834.8,50921.2,50327.2,48947.4,46742.4,43736.9,39840.5,34647.1)
# df <- data.frame(matched_ages, x, pars$Na[data_rows[1:16]])
# 
# ggplot(data=df, aes(x=matched_ages, y=x)) + 
#   geom_line() + 
#   geom_point(data=df, mapping=aes(x=matched_ages, y=pars$Na[data_rows[1:16]]), size=.2, col = "grey") + 
#   scale_x_continuous(breaks=c(0,max(matched_ages), seq(0, max(matched_ages), 20)), limits=c(0,max(matched_ages)), expand=c(0,0)) + 
#   scale_y_continuous(limits=c(0,70000), expand=c(0,0)) +
#   ylab("Population size") + 
#   xlab("Age (years)") + 
#   theme_light(base_family="Times", base_size = 12, base_line_size = 0) + 
#   theme(plot.margin=unit(c(1,1,1,1),"cm")) +
#   theme(panel.grid.major = element_blank()) +
#   theme(panel.grid.minor = element_blank())
# #save plot
# #ggsave("plots/popsize.eps", units="in", heigh=6, width=6, family = "Times")
# 
# #reasonable lambda0 prior
# #lognormal(log(.05), .1)
# par(mfrow=c(1,1))
# hist(exp(rnorm(1000, log(.06), .3)), xlim=c(0,0.1))
# 
# ##############
# #save to csv#
# ##############
# # 26/11/2020 #
# ##############
# ## 1 chain, 10 sampling iterations,
# #1 parameter (lambda0), 1 yr/age category, 80 age categories ##
# # fit$save_output_files(dir = ".", basename = NULL, timestamp = TRUE, random = TRUE)
# fit_mod <- read.csv("data/stanfit-261120.csv")   
# 

##view contents of .Rout file

