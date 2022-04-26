################################
## Linear mixed effect models ##
################################

###################
## Load packages ##
###################
library(lme4)
library(ggplot2)
library(RColorBrewer)
library(binom)

##################
## Load scripts ##
##################
source("R files/diagnostics.R")
source("R files/funcs.R")

#############
# Set fonts #
#############
if(.Platform$OS.type == "windows") { # set Times New Roman font on Windows
  # library(extrafont)
  # font_import()
  # loadfonts(device = "win")
  windowsFonts(Times=windowsFont("TT Times New Roman")) 
}

###############
## Load data ##
###############
df <- readRDS("data/global_data.rds")

###################################################
## Adjust observed prevalence to true prevalence ##
###################################################
# calculate observed prevalence
df$prev <- df$k/df$n

# adjust prevalence by se & sp
se_sp <- find_diagnostic_values(df, assays)  #find se/sp values
df$tp <- (df$prev + (se_sp$sp - 1)) / (se_sp$sp + (se_sp$se - 1))  #adjust to true prevalence

# calculate true value of k (no. positive individuals)
k_seq <- seq(0, 15000) #vector to take guess at k

for (i in 1:nrow(df)){
  k_index <- which.min(abs(k_seq/df$n[i] - df$tp[i]))
  df$tk[i] <- k_seq[k_index]  #"true" k
}

# plot(df$tk/df$n, df$tp)  #check concordance

# calculate 95% CIs around true prevalence
cis <- binom.confint(x=df$tk, n=df$n, conf.level=0.95, methods="exact")

# Store lower and upper estimates
df$tp_ci_lo <- cis$lower
df$tp_ci_up <- cis$upper


#################################
## Linear mixed-effects models ##
#################################

# Using true prevalence
# standardise year
df$yrstd <- (df$year-mean(df$year))/sd(df$year)  #intercept where year is avg year rather than =0

## complex model (country as random intercept + random gradient)
m0 <- lmer(tp~yrstd+(yrstd|country),weights=n/(I(max(n))), data=df)

## simple model (country as random intercept but not random gradient)
m1 <- lmer(tp~yrstd+(1|country),weights=n/(I(max(n))), data=df)

## compare complex vs. simple model
anova(m0, m1) #no evidence to support the more complex model

## estimate annual declines in seroprevalence (simple model with un-standardised year)
m1 <-  lmer(tp~year+(1|country),weights=n/(I(max(n))), data=df)

## Calculate 95% CIs for annual seroprevalence decline (%)
confint(m1, method="Wald", parm = "year")*100 

# Fitted values
df$yfit1 <- fitted(m1)

## Plot fit ##
# rename Iran to have shorter name for plotting
df$country <- factor(df$country, labels= c(levels(df$country)[1:6], "Iran", levels(df$country)[8:11]))

# create custom color scale
myColors        <- brewer.pal(11,"PuOr")
myColors[6]     <- brewer.pal(11,"PRGn")[5]   #make pale colour bolder
names(myColors) <- levels(df$country)
colScale        <- scale_colour_manual(name = "plotcols",values = myColors)

# make seroprevalence as %
df$prev  <- df$tp*100
df$ci_lo <- df$tp_ci_lo*100
df$ci_up <- df$tp_ci_up*100
df$yfit1 <- df$yfit1*100  #simple mod
df$yfit0 <- df$yfit0*100  #complex mod

# plot simple model
p1 <- ggplot(data=df, aes(y=yfit1, x=year, colour = country)) +
  geom_line() +
  geom_point(data=df, aes(y=prev, x=year), size=.8) +
  geom_errorbar(aes(width=.3, x=year, ymin = ci_lo, ymax = ci_up), size=0.2) +
  facet_wrap(~country, ncol=3) +
  scale_x_continuous(limits=c(1980, 2020), breaks=seq(1983, 2016, 11)) + 
  ylab("Seroprevalence (%)") +
  xlab("Median sampling year") +
  theme_light(base_size = 12, base_line_size = 0, base_family = "Times") +
  theme(strip.background = element_rect(fill="white")) +
  theme(strip.text = element_text(colour = 'black')) +
  theme(legend.position = "none", panel.grid = element_blank(), axis.ticks = element_blank())

#add colours
p1 + colScale

# #PNG
ggsave(filename = "plots/lm-adjusted-seroprev.png",
       dpi=600, height = 6, width = 6, units = "in")
