############################################################################
## Linear models & plotting of relationship between seroprevalence & time ##
############################################################################

## Fit a hierarchical mixed effects model
library(lme4)
library(ggplot2)
library(RColorBrewer)

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

# calculate prevalence
df$prev <- df$k/df$n

# standardise year
df$yrstd <- (df$year-mean(df$year))/sd(df$year)  #intercept where year is avg year rather than =0

## complex model (country as random intercept + random gradient)
m0 <- lmer(prev~yrstd+(yrstd|country),weights=n/(I(max(n))), data=df) #Martin's LME

## simple model (country as random intercept but not random gradient)
m1 <- lmer(prev~yrstd+(1|country),weights=n/(I(max(n))), data=df)

## compare complex vs. simple model
anova(m0, m1) #no evidence to support the more complex model

## estimate accuracy of publication year as proxy for sampling year
mean(df$year_published) - mean(df$year)
t.test(df$year_published, df$year, paired=T)

## estimate annual declines in seroprevalence (model with un-standardised year)
m1 <-  lmer(prev~year+(1|country),weights=n/(I(max(n))), data=df)

# Fitted values
df$yfit1 <- fitted(m1)

## Plot fit of simpler linear model to the data ##
# rename Iran to have shorter name for plotting
df$country <- factor(df$country, labels= c(levels(df$country)[1:6], "Iran", levels(df$country)[8:11]))

# create custom color scale
myColors        <- brewer.pal(11,"PuOr")
myColors[6]     <- brewer.pal(11,"PRGn")[5]   #make pale colour bolder
names(myColors) <- levels(df$country)
colScale        <- scale_colour_manual(name = "plotcols",values = myColors)

# make seroprevalence as %
df$prev  <- df$prev*100
df$ci_lo <- df$ci_lo*100
df$ci_up <- df$ci_up*100
df$yfit1 <- df$yfit1*100

# plot
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

#Save
#PDF
ggsave(filename = "plots/lm-global-seroprev.pdf",
       device = cairo_pdf, height = 6, width = 6, units = "in")

#PNG
ggsave(filename = "plots/lm-global-seroprev.png",
       dpi=600, height = 6, width = 6, units = "in")
