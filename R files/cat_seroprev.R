## Exploratory analysis for Discussion: Seroprevalence in cats
library(ggplot2)

#Source: Toxoplasma gondii infection in domestic and wild felids as public health concerns: a systematic review and meta-analysis
dat <- read.csv("data/cat_seroprev.csv")

dat$seroprev <- dat$k/dat$n

dat <- subset(dat, !is.na(year1)) #omit NAs

# Calculate sampling year
for(i in 1:nrow(dat)){
  if ( !is.na(dat$year2[i] )) { #if data sampled in a range of years
    dat$year[i] <- mean(c(dat$year1[i], dat$year2[i])) #calculate simple mean year
    
  } else if (is.na(dat$year2)[i]) { #if data sampled in one year
    dat$year[i] <- dat$year1[i]
  }
  
}

plot(dat$year, dat$seroprev)

# omit outlier (year 1880) to see any trends more clearly
which(dat$year < 1900)
dat <- dat[-113,]

cat_seroprev <- plot(dat$year, dat$seroprev, ylab = "Seroprevalence", xlab = "Year")
mod <- lm(dat$seroprev ~ dat$year)
confint(mod) # CIs overlapping 0?
abline(a=mod$coefficients[1], b=mod$coefficients[2], col="red")


ggplot(data=dat, aes(x=year, y=seroprev)) + 
  geom_point() + 
  geom_line(aes(x=year, y=mod$fitted.values), col="red") + 
  xlab("Year")+
  ylab("Seroprevalence") + 
  theme_light(base_size = 12, base_line_size = 0, base_family = "Times") + 
  theme(legend.position = "none", panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), plot.margin=unit(c(rep(0.4,4)),"cm"))

#save plot
ggsave(filename = "plots/cat_seroprev.png", dpi=600,
       height = 6, width = 6, units = "in")
