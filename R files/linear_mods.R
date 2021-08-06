############
## Linear models
############

## Fit a hierarchical mixed effects model for all countries together
library(lme4)
library(ggplot2)

# read in data
df <- readRDS("data/global_data.rds")

# calculate prevalence
df$prev <- df$k/df$n

# standardise year
df$yrstd <- (df$year-mean(df$year))/sd(df$year)  #intercept where year is avg year rather than =0

## over-fitted model
m0 <- lmer(prev~yrstd+(yrstd|country), data=df)

df$yfit0 <- fitted(m0)

p0 <- ggplot(data=df, aes(y=yfit0, x=year, group=country, col=country)) +
  geom_line(show.legend = F) +
  geom_point(data=df, aes(y=prev, x=year, col=country), show.legend = F) +
  facet_wrap(~country) + 
  ylab("Seroprevalence") + 
  xlab("Year") + 
  theme_light()

## simpler model
# m1<-lmer(prev~yrstd+(1|country), data=df)

## estimate annual declines in seroprevalence
# model with un-standardised year
m1<-lmer(prev~year+(1|country), data=df)
df$yfit1 <- fitted(m1)

#CIs around estimate
confint(m1)


# rename Iran to have shorter name for plotting
df$country <- factor(df$country, labels= c(levels(df$country)[1:6], "Iran", levels(df$country)[8:11]))

#create custom color scale
library(RColorBrewer)
# myColors <- brewer.pal(11,"PuOr")
# myColors[6] <- brewer.pal(11,"PRGn")[5]   #make pale colours bolder

display.brewer.all(n=11,type="div",exact.n=TRUE)
myColors <- brewer.pal(11,"Spectral")


names(myColors) <- levels(df$country)
colScale <- scale_colour_manual(name = "plotcols",values = myColors)


#plot
p1 <- ggplot(data=df, aes(y=yfit1, x=year, colour = country)) +
  geom_line() +
  geom_point(data=df, aes(y=prev, x=year)) +
  facet_wrap(~country, ncol=3) +
  scale_x_continuous(limits=c(1980, 2020), breaks=seq(1983, 2016, 11)) + 
  ylab("Seroprevalence") +
  xlab("Median sampling year") +
  theme_light(base_size = 12, base_line_size = 0, base_family = "Times") +
  theme(strip.background = element_rect(fill="white")) +
  theme(strip.text = element_text(colour = 'black')) +
  theme(legend.position = "none", panel.grid = element_blank(), axis.ticks = element_blank()#, 
        # strip.background = element_blank()
        )

#add colours
p1 + colScale

#Save
ggsave(filename = "plots/lm-global-seroprev.pdf",
       device = cairo_pdf, height = 6, width = 6, units = "in")

## no evidence to support the more complex model
anova(m0, m1)


## estimate accuracy of publication year as proxy for sampling year

mean(df$year_published) - mean(df$year)
t.test(df$year_published, df$year, paired=T)
