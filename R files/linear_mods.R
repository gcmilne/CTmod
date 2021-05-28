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
df$yrstd <- (df$year-mean(df$year))/sd(df$year)

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
m1<-lmer(prev~yrstd+(1|country), data=df)
df$yfit1 <- fitted(m1)

#CIs around estimate
confint(m1)

#create custom color scale
library(RColorBrewer)
myColors <- brewer.pal(11,"PuOr")
myColors[6] <- brewer.pal(11,"PRGn")[5]   #make pale colours bolder
names(myColors) <- levels(df$country)
colScale <- scale_colour_manual(name = "plotcols",values = myColors)

#plot
p1 <- ggplot(data=df, aes(y=yfit1, x=year, colour = country)) +
  geom_line() +
  geom_point(data=df, aes(y=prev, x=year)) +
  facet_wrap(~country) +
  ylab("Seroprevalence") +
  xlab("Year") +
  theme_light(base_size = 12, base_line_size = 0, base_family = "Times") +
  theme(strip.background = element_rect(fill=myColors[9])) +
  # theme(strip.text = element_text(colour = 'black')) +
  theme(legend.position = "none", panel.grid = element_blank(), axis.ticks = element_blank())

#add colours
p1 + colScale

# save plot
# ggsave(filename = "plots/lm-global-seroprev.png", width = 6, height = 6, 
       # units = "in", family = "Times")

## no evidence to support the more complex model
anova(m0, m1)
