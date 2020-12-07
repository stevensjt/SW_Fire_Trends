###########################################################################################################
# 
# This script examines high-severity percentage of like adjacencies as a function of time for managed fires
# Author: Megan Singleton 
# Updated: 12/04/2020
#
###########################################################################################################

# load libraries
library(forecast)
library(tseries)
library(dplyr) #formating
library(bbmle) #AICctab 
library(car) #durbinwatson test
library(ggplot2) #graphics
library(cowplot) #graphics
library(extrafont) #fonts
library(ggfortify) #autoplot
library(trend) #mktest
library(sjstats) #rmse
library(gvlma) #model diagnostics
library(nlme) #gls models
library(lmtest) #lr test
font_import() #register fonts, must do once per session
loadfonts()
loadfonts(device = "win")

# set working directory 
setwd("D:/Patch_Metrics/Patch_Metrics_Complete_Analysis/R/Time_Series/")

# check to see if there are any unused packages 
funchir::stale_package_check('PLADJ_ts_Managed.R')

# disable scientific notation
options(scipen=999)
options(digits=6)

# attach data
ts_patch_metrics <- read.csv(file="GEE_1984_2017_final_patch_metrics_709.csv", header=TRUE, sep=",")
names(ts_patch_metrics)

###########################################################################################################
# MEAN PLADJ by YEAR - MANAGED
###########################################################################################################

# filter dataset by managed fires only
filter.df <- filter(ts_patch_metrics, FIRE_TYPE2 == "OTHER")

##aggregate data by PLADJ by year then calculate mean 
mean_pladj.df <- aggregate(filter.df[, 44], list(filter.df[,2]), mean)
##rename columns
colnames(mean_pladj.df)[colnames(mean_pladj.df)=="Group.1"] <- "YEAR"
colnames(mean_pladj.df)[colnames(mean_pladj.df)=="x"] <- "PLADJ"
##convert character data to numeric 
mean_pladj.df <- mutate_all(mean_pladj.df, function(x) as.numeric(as.character(x)))
##format to number of digits
mean_pladj.df <- mean_pladj.df %>% 
  mutate_if(is.numeric, round, digits = 5)
mean_pladj.df

#number of fires by year
for(r in 1:nrow(mean_pladj.df)){
  y <- mean_pladj.df$YEAR[r]
  mean_pladj.df$n[r] <- 
    length(which(filter.df$YEAR==y))
}

# save
write.csv(mean_pladj.df, "./mean_pladj_man.csv")

# remove years with only 1 fire
mean_pladj.df <- mean_pladj.df[-c(1),]

# time series analysis 
myts <- ts((mean_pladj.df$PLADJ), start=c(1995:2001, 2003:2014, 2016:2017))
myts.assqt <- ts(asin(sqrt(mean_pladj.df$PLADJ/100)), start=c(1995:2001, 2003:2014, 2016:2017))
plot(myts)
list(myts)
plot(myts.assqt)
list(myts.assqt)
mk.test(myts, alternative = c("two.sided"),
        continuity = TRUE)

# test if 1st lag is autocorrelated using Durbin-Watson test 
x <- mean_pladj.df$YEAR
y <- mean_pladj.df$PLADJ
yarc <- asin(sqrt(y/100))
mylm <- lm(y~x)
myarclm <- lm(yarc~x)
set.seed(123)
durbinWatsonTest(mylm)

#examine ACF and PACF plots for autocorrelation as well
acf(myts)
pacf(myts)
acf(myts.assqt)
pacf(myts.assqt)
auto.arima(myts)
auto.arima(myts.assqt)

###########################################################################################################
# Model Selection - Beta Regression and OLS
###########################################################################################################

# set x and y 
x <- mean_pladj.df$YEAR
y <- mean_pladj.df$PLADJ
yarc <- asin(sqrt(y/100))
yper <- y/100

# create and fit models
fit1 <- lm(y ~ x)
fit2 <- lm(yarc ~ x)
fit3 <- betareg(yper ~ x, data = mean_pladj.df) #****

# mean center to mitigate collinearity in quadratic term
time.lm <- mean_pladj.df$YEAR
time <- seq(from = 1, to = 21, by = 1)
time <- time - mean(time)
time2 <- time^2
fit4 <- lm(y ~ time + time2)
fit5 <- lm(yarc ~ time + time2)
fit6 <- betareg(yper ~ time + time2, data = mean_pladj.df) #logit link

# compare competing models 
bbmle::AICctab(fit3,fit6,base=TRUE,logLik=TRUE, weights=TRUE)

# plot model diagnostics 
autoplot(fit1)
autoplot(fit2)
plot(fit3) #****
autoplot(fit4)
autoplot(fit5)
plot(fit6)

# summary statistics
summary(fit1) 
summary(fit2)
summary(fit3) #****
summary(fit4)
summary(fit5)
summary(fit6)

# normality test
shapiro.test(resid(fit1))
shapiro.test(resid(fit2))
shapiro.test(resid(fit3)) #****
shapiro.test(resid(fit4))
shapiro.test(resid(fit5))
shapiro.test(resid(fit6))

# all OLS diagnostics
gvlma(fit1)
gvlma(fit2)
gvlma(fit4)
gvlma(fit5)

#RMSE/AIC/SUMMARY
rmse(fit3)
AIC(fit3, k=2)
lrtest(fit3)

###########################################################################################################
# Plot Untransformed Model 
###########################################################################################################
x <- mean_pladj.df$YEAR
y <- mean_pladj.df$PLADJ
fit3 <- betareg(yper ~ x, data = mean_pladj.df) #****
yfit <- fitted(fit3)
yfit <- yfit*100

cols <- c("MAN" = "turquoise3")

pladj_ts <- ggplot(data = mean_pladj.df, aes(x = x, y = y)) + 
  geom_point(data = mean_pladj.df, aes(x = x, y = y, colour = "MAN"), size = 2) +
  geom_line(color = "turquoise3", size = 0.7) +
  geom_line(aes(y = yfit), color = "turquoise3", linetype = 1, size = 1) +
  theme_bw() + 
  ggtitle("\nPercent Like Adjacencies (PLADJ) \nin all Managed Fires") + 
  labs(x="\nYEAR", y = "\nPercent") +
  scale_y_continuous(breaks=seq(50, 80, 10)) +
  #scale_y_continuous(breaks=seq(0.80, 1.20, 0.05)) +
  scale_x_continuous(breaks=seq(1984,2017,5)) +
  panel_border() + 
  scale_color_manual(name = "class",
                     labels = c("SUP"),
                     values = cols)

pladj_ts <- pladj_ts + theme(plot.title = element_text(family = "Calibri", color = "black", size =16, hjust = 0.5, lineheight = 0.8, face="bold"),
                             axis.title.x = element_text(family = "Calibri", color="black", size=12, face = "bold"),
                             axis.text.x = element_text(family = "Calibri", color="black", size=12),
                             axis.title.y = element_text(family = "Calibri", color="black", size=12, face = "bold"),
                             axis.text.y = element_text(family = "Calibri", color="black", size=12),
                             #plot.margin = unit(c(1, 1.5, 1, 1.5), "cm"), #top, left, bottom, right
                             panel.border = element_rect(colour = "black", fill=NA, size=0.5))
pladj_ts

# save
OutPath <- "./PLADJ_ts_Managed.pdf"
ggsave(OutPath, plot=pladj_ts, device = cairo_pdf)

