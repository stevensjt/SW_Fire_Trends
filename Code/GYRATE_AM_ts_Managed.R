###########################################################################################################
# 
# This script examines high-severity patch radius of gyration as a function of time for managed fires
# Author: Megan Singleton 
# Updated: 12/04/2020
#
###########################################################################################################

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
funchir::stale_package_check('GYRATE_AM_ts_Managed.R')

# disable scientific notation
options(scipen=999)
options(digits=6)

# attach data
ts_patch_metrics <- read.csv(file="GEE_1984_2017_final_patch_metrics_709.csv", header=TRUE, sep=",")
names(ts_patch_metrics)

###########################################################################################################
# MEAN GYRATE_AM by YEAR - MANAGED
###########################################################################################################

# filter dataset by managed fires only
filter.df <- filter(ts_patch_metrics, FIRE_TYPE2 == "OTHER")

##aggregate data by GYRATE_AM by year then calculate mean 
mean_gyrate.df <- aggregate(filter.df[, 16], list(filter.df[,2]), mean)
##rename columns
colnames(mean_gyrate.df)[colnames(mean_gyrate.df)=="Group.1"] <- "YEAR"
colnames(mean_gyrate.df)[colnames(mean_gyrate.df)=="x"] <- "GYRATE_AM"
##convert character data to numeric 
mean_gyrate.df <- mutate_all(mean_gyrate.df, function(x) as.numeric(as.character(x)))
##format to number of digits
mean_gyrate.df <- mean_gyrate.df %>% 
  mutate_if(is.numeric, round, digits = 5)
mean_gyrate.df

# add number of fires per year
for(r in 1:nrow(mean_gyrate.df)){
  y <- mean_gyrate.df$YEAR[r]
  mean_gyrate.df$n[r] <- 
    length(which(filter.df$YEAR==y))
}

# save data 
write.csv(mean_gyrate.df, "./mean_gyrate_man.csv")

# remove years with only 1 fire
mean_gyrate.df <- mean_gyrate.df[-c(1),]

# time series analysis 
myts <- ts((mean_gyrate.df$GYRATE_AM), start=c(1995:2001, 2003:2014, 2016:2017))
myts.log <- ts(log(mean_gyrate.df$GYRATE_AM), start=c(1995:2001, 2003:2014, 2016:2017))
mk.test(myts, alternative = c("two.sided"),
        continuity = TRUE)

# test if autocorrelated using Durbin Watson test 
x <- mean_gyrate.df$YEAR
y <- mean_gyrate.df$GYRATE_AM
mylm <- lm(y~x)
set.seed(123)
durbinWatsonTest(mylm)

# examine ACF and PACF plots for autocorrelation as well
acf(myts)
pacf(myts)
acf(myts.log)
pacf(myts.log)
acf(diff(myts.log))
auto.arima(myts)
auto.arima(myts.log)

###########################################################################################################
# If No Autocorrelation - Test OLS Assumptions
###########################################################################################################

# set x and y
x <- mean_gyrate.df$YEAR
y <- mean_gyrate.df$GYRATE_AM
ylog <- log(y)
ysqrt <- sqrt(y)

# mean center to mitigate collinearity
time.lm <- mean_gyrate.df$YEAR
time <- seq(from = 1, to = 21, by = 1)
time <- time - mean(time)
time2 <- time^2

# fit models
fit1 <- lm(y ~ x)
fit2 <- lm(ylog ~ x)
fit3 <- lm(ysqrt ~ x)
fit4 <- lm(y ~ time + time2)
fit5 <- lm(ylog ~ time + time2) #***
fit6 <- lm(ysqrt ~ time + time2)

# diagnostic plots
autoplot(fit1)
autoplot(fit2)
autoplot(fit3)
autoplot(fit4)
autoplot(fit5) #*
autoplot(fit6)

# normality test
shapiro.test(resid(fit1))
shapiro.test(resid(fit2))
shapiro.test(resid(fit3))
shapiro.test(resid(fit4))
shapiro.test(resid(fit5)) #*
shapiro.test(resid(fit6))

# all diagnositics
gvlma(fit1)
gvlma(fit2)
gvlma(fit3)
gvlma(fit4)
gvlma(fit5) #*
gvlma(fit6)

# summary
summary(fit1)
summary(fit2)
summary(fit3)
summary(fit4)
summary(fit5) #***
summary(fit6)

# RMSE/AIC/SUMMARY
lrtest(fit5)
summary(fit5)
rmse(fit5)
AIC(fit5, k=2)
bbmle::AICctab(fit2,fit5,base=TRUE,logLik=TRUE, weights=TRUE)

# histogram
ggplot(mean_gyrate.df, aes(x=ylog)) + 
  geom_histogram(color="black", fill="white", bins = 15)

###########################################################################################################
# Plot Best Model
###########################################################################################################

fit5 <- lm(ylog ~ time + time2) #***
yfit <- fitted(fit5)
cols <- c("MAN" = "turquoise3")

gyrate_ts <- ggplot(data = mean_gyrate.df, aes(x = YEAR, y = ylog)) + 
  geom_point(data = mean_gyrate.df, aes(x = YEAR, y = ylog, colour = "MAN"), size = 2) +
  geom_line(color = "turquoise3", size = 0.7) +
  geom_line(aes(y = yfit), color = "turquoise3", linetype = 1, size = 1) +
  theme_bw() + 
  ggtitle("\nArea Weighted Mean Radius of\nGyration Distribution in All Managed Fires") + 
  labs(x="\nYEAR", y = "\nln(meters)") +
  scale_x_continuous(breaks=seq(1995,2015,5)) +
  panel_border() + 
  scale_color_manual(name = "class",
                     labels = c("SUP"),
                     values = cols)

gyrate_ts <- gyrate_ts + theme(plot.title = element_text(family = "Calibri", color = "black", size =16, hjust = 0.5, lineheight = 0.8, face="bold"),
                               axis.title.x = element_text(family = "Calibri", color="black", size=12, face = "bold"),
                               axis.text.x = element_text(family = "Calibri", color="black", size=12),
                               axis.title.y = element_text(family = "Calibri", color="black", size=12, face = "bold"),
                               axis.text.y = element_text(family = "Calibri", color="black", size=12),
                               #plot.margin = unit(c(1, 1.5, 1, 1.5), "cm"), #top, left, bottom, right
                               panel.border = element_rect(colour = "black", fill=NA, size=0.5))
gyrate_ts

OutPath <- "D:/Patch_Metrics/GEE_Severity_Dataset/R/Time_Series/GYRATE_AM_Time_Series/PDFs/GYRATE_AM_ts_managed_quad.pdf"
ggsave(OutPath, plot=gyrate_ts, device = cairo_pdf)

#########################################################################################################


