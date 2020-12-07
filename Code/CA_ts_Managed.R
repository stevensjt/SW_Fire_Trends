###########################################################################################################
# 
# This script examines high-severity class area as a function of time for managed fires
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
funchir::stale_package_check('CA_ts_Managed.R')

# disable scientific notation
options(scipen=999)
options(digits=6)

# attach data
ts_patch_metrics <- read.csv(file="GEE_1984_2017_final_patch_metrics_709.csv", header=TRUE, sep=",")
names(ts_patch_metrics)

###########################################################################################################
# MEAN CA by YEAR - MANAGED
###########################################################################################################

# filter for managed fires only
filter.df <- filter(ts_patch_metrics, FIRE_TYPE2 == "OTHER")

## aggregate data by CA by year then calculate mean 
mean_ca.df <- aggregate(filter.df[, 6], list(filter.df[,2]), mean)
## rename columns
colnames(mean_ca.df)[colnames(mean_ca.df)=="Group.1"] <- "YEAR"
colnames(mean_ca.df)[colnames(mean_ca.df)=="x"] <- "CA"
## convert character data to numeric 
mean_ca.df <- mutate_all(mean_ca.df, function(x) as.numeric(as.character(x)))
## format to number of digits
mean_ca.df <- mean_ca.df %>% 
  mutate_if(is.numeric, round, digits = 5)
mean_ca.df

# add number of fires per year
for(r in 1:nrow(mean_ca.df)){
  y <- mean_ca.df$YEAR[r]
  mean_ca.df$n[r] <- 
    length(which(filter.df$YEAR==y))
}

# save
write.csv(mean_ca.df, "./mean_ca_man.csv")

# remove years with only 1 fire
mean_ca.df <- mean_ca.df[-c(1),]

# quick plot of time series 
myts <- ts((mean_ca.df$CA), start=c(1995:2001,2003:2014,2016:2017))
myts.log <- ts(log(mean_ca.df$CA), start=c(1995:2001,2003:2014,2016:2017))
plot(myts)
list(myts)
plot(myts.log)
list(myts.log)
mk.test(myts, alternative = c("two.sided"),
        continuity = TRUE)

# test for lag-1 autocorrelation using Durbin Watson test 
x <- mean_ca.df$YEAR
y <- mean_ca.df$CA
ylog <- log(y)
mylm <- lm(y~x)
myloglm <- lm(ylog~x)
set.seed(123)
durbinWatsonTest(mylm)

# examine ACF and PACF plots for autocorrelation as well
acf(myts)
pacf(myts)
acf(myts.log)
pacf(myts.log)
acf(diff(myts.log))
pacf(diff(myts.log))
auto.arima(myts.log)
auto.arima(diff(myts.log))

# if no autocorrelation - try Mann Kendall trend test
mk.test(myts, alternative = c("two.sided"),
        continuity = TRUE)

###########################################################################################################
# If No Autocorrelation - Test OLS Assumptions
###########################################################################################################

# set x and y 
x <- mean_ca.df$YEAR
y <- mean_ca.df$CA
ylog <- log(y)
fit <- lm(ylog ~ x) #***
summary(fit)
plot(x, ylog)
# normality test
shapiro.test(resid(fit))
# variance test 
ncvTest(fit)
# all diagnostics
shapiro.test(resid(fit))
gvlma(fit)
autoplot(fit)
rmse(fit)
bbmle::AICctab(fit1,fit,base=TRUE,logLik=TRUE, weights=TRUE)

###########################################################################################################
# Plot Best Model
###########################################################################################################

# using linear model
fit <- lm(ylog ~ x) #***
yfit <- fitted(fit)

cols <- c("SUP" = "indianred3")

ca_ts <- ggplot(data = mean_ca.df, aes(x = YEAR, y = ylog)) + 
  geom_point(color= 'turquoise3', size = 3) + 
  geom_line(color = "turquoise3", size = 0.5) +
  geom_line(aes(y = yfit), color = "turquoise3", linetype = 1, size = 1) +
  #geom_smooth(method = "lm", se = TRUE, color = 'turquoise3', size = 1.5) + 
  theme_bw() + 
  ggtitle("\nHigh-Severity Area in Managed Fires") + 
  labs(x="\nYEAR", y = "\nln(ha)") +
  scale_y_continuous(breaks=seq(2, 7, 1)) +
  scale_x_continuous(breaks=seq(1999,2015,5)) +
  panel_border() + 
  scale_color_manual(name = "class",
                     labels = c("SUP"),
                     values = cols)

ca_ts <- ca_ts + theme(plot.title = element_text(family = "Calibri", color = "black", size =16, hjust = 0.5, lineheight = 0.8, face="bold"),
                       axis.title.x = element_text(family = "Calibri", color="black", size=12, face = "bold"),
                       axis.text.x = element_text(family = "Calibri", color="black", size=12),
                       axis.title.y = element_text(family = "Calibri", color="black", size=12, face = "bold"),
                       axis.text.y = element_text(family = "Calibri", color="black", size=12),
                       #plot.margin = unit(c(1, 1.5, 1, 1.5), "cm"), #top, left, bottom, right
                       panel.border = element_rect(colour = "black", fill=NA, size=0.5))
ca_ts

OutPath <- "D:/Patch_Metrics/R/Time_Series/CA_Time_Series/PDFs/ts_CA_managed.pdf"
ggsave(OutPath, plot=ca_ts, device = cairo_pdf)

#########################################################################################################
#########################################################################################################
