###########################################################################################################
# 
# This script examines the high-severity stand-decay coefficient as a function of time for managed fires
# Author: Megan Singleton 
# Updated: 12/04/2020
#
###########################################################################################################

library(forecast)
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
font_import() #register fonts, must do once per session
loadfonts()
loadfonts(device = "win")

# set working directory 
setwd("D:/Patch_Metrics/Patch_Metrics_Complete_Analysis/R/Time_Series/")

# check to see if there are any unused packages 
funchir::stale_package_check('SDC_ts_Managed.R')

# disable scientific notation
options(scipen=999)
options(digits=6)

# attach data
ts_patch_metrics <- read.csv(file="GEE_1984_2017_final_patch_metrics_709.csv", header=TRUE, sep=",")
names(ts_patch_metrics)

###########################################################################################################
# MEAN SDC by YEAR - MANAGED
###########################################################################################################

# filter dataset by managed fire onlu
filter.df <- filter(ts_patch_metrics, FIRE_TYPE2 == "OTHER", SDC > 0)

## aggregate data by SDC by year then calculate mean 
mean_sdc.df <- aggregate(filter.df[, 5], list(filter.df[,2]), mean)
## rename columns
colnames(mean_sdc.df)[colnames(mean_sdc.df)=="Group.1"] <- "YEAR"
colnames(mean_sdc.df)[colnames(mean_sdc.df)=="x"] <- "SDC"
## convert character data to numeric 
mean_sdc.df <- mutate_all(mean_sdc.df, function(x) as.numeric(as.character(x)))
## format number of digits
mean_sdc.df <- mean_sdc.df %>% 
  mutate_if(is.numeric, round, digits = 5)
mean_sdc.df

# calculate number of fires per year
for(r in 1:nrow(mean_sdc.df)){
  y <- mean_sdc.df$YEAR[r]
  mean_sdc.df$n[r] <- 
    length(which(filter.df$YEAR==y))
}

# save
write.csv(mean_sdc.df, "./mean_sdc_man.csv")

# remove years with only 1 fire
mean_sdc.df <- mean_sdc.df[-c(1:2),]

# time series analysis 
myts <- ts((mean_sdc.df$SDC), start=c(1997:2001, 2003:2014, 2016:2017))
myts.log <- ts(log(mean_sdc.df$SDC), start=c(1994:2001, 2003:2014, 2016:2017))
plot(myts)
plot(myts.log)

# test for autocorrelation at first lag using Durbin-Watson test 
x <- mean_sdc.df$YEAR
y <- mean_sdc.df$SDC
mylm <- lm(y~x)
set.seed(123)
durbinWatsonTest(mylm)

# log transform to stabilize variance 
ylog <- log(y + 1)
myloglm <- lm(ylog~x)

# man-kendall test for trends - if not autocorrelated
mk.test(myts, alternative = c("two.sided"),
        continuity = TRUE)

# examine ACF and PACF plots for autocorrelation beyond first lag
acf(myts)
pacf(myts)
acf(myts.log)
pacf(myts.log)
auto.arima(myts)
auto.arima(myts.log)

###########################################################################################################
# If No Autocorrelation - Test OLS Assumptions
###########################################################################################################

# set x and y
x <- mean_sdc.df$YEAR
y <- mean_sdc.df$SDC

# mean-center for collinearity of time
time <- seq(from = 1, to = 19, by = 1)
time <- time - mean(time)
time2 <- time^2
ylog <- log(y)

# create candidate models
fit1 <- lm(y ~ x)
fit2 <- lm(ylog ~ x)
fit3 <- lm(y ~ time + time2)
fit4 <- lm(ylog ~ time + time2) #****

# test competing models
bbmle::AICctab(fit2,fit4,base=TRUE,logLik=TRUE, weights=TRUE)

# normality test
shapiro.test(resid(fit4))

# model diagnostics
autoplot(fit4)
gvlma(fit4)

# AIC/RMSE/summary
summary(fit4)
AIC(fit4, k=2)
rmse(fit4)

###########################################################################################################
# Plot Best Model
###########################################################################################################

cols <- c("MAN" = "turquoise3")
x2 <- x^2
fit3 <- lm(y ~ x + x2)
yfit <- fitted(fit3)

sdc_ts <- ggplot(data = mean_sdc.df, aes(x = x, y = y)) +
  geom_point(data = mean_sdc.df, aes(x = x, y = y, colour = "MAN"), size = 2) +
  geom_line(data = mean_sdc.df, aes(x = x, y = y), color = "turquoise3", size = 1) + 
  geom_line(data = mean_sdc.df, aes(y = yfit), color = "turquoise3", size = 1) +
  theme_bw() + 
  ggtitle("\nSDC in all Managed Fires") + 
  labs(x="\nYEAR", y = "\nln(SDC)") +
  scale_x_continuous(breaks=seq(1985,2017,5)) +
  panel_border() + 
  scale_color_manual(name = "class",
                     labels = c("MAN"),
                     values = cols)


sdc_ts <- sdc_ts + theme(plot.title = element_text(family = "Calibri", color = "black", size =16, hjust = 0.5, lineheight = 0.8, face="bold"),
                         axis.title.x = element_text(family = "Calibri", color="black", size=12, face = "bold"),
                         axis.text.x = element_text(family = "Calibri", color="black", size=12),
                         axis.title.y = element_text(family = "Calibri", color="black", size=12, face = "bold"),
                         axis.text.y = element_text(family = "Calibri", color="black", size=12),
                         #plot.margin = unit(c(1, 1.5, 1, 1.5), "cm"), #top, left, bottom, right
                         panel.border = element_rect(colour = "black", fill=NA, size=0.5))
sdc_ts

OutPath <- ("./SDC_ts_MAN.pdf")
ggsave(OutPath, plot=sdc_ts, device = cairo_pdf)
