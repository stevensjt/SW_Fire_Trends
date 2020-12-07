###########################################################################################################
# 
# This script examines high-severity class area as a function of time for suppression fires
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
funchir::stale_package_check('CA_ts_Suppression.R')

# disable scientific notation
options(scipen=999)
options(digits=6)

# attach data
ts_patch_metrics <- read.csv(file="GEE_1984_2017_final_patch_metrics_709.csv", header=TRUE, sep=",")
names(ts_patch_metrics)

###########################################################################################################
# MEAN CA by YEAR - SUPPRESSION
###########################################################################################################

# filter dataset for suppression fires only
filter.df <- filter(ts_patch_metrics, FIRE_TYPE2 == "SUPPRESSION")

# aggregate data by CA by year then calculate mean 
mean_ca.df <- aggregate(filter.df[, 6], list(filter.df[,2]), mean)
# rename columns
colnames(mean_ca.df)[colnames(mean_ca.df)=="Group.1"] <- "YEAR"
colnames(mean_ca.df)[colnames(mean_ca.df)=="x"] <- "CA"
# convert character data to numeric 
mean_ca.df <- mutate_all(mean_ca.df, function(x) as.numeric(as.character(x)))
# format to number of digits
mean_ca.df <- mean_ca.df %>% 
  mutate_if(is.numeric, round, digits = 5)
mean_ca.df

# add the number of fires per year
for(r in 1:nrow(mean_ca.df)){
  y <- mean_ca.df$YEAR[r]
  mean_ca.df$n[r] <- 
    length(which(filter.df$YEAR==y))
}

# save
write.csv(mean_ca.df, "./mean_ca_sup.csv")

# remove years with only 1 fire
mean_ca.df <- mean_ca.df[-c(1),]

# quick plot of time series
myts <- ts((mean_ca.df$CA), start=c(1985:2017))
myts.log <- ts(log(mean_ca.df$CA), start=c(1985:2017))
plot(myts)
list(myts)
plot(myts.log)
list(myts.log)

# test for lag-1 autocorrelation using Durbin Watson test 
x <- mean_ca.df$YEAR
y <- mean_ca.df$CA
ylog <- log(y)
mylm <- lm(y~x)
myloglm <- lm(ylog~x)
set.seed(123)
durbinWatsonTest(mylm)
set.seed(123)
durbinWatsonTest(myloglm)

# examine ACF and PACF plots for autocorrelation as well
acf(myts)
pacf(myts)
acf(myts.log)
pacf(myts.log)
acf(diff(myts.log))
pacf(diff(myts.log))
auto.arima(diff(myts.log))
auto.arima(myts.log)

# if no autocorrelation - try mann Kendall trend test
mk.test(myts, alternative = c("two.sided"),
        continuity = TRUE)

###########################################################################################################
# If No Autocorrelation - Test OLS Assumptions
###########################################################################################################

# set x and y 
x <- mean_ca.df$YEAR
y <- mean_ca.df$CA
ylog <- log(y)
fit1 <- lm(y ~ x)
fit2 <- lm(ylog ~ x)
x2 <- x^2
fit3 <- lm(y ~ x + x2)
time <- seq(from = 1, to = 33, by = 1)
time <- time - mean(time)
time2 <- time^2
fit4 <- lm(ylog ~ time + time2)

# normality test
shapiro.test(resid(fit4))
# all diagnositics
autoplot(fit4)
gvlma(fit4)
summary(fit4)
AIC(fit4, k=2)
rmse(fit4)

bbmle::AICctab(fit2,fit4,base=TRUE,logLik=TRUE, weights=TRUE)

###########################################################################################################
# If Autocorrelated - Use GLS Model
###########################################################################################################

# determine if time series is stationary before adding AR or MA components
# log of data is used to stabilize variance 
# differencing is used to stabilize mean
adf.test(myts)
adf.test(myts.log)

# set x and y 
x <- mean_ca.df$YEAR
y <- mean_ca.df$CA

# mean-center for collinearity of time
time <- seq(from = 1, to = 33, by = 1)
time <- time - mean(time)
time2 <- time^2

# create gls models
gls.fit1 <- gls(myts ~ x, correlation = corARMA(p=1), method = "ML") 
gls.fit2 <- gls(myts.log ~ x, correlation = corARMA(p=1), method = "ML") #****
gls.fit3 <- gls(myts ~ time + time2, correlation = corARMA(p=1), method = "ML") 
gls.fit4 <- gls(myts.log ~ time + time2, correlation = corARMA(p=1), method = "ML") #*****
gls.fit5 <- gls(myts.log ~ time + time2, correlation = corARMA(q=1), method = "ML") 

# compare models
bbmle::AICctab(gls.fit2,gls.fit4,base=TRUE,logLik=TRUE, weights=TRUE)

# summary
summary(gls.fit2)
summary(gls.fit4) #**

# liklihood ratio test 
lrtest(gls.fit2)
lrtest(gls.fit4) #**

# pseudo r squared
rcompanion::nagelkerke(gls.fit3)
rcompanion::nagelkerke(gls.fit4)

# plot residuals
plot(gls.fit2) #**
shapiro.test(resid(gls.fit2))

# goodness of fit
rmse(gls.fit2)
AIC(gls.fit2, k=2)

###########################################################################################################
# Plot Best Model -Transformed
###########################################################################################################

# using linear model since there is an upward trend 
gls.fit2 <- gls(myts.log ~ x, correlation = corARMA(p=1), method = "ML") #****
yfit <- fitted(gls.fit2)

cols <- c("SUP" = "indianred3")

ca_ts <- ggplot(data = mean_ca.df, aes(x = YEAR, y = ylog)) + 
  geom_point(color= 'indianred3', size = 3) + 
  geom_line(color = "indianred3", size = 0.7) +
  geom_line(aes(y = yfit), color = "indianred3", linetype = 1, size = 1) +
  theme_bw() + 
  ggtitle("\nHigh-Severity Area in Suppression Fires") + 
  labs(x="\nYEAR", y = "\nln(ha)") +
  scale_y_continuous(breaks=seq(1, 8, 1)) +
  scale_x_continuous(breaks=seq(1985,2017,5)) +
  #ggpubr::stat_cor(label.x = 1990, label.y = 0.20) + 
  panel_border() + 
  scale_color_manual(name = "class",
                     labels = c("SUP"))

ca_ts <- ca_ts + theme(plot.title = element_text(family = "Calibri", color = "black", size =16, hjust = 0.5, lineheight = 0.8, face="bold"),
                       axis.title.x = element_text(family = "Calibri", color="black", size=12, face = "bold"),
                       axis.text.x = element_text(family = "Calibri", color="black", size=12),
                       axis.title.y = element_text(family = "Calibri", color="black", size=12, face = "bold"),
                       axis.text.y = element_text(family = "Calibri", color="black", size=12),
                       #plot.margin = unit(c(1, 1.5, 1, 1.5), "cm"), #top, left, bottom, right
                       panel.border = element_rect(colour = "black", fill=NA, size=0.5))
ca_ts

OutPath <- "./CA_ts_suppression.pdf"
ggsave(OutPath, plot=ca_ts, device = cairo_pdf)

##################################################################################
##################################################################################