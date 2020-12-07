###########################################################################################################
# 
# This script examines high-severity patch perimeter:area as a function of time for suppression fires
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
funchir::stale_package_check('PARA_AM_ts_Suppression.R')

# disable scientific notation
options(scipen=999)
options(digits=6)

# attach data
ts_patch_metrics <- read.csv(file="GEE_1984_2017_final_patch_metrics_709.csv", header=TRUE, sep=",")
names(ts_patch_metrics)

###########################################################################################################
# MEAN PARA_AM by YEAR - SUPPRESSION
###########################################################################################################

# filter by supppression fires only
filter.df <- filter(ts_patch_metrics, FIRE_TYPE2 == "SUPPRESSION")

##aggregate data by CA by year then calculate mean 
mean_para.df <- aggregate(filter.df[, 25], list(filter.df[,2]), mean)
##rename columns
colnames(mean_para.df)[colnames(mean_para.df)=="Group.1"] <- "YEAR"
colnames(mean_para.df)[colnames(mean_para.df)=="x"] <- "PARA_AM"
##convert character data to numeric 
mean_para.df <- mutate_all(mean_para.df, function(x) as.numeric(as.character(x)))
##format to number of digits
mean_para.df <- mean_para.df %>% 
  mutate_if(is.numeric, round, digits = 5)
mean_para.df

#number of fires by year
for(r in 1:nrow(mean_para.df)){
  y <- mean_para.df$YEAR[r]
  mean_para.df$n[r] <- 
    length(which(filter.df$YEAR==y))
}

# save 
write.csv(mean_para.df, "./mean_para_sup.csv")

# remove years with one fire
mean_para.df <- mean_para.df[-c(1),]

#time series analysis 
myts <- ts(mean_para.df$PARA_AM, start=1985, end= 2017)
myts.log <- ts(log(mean_para.df$PARA_AM), start=1985, end= 2017)
plot(myts)
plot(myts.log)

# if no autocorrelation - Mann Kendall trend test
mk.test(myts, alternative = c("two.sided"),
        continuity = TRUE)

# test if autocorrelated using Durbin-Watson test 
x <- mean_para.df$YEAR
y <- mean_para.df$PARA_AM
ylog <- log(y)
mylm <- lm(y~x)
set.seed(123)
durbinWatsonTest(mylm)

#examine ACF and PACF plots for autocorrelation as well
acf(myts)
pacf(myts)
acf(myts.log)
pacf(myts.log)
auto.arima(myts)

###########################################################################################################
# If No Autocorrelation - Test OLS Assumptions
###########################################################################################################

x <- mean_para.df$YEAR
y <- mean_para.df$PARA_AM
ylog <- log(y)
fit1 <- lm(y ~ x)
fit2 <- lm(ylog ~ x)
x2 <- x^2
fit3 <- lm(y ~ x + x2)
time <- seq(from = 1, to = 33, by = 1)
time <- time - mean(time)
time2 <- time^2
fit4 <- lm(ylog ~ time + time2) #*****

#normality test
shapiro.test(resid(fit4))
#all diagnositics
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
x_para <- mean_para.df$YEAR
y_para <- mean_para.df$PARA_AM

# mean-center to mitigate collinearity of time
time <- seq(from = 1, to = 33, by = 1)
time <- time - mean(time)
time2 <- time^2

# create gls models
gls.fit1 <- gls(myts.log ~ x, correlation = corARMA(p=1), method = "ML")
gls.fit2 <- gls(myts.log ~ x, correlation = corARMA(q=1), method = "ML")
gls.fit3 <- gls(myts.log ~ time + time2, correlation = corARMA(p=1), method = "ML")
gls.fit4 <- gls(myts.log ~ time + time2, correlation = corARMA(q=1), method = "ML") #***
gls.fit5 <- gls(myts ~ time + time2, correlation = corARMA(q=1), method = "ML")

# compare models
bbmle::AICctab(gls.fit3,gls.fit4,base=TRUE,logLik=TRUE, weights=TRUE)

# summary
summary(gls.fit4) #**
summary(gls.fit5) 

# liklihood ratio test 
lrtest(gls.fit3)
lrtest(gls.fit4) #**
lrtest(gls.fit5) 

# pseudo r squared
rcompanion::nagelkerke(gls.fit3)
rcompanion::nagelkerke(gls.fit4) #**

# plot residuals
plot(gls.fit3) 
plot(gls.fit4) #**
plot(gls.fit5) 
shapiro.test(resid(gls.fit4))

# goodness of fit
rmse(gls.fit4)
AIC(gls.fit4, k=2)


###########################################################################################################
# Plot Best Model 
###########################################################################################################
yfit <- fitted(gls.fit4)
cols <- c("SUP" = "indianred3")

para_ts <- ggplot(data = mean_para.df, aes(x = x_para, y = log(y_para))) + 
  geom_point(color= 'indianred3', size = 3) + 
  geom_line(color = "indianred3", size = 0.7) +
  geom_line(aes(y = yfit), color = "indianred3", linetype = 1, size = 1) +
  theme_bw() + 
  ggtitle("\nArea-Weighted Mean Perimeter:Area in all Suppression Fires") + 
  labs(x="\nYEAR", y = "\nPARA_AM") +
  #scale_y_continuous(breaks=seq(200, 800, 150)) +
  scale_x_continuous(breaks=seq(1985,2017,5)) +
  panel_border()

para_ts <- para_ts + theme(plot.title = element_text(family = "Calibri", color = "black", size =16, hjust = 0.5, lineheight = 0.8, face="bold"),
                           axis.title.x = element_text(family = "Calibri", color="black", size=12, face = "bold"),
                           axis.text.x = element_text(family = "Calibri", color="black", size=12),
                           axis.title.y = element_text(family = "Calibri", color="black", size=12, face = "bold"),
                           axis.text.y = element_text(family = "Calibri", color="black", size=12),
                           plot.margin = unit(c(1, 1.5, 1, 1.5), "cm"), #top, left, bottom, right
                           panel.border = element_rect(colour = "black", fill=NA, size=0.5))
para_ts

OutPath <- "./PARA_AM_ts_supp.pdf"
ggsave(OutPath, plot=para_ts, device = cairo_pdf)
