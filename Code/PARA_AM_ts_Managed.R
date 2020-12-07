###########################################################################################################
# 
# This script examines high-severity patch perimeter:area as a function of time for managed fires
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
funchir::stale_package_check('PARA_AM_ts_Managed.R')

# disable scientific notation
options(scipen=999)
options(digits=6)

# attach data
ts_patch_metrics <- read.csv(file="GEE_1984_2017_final_patch_metrics_709.csv", header=TRUE, sep=",")
names(ts_patch_metrics)

###########################################################################################################
# MEAN PARA_AM by YEAR - MANAGED
###########################################################################################################

# filter dataset by managed fires only
filter.df <- filter(ts_patch_metrics, FIRE_TYPE2 == "OTHER")

## aggregate data by CA by year then calculate mean 
mean_para.df <- aggregate(filter.df[, 25], list(filter.df[,2]), mean)
## rename columns
colnames(mean_para.df)[colnames(mean_para.df)=="Group.1"] <- "YEAR"
colnames(mean_para.df)[colnames(mean_para.df)=="x"] <- "PARA_AM"
## convert character data to numeric 
mean_para.df <- mutate_all(mean_para.df, function(x) as.numeric(as.character(x)))
## format to number of digits
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
write.csv(mean_para.df, "./mean_para_man.csv")

# remove years with only 1 fire
mean_para.df <- mean_para.df[-c(1),]

#time series analysis 
myts <- ts(mean_para.df$PARA_AM, start=c(1995:2001, 2003:2014, 2016:2017))
myts.log <- ts(log(mean_para.df$PARA_AM), start=c(1995:2001, 2003:2014, 2016:2017))
plot(myts)
plot(myts.log)
mk.test(myts, alternative = c("two.sided"),
        continuity = TRUE)

# test if autocorrelated using Durbin Watson test 
x <- mean_para.df$YEAR
y <- mean_para.df$PARA_AM
ylog <- log(y)
mylm <- lm(y~x)
set.seed(123)
durbinWatsonTest(mylm)

# examine ACF and PACF plots for autocorrelation as well
acf(myts)
pacf(myts)
acf(myts.log)
pacf(myts.log)
acf(diff(myts.log))
pacf(diff(myts.log))
auto.arima(diff(myts.log))

###########################################################################################################
# If No Autocorrelation - Test OLS Assumptions
###########################################################################################################

# set x and y 
x <- mean_para.df$YEAR
y <- mean_para.df$PARA_AM

# create models
ylog <- log(y)
fit1 <- lm(y ~ x)
fit2 <- lm(ylog ~ x) #*****
x2 <- x^2
fit3 <- lm(y ~ x + x2)

# mean center time to mitigate collinearity
time <- seq(from = 1, to = 21, by = 1)
time <- time - mean(time)
time2 <- time^2
fit4 <- lm(ylog ~ time + time2) 

# compare models
bbmle::AICctab(fit2,fit4,base=TRUE,logLik=TRUE, weights=TRUE)

# normality test
shapiro.test(resid(fit1)) 
shapiro.test(resid(fit2))
shapiro.test(resid(fit3))
shapiro.test(resid(fit4))

# plot model diagnostics
autoplot(fit1) 
autoplot(fit2) #**
autoplot(fit3)
autoplot(fit4)

# all diagnostics
gvlma(fit1) #***
gvlma(fit2)
gvlma(fit3)
gvlma(fit4)

# summary of models
summary(fit1) 
summary(fit2) #***
summary(fit3)
summary(fit4)

# RMSE/AIC/SUMMARY
summary(fit2)
rmse(fit2)
AIC(fit2, k=2)


###########################################################################################################
# Plot Best Model
###########################################################################################################

fit2 <- lm(ylog ~ x) #*****
yfit <- fitted(fit2)
cols <- c("MAN" = "turquoise3")

para_ts <- ggplot(data = mean_para.df, aes(x = YEAR, y = ylog)) + 
  geom_point(color= 'turquoise3', size = 3) + 
  geom_line(color = "turquoise3", size = 0.7) +
  geom_line(aes(y = yfit), color = "turquoise3", linetype = 1, size = 1) +
  #geom_smooth(method = "lm", se = TRUE, color = 'turquoise3', size = 1.5) + 
  theme_bw() + 
  ggtitle("\nArea-Weighted Mean Perimeter:Area in all Managed Fires") + 
  labs(x="\nYEAR", y = "\nPARA_AM") +
  #scale_y_continuous(breaks=seq(200, 800, 150)) +
  scale_x_continuous(breaks=seq(1985,2017,5)) +
  #ggpubr::stat_cor(label.x = 1990, label.y = 0.20) + 
  panel_border() + 
  scale_color_manual(name = "class",
                     labels = c("MAN"),
                     values = cols) 

para_ts <- para_ts + theme(plot.title = element_text(family = "Calibri", color = "black", size =16, hjust = 0.5, lineheight = 0.8, face="bold"),
                           axis.title.x = element_text(family = "Calibri", color="black", size=12, face = "bold"),
                           axis.text.x = element_text(family = "Calibri", color="black", size=12),
                           axis.title.y = element_text(family = "Calibri", color="black", size=12, face = "bold"),
                           axis.text.y = element_text(family = "Calibri", color="black", size=12),
                           #plot.margin = unit(c(1, 1.5, 1, 1.5), "cm"), #top, left, bottom, right
                           panel.border = element_rect(colour = "black", fill=NA, size=0.5))
para_ts

OutPath <- "./PARA_AM_ts_man.pdf"
ggsave(OutPath, plot=para_ts, device = cairo_pdf)
