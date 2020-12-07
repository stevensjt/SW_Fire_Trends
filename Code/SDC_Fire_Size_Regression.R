###########################################################################################################
# 
# This script examines fire size as a function of the high-severity stand decay coefficient
# Author: Megan Singleton 
# Updated: 12/2/2020
#
###########################################################################################################

# load libraries
library(bbmle)
library(car) 
library(dplyr)
library(EnvStats)
library(gvlma)
library(MASS) 
library(sjstats)
library(extrafont)
loadfonts()

# set working directory 
setwd("D:/Patch_Metrics/Patch_Metrics_Complete_Analysis/R/Fire_Size_Regression/")

# check to see if there are any unused packages 
funchir::stale_package_check('SDC_Fire_Size_Regression.R')

# disable scientific notation
options(scipen=999)

# attach data
patch_metrics <- read.csv(file="GEE_1984_2017_final_patch_metrics_709.csv", header=TRUE, sep=",")
names(patch_metrics)

#create empty dataframe to input regression results in
sdc_results.df <- data.frame(matrix(vector(), 1, 6,
                                     dimnames=list(c(), 
                                                   c("B0", "B1", "P", "RMSE", "R2", "AIC"))),
                              stringsAsFactors=F)

###########################################################################################################
# Analysis of SDC_AM by Fire Size - SUPPRESSION
###########################################################################################################
dev.off()
# filter by suppression fires only
filter.df <- filter(patch_metrics, FIRE_TYPE2 == "SUPPRESSION",SDC >0)
names(filter.df)

# set x and y 
x_sdc <- filter.df$FORESTED_BURNED
y_sdc <- filter.df$SDC
sdc.df <- data.frame(x_sdc,y_sdc)
plot(x_sdc, y_sdc)

# log transform
x_sdc_log <- log(x_sdc)
y_sdc_log <- log(y_sdc)
plot(x_sdc_log, y_sdc_log)

# create linear model
lm.sdc <- lm(y_sdc ~ x_sdc, data = sdc.df)
summary(lm.sdc)

# test for appropriate transformation
car::boxCox(lm.sdc)

# create candidate models
lm.sdc_log_x <- lm(y_sdc ~ log(x_sdc), data = sdc.df)
lm.sdc_log_y <- lm(log(y_sdc) ~ x_sdc, data = sdc.df)
lm.sdc_log_log <- lm(log(y_sdc) ~ log(x_sdc), data = sdc.df)
lm.sdc_sqrt_log <- lm(sqrt(y_sdc) ~ log(x_sdc), data = sdc.df)
lm.sdc_sqrt <- lm(sqrt(y_sdc) ~ sqrt(x_sdc), data = sdc.df)
x2 <- x_sdc^2
lm.sdc_quad <- lm(y_sdc ~ x_sdc + x2, data = sdc.df)
lm.sdc_log_quad <- lm(log(y_sdc) ~ x_sdc + x2, data = sdc.df)
lm.sdc_log_log_quad <- lm(log(y_sdc) ~ log(x_sdc) + log(x2), data = sdc.df)
lm.sdc_sqrt_log <- lm(sqrt(y_sdc) ~ log(x_sdc), data = sdc.df)

# test models
AICctab(lm.sdc_log_log, lm.sdc_log_log_quad, base = T, delta = T, weights = T)

# plot model diagnostics
par(mfrow = c(2, 2))  # Split the plotting panel into a 2 x 2 grid
plot(lm.sdc)
plot(lm.sdc_log_y)
plot(lm.sdc_log_log)
plot(lm.sdc_sqrt_log)

# normality test 
shapiro.test(resid(lm.sdc))
shapiro.test(resid(lm.sdc_log_x))
shapiro.test(resid(lm.sdc_log_log))
shapiro.test(resid(lm.sdc_sqrt_log))
fit_sdc <- studres(lm.sdc_log_log)
fit_sdc <- studres(lm.sdc_sqrt_log)
hist(fit_sdc)

# testing homoscedasticity
ncvTest(lm.sdc_log_x)
ncvTest(lm.sdc_log_log)

# all diagnostics
gvlma(lm.sdc)
gvlma(lm.sdc_log_y)
gvlma(lm.sdc_log_log)
gvlma(lm.sdc_sqrt_log)

# RMSE/SUMMARY/AIC
summary(lm.sdc_log_log)
rmse(lm.sdc_log_log)
AIC(lm.sdc_log_log)

# retrieve summary statistics and put into dataframe
lm.sdc.coef1 <- lm.sdc_log_log$coefficients[1]
lm.sdc.coef2 <- lm.sdc_log_log$coefficients[2]
lm.sdc.p <- summary(lm.sdc_log_log)$coefficients[2, 4]
lm.sdc.rmse <- rmse(lm.sdc_log_log)
lm.sdc.adr <- summary(lm.sdc_log_log)$r.squared
lm.sdc.aic <- AIC(lm.sdc_log_log)

sdc_results.df <- rbind.data.frame(sdc_results.df, "SUPP" = c(lm.sdc.coef1, lm.sdc.coef2, lm.sdc.p, 
                                                                lm.sdc.rmse, lm.sdc.adr, lm.sdc.aic)) 
sdc_results.df <- sdc_results.df[-c(1),]; row.names(sdc_results.df) <- 1:nrow(sdc_results.df)
sdc_results.df
summary(lm.sdc_log_log)
summary(lm.sdc_sqrt_log)

###########################################################################################################
# Analysis of SDC_AM by Fire Size - MANAGED
###########################################################################################################

# filter by managed fires only
filter.df1 <- filter(patch_metrics, FIRE_TYPE2 == "OTHER", SDC > 0)
names(filter.df1)

# set x and y 
x_sdc1 <- filter.df1$FORESTED_BURNED
y_sdc1 <- filter.df1$SDC

# create new dataframe and plot to get general idea of data
sdc.df1 <- data.frame(x_sdc1,y_sdc1)
dev.off()
plot(sdc.df1)

# create linear model
lm.sdc1 <- lm(y_sdc1 ~ x_sdc1, data = sdc.df1)
summary(lm.sdc1)

# test for appropriate transformation
# lambda close to 0 suggest log transformation
car::boxCox(lm.sdc1)

# create candidate models
lm.sdc_log_x1 <- lm(y_sdc1 ~ log(x_sdc1), data = sdc.df1)
lm.sdc_log_y1 <- lm(log(y_sdc1) ~ x_sdc1, data = sdc.df1)
lm.sdc_log_log1 <- lm(log(y_sdc1) ~ log(x_sdc1), data = sdc.df1)
lm.sdc_sqrt1 <- lm(sqrt(y_sdc1) ~ x_sdc1, data = sdc.df1)
lm.sdc_sqrt_sqrt1 <- lm(sqrt(y_sdc1) ~ sqrt(x_sdc1), data = sdc.df1)
lm.sdc_log_sqrt1 <- lm(log(y_sdc1) ~ sqrt(x_sdc1), data = sdc.df1)
lm.sdc_sqrt_log1 <- lm(sqrt(y_sdc1) ~ log(x_sdc1), data = sdc.df1)

# plot model diagnostics
par(mfrow = c(2, 2))  # Split the plotting panel into a 2 x 2 grid
plot(lm.sdc1)
plot(lm.sdc_log_x1)
plot(lm.sdc_log_log1)
plot(lm.sdc_sqrt1)
plot(lm.sdc_sqrt_sqrt1)
plot(lm.sdc_log_sqrt1)
plot(lm.sdc_sqrt_log1)

# normality test 
shapiro.test(resid(lm.sdc1))
shapiro.test(resid(lm.sdc_log_x1))
shapiro.test(resid(lm.sdc_log_log1))
shapiro.test(resid(lm.sdc_sqrt_sqrt1))
shapiro.test(resid(lm.sdc_sqrt1))
shapiro.test(resid(lm.sdc_sqrt_log1))
fit_sdc1 <- studres(lm.sdc_log_log1)
fit_sdc1 <- studres(lm.sdc_sqrt1)
hist(fit_sdc1)

# testing homoscedasticity
ncvTest(lm.sdc_log_x1)
ncvTest(lm.sdc_log_log1)

# all diagnositics
gvlma(lm.sdc_log_log1)
gvlma(lm.sdc_sqrt_log1)

# RMSE/SUMMARY/AIC
summary(lm.sdc_log_log1)
rmse(lm.sdc_log_log1)
AIC(lm.sdc_log_log1)

# retrieve summary statistics and put into dataframe
lm.sdc.coef1.1 <- lm.sdc_log_log1$coefficients[1]
lm.sdc.coef2.1 <- lm.sdc_log_log1$coefficients[2]
lm.sdc.p1 <- summary(lm.sdc_log_log1)$coefficients[2, 4]
lm.sdc.rmse1 <- rmse(lm.sdc_log_log1)
lm.sdc.adr1 <- summary(lm.sdc_log_log1)$r.squared
lm.sdc.aic1 <- AIC(lm.sdc_log_log1)

sdc_results.df1 <- rbind.data.frame(sdc_results.df, "MAN" = c(lm.sdc.coef1.1, lm.sdc.coef2.1, lm.sdc.p1, 
                                                                lm.sdc.rmse1, lm.sdc.adr1, lm.sdc.aic1))
summary(lm.sdc_sqrt_log1)
summary(lm.sdc_log_log1)

#add row names and write data table as csv
rownames(sdc_results.df1) <- c("SUPPRESSION", "MANAGED")
sdc_results.df1
write.csv(sdc_results.df, "./sdc_regression_results.csv")

###########################################################################################################
###########################################################################################################
