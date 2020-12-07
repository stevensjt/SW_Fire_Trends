###########################################################################################################
# 
# This script examines proportion high-severity fire as a function of high-severity patch radius of gyration
# Author: Megan Singleton 
# Updated: 12/2/2020
#
###########################################################################################################

# load libraries 
library(ggplot2)
library(dplyr)
library(MASS)
library(bbmle)
library(car)
library(gvlma)
library(bbmle)
library(sjstats)

# set working directory 
setwd("D:/Patch_Metrics/Patch_Metrics_Complete_Analysis/R/Proportion_Regression/")

# check to see if there are any unused packages 
funchir::stale_package_check('GYRATE_AM_Proportion_HS_Regression.R')

# disable scientific notation
options(scipen=999)
options(digits=6)

# attach data
patch_metrics <- read.csv(file="GEE_1984_2017_final_patch_metrics_709.csv", header=TRUE, sep=",")
names(patch_metrics)

#create empty dataframe to input regression results in
gyrate_results.df <- data.frame(matrix(vector(), 1, 6,
                                   dimnames=list(c(), 
                                                 c("B0", "B1", "P", "RMSE", "R2", "AIC"))),
                            stringsAsFactors=F)


###########################################################################################################
# Analysis of GYRATE_AM by Proportion High-Severity - SUPPRESSION
###########################################################################################################
dev.off()
# filter by suppression fires only
filter.df <- filter(patch_metrics, FIRE_TYPE2 == "SUPPRESSION")

# set x and y 
x_gyrate <- filter.df$PLAND
y_gyrate <- filter.df$GYRATE_AM
gyrate.df <- data.frame(x_gyrate, y_gyrate)

# transform data and put into dataframe
x_logit_gyrate <- car::logit(x_gyrate)
x_log_gyrate <- log(x_gyrate)
y_log_gyrate <- log(y_gyrate)
plot(x_logit_gyrate, y_log_gyrate)

# create candidate models 
lm.gyrate <- lm(y_gyrate ~ x_gyrate, data = gyrate.df)
logit.model1 <- lm(y_gyrate ~ x_gyrate + I(x_gyrate^2), data = gyrate.df)
logit.model2 <- lm(y_log_gyrate ~ x_log_gyrate + I(x_log_gyrate^2), data = gyrate.df)
logit.model3 <- lm(y_log_gyrate ~ x_logit_gyrate, data = gyrate.df)   
logit.model4 <- lm(y_log_gyrate ~ x_logit_gyrate + I(x_logit_gyrate^2), data = gyrate.df)

# mean centered
x_logit_center <- scale(x_logit_gyrate, center = TRUE, scale = FALSE)
x_logit_center2 <- x_logit_center^2
logit.model5 <- lm(y_log_gyrate ~ x_logit_center + x_logit_center2, data = gyrate.df) #####******
summary(logit.model5)

# test models
AICctab(logit.model1, logit.model2, base = T, delta = T, weights = T)
AICctab(logit.model4, logit.model5, base = T, delta = T, weights = T)

# normality test 
shapiro.test(resid(lm.gyrate))
shapiro.test(resid(logit.model1))
shapiro.test(resid(logit.model2))
shapiro.test(resid(logit.model3))
shapiro.test(resid(logit.model4))
shapiro.test(resid(logit.model5))

# histogram
fit_gyrate <- studres(logit.model2)
hist(fit_gyrate)

# plot model diagnostics
par(mfrow = c(2, 2))  # Split the plotting panel into a 2 x 2 grid
plot(lm.gyrate)
plot(logit.model1)
plot(logit.model2)
plot(logit.model3)
plot(logit.model4)

# all diagnostics
gvlma(lm.gyrate)
gvlma(logit.model1)
gvlma(logit.model2)
gvlma(logit.model3)

# RMSE/SUMMARY/AIC
summary(logit.model2)
rmse(logit.model2)
AIC(logit.model2)

# retrieve summary statistics and put into dataframe
lm.gyrate.coef1 <- logit.model5$coefficients[1]
lm.gyrate.coef2 <- logit.model5$coefficients[2]
lm.gyrate.p <- summary(logit.model5)$coefficients[2, 4]
lm.gyrate.rmse <- rmse(logit.model5)
lm.gyrate.adr <- summary(logit.model5)$r.squared
lm.gyrate.aic <- AIC(logit.model5)

gyrate_results.df <- rbind.data.frame(gyrate_results.df, "SUPP" = c(lm.gyrate.coef1, lm.gyrate.coef2, lm.gyrate.p, 
                                                            lm.gyrate.rmse, lm.gyrate.adr, lm.gyrate.aic)) 
gyrate_results.df <- gyrate_results.df[-c(1),]; row.names(gyrate_results.df) <- 1:nrow(gyrate_results.df)
gyrate_results.df

###########################################################################################################
# Analysis of GYRATE_AM by Proportion High-Severity - Managed Fire
###########################################################################################################
dev.off()
# filter by managed fires only
filter.df1 <- filter(patch_metrics, FIRE_TYPE2 == "OTHER")

# set x and y 
x_gyrate1 <- filter.df1$PLAND
y_gyrate1 <- filter.df1$GYRATE_AM

# transform data and put into dataframe
x_logit_gyrate1 <- car::logit(x_gyrate1)
x_log_gyrate1 <- log(x_gyrate1)
y_log_gyrate1 <- log(y_gyrate1)
gyrate.df1 <- data.frame(x_gyrate1,y_log_gyrate1)
plot(x_logit_gyrate1, y_log_gyrate1)

# create candidate models 
lm.gyrate1 <- lm(y_gyrate1 ~ x_gyrate1, data = gyrate.df1)
logit.model6 <- lm(y_gyrate1 ~ x_gyrate1 + I(x_gyrate1^2), data = gyrate.df1)
logit.model7 <- lm(y_log_gyrate1 ~ x_log_gyrate1 + I(x_log_gyrate1^2), data = gyrate.df1)
logit.model8 <- lm(y_log_gyrate1 ~ x_logit_gyrate1, data = gyrate.df1)   
logit.model9 <- lm(y_log_gyrate1 ~ x_logit_gyrate1 + I(x_logit_gyrate1^2), data = gyrate.df1)

# mean centered
x_logit_center1 <- scale(x_logit_gyrate1, center = TRUE, scale = FALSE)
x_logit_center12 <- x_logit_center1^2
logit.model10 <- lm(y_log_gyrate1 ~ x_logit_center1 + x_logit_center12, data = gyrate.df1) #####******
summary(logit.model10)

# test models
AICctab(logit.model8, logit.model10, base = T, delta = T, weights = T)
AICctab(logit.model9, logit.model10, base = T, delta = T, weights = T)

# normality test 
shapiro.test(resid(lm.gyrate1))
shapiro.test(resid(logit.model6)) 
shapiro.test(resid(logit.model7))
shapiro.test(resid(logit.model8))
shapiro.test(resid(logit.model9))

# histogram
fit_gyrate1 <- studres(logit.model6)
hist(fit_gyrate1)

# plot model diagnostics
par(mfrow = c(2, 2))  # Split the plotting panel into a 2 x 2 grid
plot(lm.gyrate1)
plot(logit.model7)
plot(logit.model8)
plot(logit.model9)

# all diagnositics
gvlma(lm.gyrate1)
gvlma(logit.model7)
gvlma(logit.model8)
gvlma(logit.model9)

# RMSE/SUMMARY/AIC
summary(logit.model10)
rmse(logit.model10)
AIC(logit.model10)

# retrieve summary statistics and put into dataframe
lm.gyrate.coef1.1 <- logit.model10$coefficients[1]
lm.gyrate.coef2.1 <- logit.model10$coefficients[2]
lm.gyrate.p1 <- summary(logit.model10)$coefficients[2, 4]
lm.gyrate.rmse1 <- rmse(logit.model10)
lm.gyrate.adr1 <- summary(logit.model10)$r.squared
lm.gyrate.aic1 <- AIC(logit.model10)

gyrate_results.df1 <- rbind.data.frame(gyrate_results.df, "MAN" = c(lm.gyrate.coef1.1, lm.gyrate.coef2.1, lm.gyrate.p1, 
                                                            lm.gyrate.rmse1, lm.gyrate.adr1, lm.gyrate.aic1)) 
gyrate_results.df1

#add row names and write data table as csv
rownames(gyrate_results.df1) <- c("SUPPRESSION", "MANAGED")
gyrate_results.df1
write.csv(gyrate_results.df, "./gyrate_prop_regression.csv")
