###########################################################################################################
# 
# This script examines fire size as a function of high-severity patch proportion of like adjacencies (PLADJ)
# Author: Megan Singleton 
# Updated: 12/2/2020
#
###########################################################################################################

# load libraries
library(ggplot2)
library(dplyr)
library(bbmle)
library(betareg)
library(rcompanion)

# set working directory 
setwd("D:/Patch_Metrics/Patch_Metrics_Complete_Analysis/R/Fire_Size_Regression/")

# check to see if there are any unused packages 
funchir::stale_package_check('PLADJ_Fire_Size_Regression.R')

# disable scientific notation
options(scipen=999)

# attach data
patch_metrics <- read.csv(file="GEE_1984_2017_final_patch_metrics_709.csv", header=TRUE, sep=",")
names(patch_metrics)

###########################################################################################################
# Analysis of PLADJ by Fire Size - SUPPRESSION
###########################################################################################################

dev.off()
# filter by suppression fires only
filter.df <- filter(patch_metrics, FIRE_TYPE2 == "SUPPRESSION")
names(filter.df)

# set x and y
x_pladj <- filter.df$FORESTED_BURNED
y_pladj <- filter.df$PLADJ
y_per <- y_pladj/100
pladj.df <- data.frame(x_pladj,y_pladj)
plot(pladj.df)

# create candidate models
lm.pladj <- lm(y_pladj ~ x_pladj, data = pladj.df)
lm.pladj_log_asin <- lm((asin(sqrt(y_pladj/100)) ~ log(x_pladj)), data = pladj.df)
lm.pladj_asin <- lm((asin(sqrt(y_pladj/100)) ~ x_pladj), data = pladj.df)
lm.pladj_log <- lm(y_pladj ~ log(x_pladj), data = pladj.df)
beta.pladj <- betareg(y_per ~ x_pladj, data = pladj.df)

# plot
autoplot(lm.pladj)
autoplot(lm.pladj_asin)
autoplot(lm.pladj_log)
autoplot(lm.pladj_log_asin)

# plot model diagnostics
par(mfrow = c(2, 2))  # Split the plotting panel into a 2 x 2 grid
plot(lm.pladj)
plot(lm.pladj_asin)
plot(lm.pladj_log)
plot(lm.pladj_log_asin)
plot(beta.pladj)
plot(beta.pladj, which = 1:4, type = "pearson")

# normality test 
shapiro.test(resid(lm.pladj_log))
shapiro.test(resid(lm.pladj_log_asin))
shapiro.test(resid(beta.pladj))

#likelihood ratio test for p-value
lrtest(beta.pladj)

# RMSE/SUMMARY/AIC
summary(beta.pladj) #use R2 here
rmse(beta.pladj)
AIC(beta.pladj)
rcompanion::nagelkerke(beta.pladj)

summary(beta.pladj) #use R2 here
exp(coef(beta.pladj))[2]
coef(beta.pladj)
exp(0.000009967)

###########################################################################################################
# Analysis of PLADJ by Fire Size - MANAGED
###########################################################################################################

dev.off()
# filter by managed fires only
filter.df1 <- filter(patch_metrics, FIRE_TYPE2 == "OTHER")
names(filter.df1)

# set x and y
x_pladj1 <- filter.df1$FORESTED_BURNED
y_pladj1 <- filter.df1$PLADJ
y_per1 <- y_pladj1/100
pladj.df1 <- data.frame(x_pladj1,y_pladj1)
plot(pladj.df1)

# create candidate models
lm.pladj1 <- lm(y_pladj1 ~ x_pladj1, data = pladj.df1)
lm.pladj_log_asin1 <- lm((asin(sqrt(y_pladj1/100)) ~ log(x_pladj1)), data = pladj.df1)
lm.pladj_asin1 <- lm((asin(sqrt(y_pladj1/100)) ~ x_pladj1), data = pladj.df1)
lm.pladj_log1 <- lm(y_pladj1 ~ log(x_pladj1), data = pladj.df1)
beta.pladj1 <- betareg(y_per1 ~ x_pladj1, data = pladj.df1)

# plot model diagnostics
par(mfrow = c(2, 2))  # Split the plotting panel into a 2 x 2 grid
plot(lm.pladj1)
plot(lm.pladj_asin1)
plot(lm.pladj_log1)
plot(lm.pladj_log_asin1)
plot(beta.pladj1)
plot(beta.pladj1, which = 1:4, type = "pearson")

# normality test 
shapiro.test(resid(lm.pladj_log1))
shapiro.test(resid(lm.pladj_log_asin1))
shapiro.test(resid(beta.pladj1))

#likelihood ratio test for p-value
lrtest(beta.pladj1)

# RMSE/SUMMARY/AIC
summary(beta.pladj1) #use R2 here
rmse(beta.pladj1)
AIC(beta.pladj1)

summary(beta.pladj1)
exp(coef(beta.pladj1))[2]
coef(beta.pladj1)
exp(0.000009967)


