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
library(MASS)
library(bbmle)
library(sjstats)
library(betareg)
library(bbmle)
library(lmtest)

# set working directory 
setwd("D:/Patch_Metrics/Patch_Metrics_Complete_Analysis/R/Proportion_Regression/")

# check to see if there are any unused packages 
funchir::stale_package_check('PLADJ_Proportion_HS_Regression.R')

# disable scientific notation
options(scipen=999)

# attach data
patch_metrics <- read.csv(file="GEE_1984_2017_final_patch_metrics_709.csv", header=TRUE, sep=",")
names(patch_metrics)

###########################################################################################################
# Analysis of PLADJ by Proportino High-Severity - SUPPRESSION
###########################################################################################################

dev.off()
# filter by suppression fires only
filter.df <- filter(patch_metrics, FIRE_TYPE2 == "SUPPRESSION")
names(filter.df)

# set x and y
x_pladj <- filter.df$PLAND
y_pladj <- filter.df$PLADJ
y_per <- y_pladj/100
pladj.df <- data.frame(x_pladj,y_pladj)
plot(pladj.df)

# create candidate models
lm.pladj <- lm(y_pladj ~ x_pladj, data = pladj.df)
lm.pladj_log_asin <- lm((asin(sqrt(y_pladj/100)) ~ log(x_pladj)), data = pladj.df)
lm.pladj_asin <- lm((asin(sqrt(y_pladj/100)) ~ x_pladj), data = pladj.df)
lm.pladj_log <- lm(y_pladj ~ log(x_pladj), data = pladj.df)
beta.pladj1 <- betareg(y_per ~ x_pladj, data = pladj.df, link = "logit") 
beta.pladj2 <- betareg(y_per ~ x_pladj, data = pladj.df, link = "loglog") #####******
beta.pladj3 <- betareg(y_per ~ x_pladj, data = pladj.df, link = "probit") 

# test models
AICctab(beta.pladj1, beta.pladj2, base = T, delta = T, weights = T)

# plot model diagnostics
par(mfrow = c(2, 2))  # Split the plotting panel into a 2 x 2 grid
plot(lm.pladj)
plot(lm.pladj_asin)
plot(lm.pladj_log)
plot(lm.pladj_log_asin)
plot(beta.pladj2)
plot(beta.pladj, which = 1:4, type = "pearson")

#likelihood ratio test for p-value
lrtest(beta.pladj2)

# RMSE/SUMMARY/AIC
summary(beta.pladj2) #use R2 here
rmse(beta.pladj2)
AIC(beta.pladj2)
rcompanion::nagelkerke(beta.pladj2)

summary(beta.pladj2) #use R2 here
exp(coef(beta.pladj2))[2]
coef(beta.pladj2)

###########################################################################################################
# Analysis of PLADJ by Fire Size - MANAGED
###########################################################################################################

dev.off()
# filter by managed fires only
filter.df1 <- filter(patch_metrics, FIRE_TYPE2 == "OTHER")
names(filter.df1)

# set x and y
x_pladj1 <- filter.df1$PLAND
y_pladj1 <- filter.df1$PLADJ
y_per1 <- y_pladj1/100
pladj.df1 <- data.frame(x_pladj1,y_pladj1)
plot(pladj.df1)

# create candidate models
lm.pladj1 <- lm(y_pladj1 ~ x_pladj1, data = pladj.df1)
lm.pladj_log_asin1 <- lm((asin(sqrt(y_pladj1/100)) ~ log(x_pladj1)), data = pladj.df1)
lm.pladj_asin1 <- lm((asin(sqrt(y_pladj1/100)) ~ x_pladj1), data = pladj.df1)
lm.pladj_log1 <- lm(y_pladj1 ~ log(x_pladj1), data = pladj.df1)
beta.pladj4 <- betareg(y_per1 ~ x_pladj1, data = pladj.df, link = "logit") 
beta.pladj5 <- betareg(y_per1 ~ x_pladj1, data = pladj.df, link = "loglog") #####******
beta.pladj6 <- betareg(y_per1 ~ x_pladj1, data = pladj.df, link = "probit") 

# plot model diagnostics
par(mfrow = c(2, 2))  # Split the plotting panel into a 2 x 2 grid
plot(lm.pladj1)
plot(lm.pladj_asin1)
plot(lm.pladj_log1)
plot(lm.pladj_log_asin1)
plot(beta.pladj5)
plot(beta.pladj1, which = 1:4, type = "pearson")

#likelihood ratio test for p-value
lrtest(beta.pladj5)

# test models
AICctab(beta.pladj5, beta.pladj6, base = T, delta = T, weights = T)

# RMSE/SUMMARY/AIC
summary(beta.pladj5) #use R2 here
rmse(beta.pladj5)
AIC(beta.pladj5)

exp(coef(beta.pladj5))[2]
coef(beta.pladj5)


