###########################################################################################################
# 
# This script examines proportion high-severity fire as a function of high-severity class area
# Author: Megan Singleton 
# Updated: 11/30/2020
#
###########################################################################################################

# load libraries 
library(ggplot2)
library(dplyr)
library(MASS)
library(bbmle)
library(car)
library(bbmle)
library(gvlma)
library(sjstats)

# set working directory 
setwd("D:/Patch_Metrics/Patch_Metrics_Complete_Analysis/R/Proportion_Regression/")

# check to see if there are any unused packages 
funchir::stale_package_check('CA_Proportion_HS_Regression.R')

# disable scientific notation
options(scipen=999)
options(digits=6)

# attach data
patch_metrics <- read.csv(file="GEE_1984_2017_final_patch_metrics_709.csv", header=TRUE, sep=",")
names(patch_metrics)

#create empty dataframe to input regression results in
ca_results.df <- data.frame(matrix(vector(), 1, 6,
                                   dimnames=list(c(), 
                                                 c("B0", "B1", "P", "RMSE", "R2", "AIC"))),
                            stringsAsFactors=F)


###########################################################################################################
# Analysis of CA by Proportion High-Severity - SUPPRESSION
###########################################################################################################

# filter by suppression fires only
filter.df <- filter(patch_metrics, FIRE_TYPE2 == "SUPPRESSION")

# set x and y 
x_ca <- filter.df$PLAND
y_ca <- filter.df$CA

# transform data and put into dataframe
x_logit_ca <- car::logit(x_ca)
y_log_ca <- log(y_ca)
ca.df <- data.frame(x_logit_ca, y_log_ca)
plot(x_logit_ca, y_log_ca)

# create candidate models
lm.ca <- lm(y_ca ~ x_ca, data = ca.df)
logit.model1 <- lm(y_log_ca ~ x_ca + I(x_ca^2), data = ca.df)
logit.model2 <- lm(y_log_ca ~ x_logit_ca, data = ca.df)             #####******      
logit.model3 <- lm(y_log_ca ~ x_logit_ca + I(x_logit_ca^2), data = ca.df) 

# test models
AICctab(logit.model1, logit.model2, base = T, delta = T, weights = T)
AICctab(logit.model2, logit.model3, base = T, delta = T, weights = T)

# normality test 
shapiro.test(resid(lm.ca))
shapiro.test(resid(logit.model1))
shapiro.test(resid(logit.model2))
shapiro.test(resid(logit.model3))

# histogram
fit_ca <- studres(logit.model2)
hist(fit_ca)

# plot model diagnostics
par(mfrow = c(2, 2))  # Split the plotting panel into a 2 x 2 grid
plot(lm.ca)
plot(logit.model1)
plot(logit.model2) #***
plot(logit.model3)

# all diagnostics
gvlma(lm.ca)
gvlma(logit.model1)
gvlma(logit.model2)
gvlma(logit.model3)

# RMSE/SUMMARY/AIC
summary(logit.model2)
rmse(logit.model2)
AIC(logit.model2)

# retrieve summary statistics and put into dataframe
lm.ca.coef1 <- logit.model2$coefficients[1]
lm.ca.coef2 <- logit.model2$coefficients[2]
lm.ca.p <- summary(logit.model2)$coefficients[2, 4]
lm.ca.rmse <- rmse(logit.model2)
lm.ca.adr <- summary(logit.model2)$r.squared
lm.ca.aic <- AIC(logit.model2)

ca_results.df <- rbind.data.frame(ca_results.df, "SUPP" = c(lm.ca.coef1, lm.ca.coef2, lm.ca.p, 
                                                            lm.ca.rmse, lm.ca.adr, lm.ca.aic)) 
ca_results.df <- ca_results.df[-c(1),]; row.names(ca_results.df) <- 1:nrow(ca_results.df)
ca_results.df

###########################################################################################################
# Analysis of CA by Proportion High-Severity - Managed Fire
###########################################################################################################
dev.off()
# filter by managed fires only
filter.df1 <- filter(patch_metrics, FIRE_TYPE2 == "OTHER")

# set x and y 
x_ca1 <- filter.df1$PLAND
y_ca1 <- filter.df1$CA

# transform data and put into dataframe
x_logit_ca1 <- car::logit(x_ca1)
y_log_ca1 <- log(y_ca1)
ca.df1 <- data.frame(x_ca1,y_log_ca1)
plot(x_logit_ca1, y_log_ca1)

# create candidate models 
lm.ca1 <- lm(y_ca1 ~ x_ca1, data = ca.df1)
logit.model5 <- lm(y_log_ca1 ~ x_ca1 + I(x_ca1^2), data = ca.df1)
logit.model6 <- lm(y_log_ca1 ~ x_logit_ca1, data = ca.df1)                   #####******
logit.model7 <- lm(y_log_ca1 ~ x_logit_ca1 + I(x_logit_ca1^2), data = ca.df1) 

# test models
AICctab(logit.model5, logit.model6, base = T, delta = T, weights = T)
AICctab(logit.model6, logit.model7, base = T, delta = T, weights = T)

# normality test 
shapiro.test(resid(lm.ca1))
shapiro.test(resid(logit.model5))
shapiro.test(resid(logit.model6)) 
shapiro.test(resid(logit.model7))

# histogram
fit_ca1 <- studres(logit.model6)
hist(fit_ca1)

# plot model diagnostics
par(mfrow = c(2, 2))  # Split the plotting panel into a 2 x 2 grid
plot(lm.ca1)
plot(logit.model5)
plot(logit.model6)
plot(logit.model7)

# all diagnositics
gvlma(lm.ca1)
gvlma(logit.model5)
gvlma(logit.model6)
gvlma(logit.model7)

# RMSE/SUMMARY/AIC
summary(logit.model6)
rmse(logit.model6)
AIC(logit.model6)

# retrieve summary statistics and put into dataframe
lm.ca.coef1.1 <- logit.model6$coefficients[1]
lm.ca.coef2.1 <- logit.model6$coefficients[2]
lm.ca.p1 <- summary(logit.model6)$coefficients[2, 4]
lm.ca.rmse1 <- rmse(logit.model6)
lm.ca.adr1 <- summary(logit.model6)$r.squared
lm.ca.aic1 <- AIC(logit.model6)

ca_results.df1 <- rbind.data.frame(ca_results.df, "MAN" = c(lm.ca.coef1.1, lm.ca.coef2.1, lm.ca.p1, 
                                                            lm.ca.rmse1, lm.ca.adr1, lm.ca.aic1)) 
ca_results.df1

#add row names and write data table as csv
rownames(ca_results.df1) <- c("SUPPRESSION", "MANAGED")
ca_results.df1
write.csv(ca_results.df, "./ca_regression_results_sup_man.csv")
