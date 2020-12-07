###########################################################################################################
# 
# This script examines proportion high-severity fire as a function of high-severity stand decay coefficent
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
library(bbmle)
library(gvlma)
library(sjstats)

# set working directory 
setwd("D:/Patch_Metrics/Patch_Metrics_Complete_Analysis/R/Proportion_Regression/")

# check to see if there are any unused packages 
funchir::stale_package_check('SDC_Proportion_HS_Regression.R')

# disable scientific notation
options(scipen=999)
options(digits=6)

# attach data
patch_metrics <- read.csv(file="GEE_1984_2017_final_patch_metrics_709.csv", header=TRUE, sep=",")
names(patch_metrics)

#create empty dataframe to input regression results in
sdc_results.df <- data.frame(matrix(vector(), 1, 6,
                                       dimnames=list(c(), 
                                                     c("B0", "B1", "P", "RMSE", "R2", "AIC"))),
                                stringsAsFactors=F)


###########################################################################################################
# Analysis of SDC by Proportion High-Severity - SUPPRESSION
###########################################################################################################
dev.off()
# filter by suppression fires only
filter.df <- filter(patch_metrics, FIRE_TYPE2 == "SUPPRESSION", SDC > 0)

# set x and y 
x_sdc <- filter.df$PLAND
y_sdc <- filter.df$SDC
sdc.df <- data.frame(x_sdc, y_sdc)

# transform data and put into dataframe
x_logit_sdc <- car::logit(x_sdc)
x_log_sdc <- log(x_sdc)
y_log_sdc <- log(y_sdc)
plot(x_log_sdc, y_log_sdc)

# create candidate models 
lm.sdc <- lm(y_sdc ~ x_sdc, data = sdc.df)
logit.model1 <- lm(y_sdc ~ x_sdc + I(x_sdc^2), data = sdc.df)
logit.model2 <- lm(y_log_sdc ~ x_log_sdc + I(x_log_sdc^2), data = sdc.df)
logit.model3 <- lm(y_log_sdc ~ x_logit_sdc, data = sdc.df)   
logit.model4 <- lm(y_log_sdc ~ x_logit_sdc + I(x_logit_sdc^2), data = sdc.df)

# mean centered
x_logit_center <- scale(x_logit_sdc, center = TRUE, scale = FALSE)
x_logit_center2 <- x_logit_center^2
logit.model5 <- lm(y_log_sdc ~ x_logit_center + x_logit_center2, data = sdc.df) #####******
summary(logit.model5)

# test models
AICctab(logit.model4, logit.model5, base = T, delta = T, weights = T)

# normality test 
shapiro.test(resid(lm.sdc))
shapiro.test(resid(logit.model1))
shapiro.test(resid(logit.model2))
shapiro.test(resid(logit.model3))
shapiro.test(resid(logit.model4))
shapiro.test(resid(logit.model5))

# histogram
fit_sdc <- studres(logit.model5)
hist(fit_sdc)

# plot model diagnostics
par(mfrow = c(2, 2))  # Split the plotting panel into a 2 x 2 grid
plot(lm.sdc)
plot(logit.model1)
plot(logit.model2)
plot(logit.model3)
plot(logit.model4)
plot(logit.model5)

# all diagnostics
gvlma(lm.sdc)
gvlma(logit.model3)
gvlma(logit.model4)
gvlma(logit.model5)

# RMSE/SUMMARY/AIC
summary(logit.model5)
rmse(logit.model5)
AIC(logit.model5)

# retrieve summary statistics and put into dataframe
lm.sdc.coef1 <- logit.model5$coefficients[1]
lm.sdc.coef2 <- logit.model5$coefficients[2]
lm.sdc.p <- summary(logit.model5)$coefficients[2, 4]
lm.sdc.rmse <- rmse(logit.model5)
lm.sdc.adr <- summary(logit.model5)$r.squared
lm.sdc.aic <- AIC(logit.model5)

sdc_results.df <- rbind.data.frame(sdc_results.df, "SUPP" = c(lm.sdc.coef1, lm.sdc.coef2, lm.sdc.p, 
                                                                    lm.sdc.rmse, lm.sdc.adr, lm.sdc.aic)) 
sdc_results.df <- sdc_results.df[-c(1),]; row.names(sdc_results.df) <- 1:nrow(sdc_results.df)
sdc_results.df

###########################################################################################################
# Analysis of SDC by Proportion High-Severity - Managed Fire
###########################################################################################################
dev.off()
# filter by managed fires only
filter.df1 <- filter(patch_metrics, FIRE_TYPE2 == "OTHER", SDC > 0)

# set x and y 
x_sdc1 <- filter.df1$PLAND
y_sdc1 <- filter.df1$SDC

# transform data and put into dataframe
x_logit_sdc1 <- car::logit(x_sdc1)
x_log_sdc1 <- log(x_sdc1)
y_log_sdc1 <- log(y_sdc1)
sdc.df1 <- data.frame(x_sdc1,y_log_sdc1)
plot(x_logit_sdc1, y_log_sdc1)

# create candidate models 
lm.sdc1 <- lm(y_sdc1 ~ x_sdc1, data = sdc.df1)
logit.model6 <- lm(y_sdc1 ~ x_sdc1 + I(x_sdc1^2), data = sdc.df1)
logit.model7 <- lm(y_log_sdc1 ~ x_log_sdc1 + I(x_log_sdc1^2), data = sdc.df1)
logit.model8 <- lm(y_log_sdc1 ~ x_logit_sdc1, data = sdc.df1)   
logit.model9 <- lm(y_log_sdc1 ~ x_logit_sdc1 + I(x_logit_sdc1^2), data = sdc.df1)

# mean centered
x_logit_center1 <- scale(x_logit_sdc1, center = TRUE, scale = FALSE)
x_logit_center12 <- x_logit_center1^2
logit.model10 <- lm(y_log_sdc1 ~ x_logit_center1 + x_logit_center12, data = sdc.df1) #####******
summary(logit.model10)

# test models
AICctab(logit.model8, logit.model10, base = T, delta = T, weights = T)
AICctab(logit.model9, logit.model10, base = T, delta = T, weights = T)

# normality test 
shapiro.test(resid(lm.sdc1))
shapiro.test(resid(logit.model6)) 
shapiro.test(resid(logit.model7))
shapiro.test(resid(logit.model8))
shapiro.test(resid(logit.model9))

# histogram
fit_sdc1 <- studres(logit.model9)
hist(fit_sdc1)

# plot model diagnostics
par(mfrow = c(2, 2))  # Split the plotting panel into a 2 x 2 grid
plot(lm.sdc1)
plot(logit.model7)
plot(logit.model8)
plot(logit.model9)

# all diagnositics
gvlma(lm.sdc1)
gvlma(logit.model7)
gvlma(logit.model8)
gvlma(logit.model9)

# RMSE/SUMMARY/AIC
summary(logit.model10)
rmse(logit.model10)
AIC(logit.model10)

# retrieve summary statistics and put into dataframe
lm.sdc.coef1.1 <- logit.model10$coefficients[1]
lm.sdc.coef2.1 <- logit.model10$coefficients[2]
lm.sdc.p1 <- summary(logit.model10)$coefficients[2, 4]
lm.sdc.rmse1 <- rmse(logit.model10)
lm.sdc.adr1 <- summary(logit.model10)$r.squared
lm.sdc.aic1 <- AIC(logit.model10)

sdc_results.df1 <- rbind.data.frame(sdc_results.df, "MAN" = c(lm.sdc.coef1.1, lm.sdc.coef2.1, lm.sdc.p1, 
                                                                    lm.sdc.rmse1, lm.sdc.adr1, lm.sdc.aic1)) 
sdc_results.df1

#add row names and write data table as csv
rownames(sdc_results.df1) <- c("SUPPRESSION", "MANAGED")
sdc_results.df1
write.csv(sdc_results.df, "./sdc_prop_regression.csv")
