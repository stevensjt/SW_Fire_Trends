###########################################################################################################
# 
# This script examines proportion high-severity fire as a function of high-severity perimetere:area 
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
funchir::stale_package_check('PARA_AM_Proportion_HS_Regression.R')

# disable scientific notation
options(scipen=999)
options(digits=6)

# attach data
patch_metrics <- read.csv(file="GEE_1984_2017_final_patch_metrics_709.csv", header=TRUE, sep=",")
names(patch_metrics)

#create empty dataframe to input regression results in
para_results.df <- data.frame(matrix(vector(), 1, 6,
                                    dimnames=list(c(), 
                                                  c("B0", "B1", "P", "RMSE", "R2", "AIC"))),
                             stringsAsFactors=F)


###########################################################################################################
# Analysis of PARA_AM by Proportion High-Severity - SUPPRESSION
###########################################################################################################
dev.off()
# filter by suppression fires only
filter.df <- filter(patch_metrics, FIRE_TYPE2 == "SUPPRESSION")

# set x and y 
x_para <- filter.df$PLAND
y_para <- filter.df$PARA_AM
para.df <- data.frame(x_para, y_para)

# transform data and put into dataframe
x_logit_para <- car::logit(x_para)
x_log_para <- log(x_para)
y_log_para <- log(y_para)
plot(x_log_para, y_log_para)

# create candidate models 
lm.para <- lm(y_para ~ x_para, data = para.df)
logit.model1 <- lm(y_para ~ x_para + I(x_para^2), data = para.df)
logit.model2 <- lm(y_log_para ~ x_log_para + I(x_log_para^2), data = para.df)
logit.model3 <- lm(y_log_para ~ x_logit_para, data = para.df)   
logit.model4 <- lm(y_log_para ~ x_logit_para + I(x_logit_para^2), data = para.df)

# mean centered
x_logit_center <- scale(x_logit_para, center = TRUE, scale = FALSE)
x_logit_center2 <- x_logit_center^2
logit.model5 <- lm(y_log_para ~ x_logit_center + x_logit_center2, data = para.df) #####******
summary(logit.model5)

# test models
AICctab(logit.model4, logit.model5, base = T, delta = T, weights = T)

# normality test 
shapiro.test(resid(lm.para))
shapiro.test(resid(logit.model1))
shapiro.test(resid(logit.model2))
shapiro.test(resid(logit.model3))
shapiro.test(resid(logit.model4))
shapiro.test(resid(logit.model5))

# histogram
fit_para <- studres(logit.model5)
hist(fit_para)

# plot model diagnostics
par(mfrow = c(2, 2))  # Split the plotting panel into a 2 x 2 grid
plot(lm.para)
plot(logit.model1)
plot(logit.model2)
plot(logit.model3)
plot(logit.model4)
plot(logit.model5)

# all diagnostics
gvlma(lm.para)
gvlma(logit.model3)
gvlma(logit.model4)
gvlma(logit.model5)

# RMSE/SUMMARY/AIC
summary(logit.model5)
rmse(logit.model5)
AIC(logit.model5)

# retrieve summary statistics and put into dataframe
lm.para.coef1 <- logit.model5$coefficients[1]
lm.para.coef2 <- logit.model5$coefficients[2]
lm.para.p <- summary(logit.model5)$coefficients[2, 4]
lm.para.rmse <- rmse(logit.model5)
lm.para.adr <- summary(logit.model5)$r.squared
lm.para.aic <- AIC(logit.model5)

para_results.df <- rbind.data.frame(para_results.df, "SUPP" = c(lm.para.coef1, lm.para.coef2, lm.para.p, 
                                                              lm.para.rmse, lm.para.adr, lm.para.aic)) 
para_results.df <- para_results.df[-c(1),]; row.names(para_results.df) <- 1:nrow(para_results.df)
para_results.df

###########################################################################################################
# Analysis of PARA_AM by Proportion High-Severity - Managed Fire
###########################################################################################################
dev.off()
# filter by managed fires only
filter.df1 <- filter(patch_metrics, FIRE_TYPE2 == "OTHER")

# set x and y 
x_para1 <- filter.df1$PLAND
y_para1 <- filter.df1$PARA_AM

# transform data and put into dataframe
x_logit_para1 <- car::logit(x_para1)
x_log_para1 <- log(x_para1)
y_log_para1 <- log(y_para1)
para.df1 <- data.frame(x_para1,y_log_para1)
plot(x_logit_para1, y_log_para1)

# create candidate models 
lm.para1 <- lm(y_para1 ~ x_para1, data = para.df1)
logit.model6 <- lm(y_para1 ~ x_para1 + I(x_para1^2), data = para.df1)
logit.model7 <- lm(y_log_para1 ~ x_log_para1 + I(x_log_para1^2), data = para.df1)
logit.model8 <- lm(y_log_para1 ~ x_logit_para1, data = para.df1)   
logit.model9 <- lm(y_log_para1 ~ x_logit_para1 + I(x_logit_para1^2), data = para.df1)

# mean centered
x_logit_center1 <- scale(x_logit_para1, center = TRUE, scale = FALSE)
x_logit_center12 <- x_logit_center1^2
logit.model10 <- lm(y_log_para1 ~ x_logit_center1 + x_logit_center12, data = para.df1) #####******
summary(logit.model10)

# test models
AICctab(logit.model8, logit.model10, base = T, delta = T, weights = T)
AICctab(logit.model9, logit.model10, base = T, delta = T, weights = T)

# normality test 
shapiro.test(resid(lm.para1))
shapiro.test(resid(logit.model6)) 
shapiro.test(resid(logit.model7))
shapiro.test(resid(logit.model8))
shapiro.test(resid(logit.model9))

# histogram
fit_para1 <- studres(logit.model9)
hist(fit_para1)

# plot model diagnostics
par(mfrow = c(2, 2))  # Split the plotting panel into a 2 x 2 grid
plot(lm.para1)
plot(logit.model7)
plot(logit.model8)
plot(logit.model9)

# all diagnositics
gvlma(lm.para1)
gvlma(logit.model7)
gvlma(logit.model8)
gvlma(logit.model9)

# RMSE/SUMMARY/AIC
summary(logit.model10)
rmse(logit.model10)
AIC(logit.model10)

# retrieve summary statistics and put into dataframe
lm.para.coef1.1 <- logit.model10$coefficients[1]
lm.para.coef2.1 <- logit.model10$coefficients[2]
lm.para.p1 <- summary(logit.model10)$coefficients[2, 4]
lm.para.rmse1 <- rmse(logit.model10)
lm.para.adr1 <- summary(logit.model10)$r.squared
lm.para.aic1 <- AIC(logit.model10)

para_results.df1 <- rbind.data.frame(para_results.df, "MAN" = c(lm.para.coef1.1, lm.para.coef2.1, lm.para.p1, 
                                                              lm.para.rmse1, lm.para.adr1, lm.para.aic1)) 
para_results.df1

# add row names and write data table as csv
rownames(para_results.df1) <- c("SUPPRESSION", "MANAGED")
para_results.df1
write.csv(para_results.df, "./para_prop_regression.csv")
