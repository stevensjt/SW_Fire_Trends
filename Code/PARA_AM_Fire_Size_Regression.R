###########################################################################################################
# 
# This script examines fire size as a function of high-severity area-weighted perimeter-area ratio
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
funchir::stale_package_check('PARA_AM_Fire_Size_Regression.R')

# disable scientific notation
options(scipen=999)

# attach data
patch_metrics <- read.csv(file="GEE_1984_2017_final_patch_metrics_709.csv", header=TRUE, sep=",")
names(patch_metrics)

#create empty dataframe to input regression results in
para_results.df <- data.frame(matrix(vector(), 1, 6,
                                       dimnames=list(c(), 
                                                     c("B0", "B1", "P", "RMSE", "R2", "AIC"))),
                                stringsAsFactors=F)

###########################################################################################################
# Analysis of PARA_AM by Fire Size - SUPPRESSION
###########################################################################################################
dev.off()
# filter by suppression fires only
filter.df <- filter(patch_metrics, FIRE_TYPE2 == "SUPPRESSION")
names(filter.df)

# set x and y 
x_para <- filter.df$FORESTED_BURNED
y_para <- filter.df$PARA_AM
para.df <- data.frame(x_para,y_para)
plot(x_para, y_para)

# log transform
x_para_log <- log(x_para)
y_para_log <- log(y_para)
plot(x_para_log, y_para_log)

# create linear model
lm.para <- lm(y_para ~ x_para, data = para.df)
summary(lm.para)

# test for appropriate transformation
car::boxCox(lm.para)

# create candidate models
lm.para_log_x <- lm(y_para ~ log(x_para), data = para.df)
lm.para_log_y <- lm(log(y_para) ~ x_para, data = para.df)
lm.para_log_log <- lm(log(y_para) ~ log(x_para), data = para.df)
lm.para_sqrt_log <- lm(sqrt(y_para) ~ log(x_para), data = para.df)
lm.para_sqrt <- lm(sqrt(y_para) ~ sqrt(x_para), data = para.df)
x2 <- x_para^2
lm.para_quad <- lm(y_para ~ x_para + x2, data = para.df)
lm.para_log_quad <- lm(log(y_para) ~ x_para + x2, data = para.df)
lm.para_log_log_quad <- lm(log(y_para) ~ log(x_para) + log(x2), data = para.df)
lm.para_sqrt_log <- lm(sqrt(y_para) ~ log(x_para), data = para.df)

# test models
AICctab(lm.para_log_log, lm.para_log_log_quad, base = T, delta = T, weights = T)

# plot model diagnostics
par(mfrow = c(2, 2))  # Split the plotting panel into a 2 x 2 grid
plot(lm.para)
plot(lm.para_log_y)
plot(lm.para_log_log)
plot(lm.para_sqrt_log)

#normality test 
shapiro.test(resid(lm.para))
shapiro.test(resid(lm.para_log_x))
shapiro.test(resid(lm.para_log_log))
shapiro.test(resid(lm.para_sqrt_log))
fit_para <- studres(lm.para_log_log)
fit_para <- studres(lm.para_sqrt_log)
hist(fit_para)

#testing homoscedasticity
ncvTest(lm.para_log_x)
ncvTest(lm.para_log_log)

# all diagnostics
gvlma(lm.para)
gvlma(lm.para_log_y)
gvlma(lm.para_log_log)
gvlma(lm.para_sqrt_log)

# RMSE/SUMMARY/AIC
summary(lm.para_log_log)
rmse(lm.para_log_log)
AIC(lm.para_log_log)

# retrieve summary statistics and put into dataframe
lm.para.coef1 <- lm.para_log_log$coefficients[1]
lm.para.coef2 <- lm.para_log_log$coefficients[2]
lm.para.p <- summary(lm.para_log_log)$coefficients[2, 4]
lm.para.rmse <- rmse(lm.para_log_log)
lm.para.adr <- summary(lm.para_log_log)$r.squared
lm.para.aic <- AIC(lm.para_log_log)

para_results.df <- rbind.data.frame(para_results.df, "SUPP" = c(lm.para.coef1, lm.para.coef2, lm.para.p, 
                                                                    lm.para.rmse, lm.para.adr, lm.para.aic)) 
para_results.df <- para_results.df[-c(1),]; row.names(para_results.df) <- 1:nrow(para_results.df)
para_results.df
summary(lm.para_log_log)
summary(lm.para_sqrt_log)

###########################################################################################################
# Analysis of PARA_AM by Fire Size - MANAGED
###########################################################################################################

# filter by managed fires only
filter.df1 <- filter(patch_metrics, FIRE_TYPE2 == "OTHER")
names(filter.df1)

# set x and y 
x_para1 <- filter.df1$FORESTED_BURNED
y_para1 <- filter.df1$PARA_AM

# create new dataframe and plot to get general idea of data
para.df1 <- data.frame(x_para1,y_para1)
dev.off()
plot(para.df1)

# create linear model
lm.para1 <- lm(y_para1 ~ x_para1, data = para.df1)
summary(lm.para1)

# test for appropriate transformation
# lambda close to 0 suggest log transformation
car::boxCox(lm.para1)

# create candidate models
lm.para_log_x1 <- lm(y_para1 ~ log(x_para1), data = para.df1)
lm.para_log_y1 <- lm(log(y_para1) ~ x_para1, data = para.df1)
lm.para_log_log1 <- lm(log(y_para1) ~ log(x_para1), data = para.df1)
lm.para_sqrt1 <- lm(sqrt(y_para1) ~ x_para1, data = para.df1)
lm.para_sqrt_sqrt1 <- lm(sqrt(y_para1) ~ sqrt(x_para1), data = para.df1)
lm.para_log_sqrt1 <- lm(log(y_para1) ~ sqrt(x_para1), data = para.df1)
lm.para_sqrt_log1 <- lm(sqrt(y_para1) ~ log(x_para1), data = para.df1)

# plot model diagnostics
par(mfrow = c(2, 2))  # Split the plotting panel into a 2 x 2 grid
plot(lm.para1)
plot(lm.para_log_x1)
plot(lm.para_log_log1)
plot(lm.para_sqrt1)
plot(lm.para_sqrt_sqrt1)
plot(lm.para_log_sqrt1)
plot(lm.para_sqrt_log1)

# normality test 
shapiro.test(resid(lm.para1))
shapiro.test(resid(lm.para_log_x1))
shapiro.test(resid(lm.para_log_log1))
shapiro.test(resid(lm.para_sqrt_sqrt1))
shapiro.test(resid(lm.para_sqrt1))
shapiro.test(resid(lm.para_sqrt_log1))
fit_para1 <- studres(lm.para_log_log1)
fit_para1 <- studres(lm.para_sqrt1)
hist(fit_para1)

# testing homoscedasticity
ncvTest(lm.para_log_x1)
ncvTest(lm.para_log_log1)

# all diagnositics
gvlma(lm.para_log_log1)
gvlma(lm.para_sqrt_log1)

# RMSE/SUMMARY/AIC
summary(lm.para_log_log1)
rmse(lm.para_log_log1)
AIC(lm.para_log_log1)

# retrieve summary statistics and put into dataframe
lm.para.coef1.1 <- lm.para_log_log1$coefficients[1]
lm.para.coef2.1 <- lm.para_log_log1$coefficients[2]
lm.para.p1 <- summary(lm.para_log_log1)$coefficients[2, 4]
lm.para.rmse1 <- rmse(lm.para_log_log1)
lm.para.adr1 <- summary(lm.para_log_log1)$r.squared
lm.para.aic1 <- AIC(lm.para_log_log1)

para_results.df1 <- rbind.data.frame(para_results.df, "MAN" = c(lm.para.coef1.1, lm.para.coef2.1, lm.para.p1, 
                                                                    lm.para.rmse1, lm.para.adr1, lm.para.aic1))
summary(lm.para_sqrt_log1)
summary(lm.para_log_log1)

#add row names and write data table as csv
rownames(para_results.df1) <- c("SUPPRESSION", "MANAGED")
para_results.df1
write.csv(para_results.df, "./para_regression_results_sup_man.csv")

###########################################################################################################
###########################################################################################################
