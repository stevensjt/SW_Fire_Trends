###########################################################################################################
# 
# This script examines fire size as a function of high-severity patch radius of gyration
# Author: Megan Singleton 
# Updated: 11/30/2020
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
funchir::stale_package_check('GYRATE_AM_Fire_Size_Regression.R')

# disable scientific notation
options(scipen=999)

# attach data
patch_metrics <- read.csv(file="GEE_1984_2017_final_patch_metrics_709.csv", header=TRUE, sep=",")
names(patch_metrics)

#create empty dataframe to input regression results in
gyrate_results.df <- data.frame(matrix(vector(), 1, 6,
                                   dimnames=list(c(), 
                                                 c("B0", "B1", "P", "RMSE", "R2", "AIC"))),
                            stringsAsFactors=F)

###########################################################################################################
# Analysis of GYRATE_AM by Fire Size - SUPPRESSION
###########################################################################################################

dev.off()
# filter by suppression fires only
filter.df <- filter(patch_metrics, FIRE_TYPE2 == "SUPPRESSION")
names(filter.df)

# set x and y 
x_gyrate <- filter.df$FORESTED_BURNED
y_gyrate <- filter.df$GYRATE_AM
plot(x_gyrate, y_gyrate)

# log transform
x_gyrate_log <- log(x_gyrate)
y_gyrate_log <- log(y_gyrate)
plot(x_gyrate_log, y_gyrate_log)

# create new dataframe to create models
gyrate.df <- data.frame(x_gyrate,y_gyrate)

# create linear model
lm.gyrate <- lm(y_gyrate ~ x_gyrate, data = gyrate.df)
summary(lm.gyrate)

# test for appropriate transformation
# lambda close to 0 suggest log transformation
boxcox(lm.gyrate)

# create candidate models
lm.gyrate_log_x <- lm(y_gyrate ~ log(x_gyrate), data = gyrate.df)
lm.gyrate_log_y <- lm(log(y_gyrate) ~ x_gyrate, data = gyrate.df)
lm.gyrate_log_log <- lm(log(y_gyrate) ~ log(x_gyrate), data = gyrate.df)
lm.gyrate_sqrt_log <- lm(sqrt(y_gyrate) ~ log(x_gyrate), data = gyrate.df)
x2 <- x_gyrate^2
lm.gyrate_quad <- lm(y_gyrate ~ x_gyrate + x2, data = gyrate.df)
lm.gyrate_log_quad <- lm(log(y_gyrate) ~ x_gyrate + x2, data = gyrate.df)
lm.gyrate_log_log_quad <- lm(log(y_gyrate) ~ log(x_gyrate) + log(x2), data = gyrate.df)

# test models
AICctab(lm.gyrate_log_log, lm.gyrate_log_log_quad, base = T, delta = T, weights = T)
AICctab(lm.gyrate_log_log, lm.gyrate_log_y, base = T, delta = T, weights = T)

# plot model diagnostics
par(mfrow = c(2, 2))  # Split the plotting panel into a 2 x 2 grid
plot(lm.gyrate)
plot(lm.gyrate_log_y)
plot(lm.gyrate_log_log)

#normality test 
shapiro.test(resid(lm.gyrate))
shapiro.test(resid(lm.gyrate_log_x))
shapiro.test(resid(lm.gyrate_log_log))
fit_gyrate <- studres(lm.gyrate_log_log)
hist(fit_gyrate)

#testing homoscedasticity
ncvTest(lm.gyrate_log_x)
ncvTest(lm.gyrate_log_log)

# all diagnostics
gvlma(lm.gyrate)
gvlma(lm.gyrate_log_y)
gvlma(lm.gyrate_log_log)

# retrieve summary statistics and put into dataframe
lm.gyrate.coef1 <- lm.gyrate_log_log$coefficients[1]
lm.gyrate.coef2 <- lm.gyrate_log_log$coefficients[2]
lm.gyrate.p <- summary(lm.gyrate_log_log)$coefficients[2, 4]
lm.gyrate.rmse <- rmse(lm.gyrate_log_log)
lm.gyrate.adr <- summary(lm.gyrate_log_log)$r.squared
lm.gyrate.aic <- AIC(lm.gyrate_log_log)

gyrate_results.df <- rbind.data.frame(gyrate_results.df, "SUPP" = c(lm.gyrate.coef1, lm.gyrate.coef2, lm.gyrate.p, 
                                                            lm.gyrate.rmse, lm.gyrate.adr, lm.gyrate.aic)) 
gyrate_results.df <- gyrate_results.df[-c(1),]; row.names(gyrate_results.df) <- 1:nrow(gyrate_results.df)
gyrate_results.df

###########################################################################################################
# Analysis of GYRATE_AM by Fire Size - MANAGED
###########################################################################################################

# filter by managed fires only
filter.df1 <- filter(patch_metrics, FIRE_TYPE2 == "OTHER")
names(filter.df1)

# set x and y 
x_gyrate1 <- filter.df1$FORESTED_BURNED
y_gyrate1 <- filter.df1$GYRATE_AM

# create new dataframe and plot to get general idea of data
gyrate.df1 <- data.frame(x_gyrate1,y_gyrate1)
dev.off()
plot(gyrate.df1)

# create linear model
lm.gyrate1 <- lm(y_gyrate1 ~ x_gyrate1, data = gyrate.df1)
summary(lm.gyrate1)

# test for appropriate transformation
# lambda close to 0 suggest log transformation
boxcox(lm.gyrate1)

# create candidate models
lm.gyrate_log_x1 <- lm(y_gyrate1 ~ log(x_gyrate1), data = gyrate.df1)
lm.gyrate_log_y1 <- lm(log(y_gyrate1) ~ x_gyrate1, data = gyrate.df1)
lm.gyrate_log_log1 <- lm(log(y_gyrate1) ~ log(x_gyrate1), data = gyrate.df1)
lm.gyrate_sqrt1 <- lm((1/(sqrt(y_gyrate1))) ~ x_gyrate1, data = gyrate.df1)

# plot model diagnostics
par(mfrow = c(2, 2))  # Split the plotting panel into a 2 x 2 grid
plot(lm.gyrate1)
plot(lm.gyrate_log_x1)
plot(lm.gyrate_log_log1)
plot(lm.gyrate_sqrt1)

# normality test 
shapiro.test(resid(lm.gyrate1))
shapiro.test(resid(lm.gyrate_log_x1))
shapiro.test(resid(lm.gyrate_log_log1))
fit_gyrate1 <- studres(lm.gyrate_log_log1)
hist(fit_gyrate1)

# testing homoscedasticity
ncvTest(lm.gyrate_log_x1)
ncvTest(lm.gyrate_log_log1)

# all diagnositics
gvlma(lm.gyrate_log_log1)

# retrieve summary statistics and put into dataframe
lm.gyrate.coef1.1 <- lm.gyrate_log_log1$coefficients[1]
lm.gyrate.coef2.1 <- lm.gyrate_log_log1$coefficients[2]
lm.gyrate.p1 <- summary(lm.gyrate_log_log1)$coefficients[2, 4]
lm.gyrate.rmse1 <- rmse(lm.gyrate_log_log1)
lm.gyrate.adr1 <- summary(lm.gyrate_log_log1)$r.squared
lm.gyrate.aic1 <- AIC(lm.gyrate_log_log1)

gyrate_results.df1 <- rbind.data.frame(gyrate_results.df, "MAN" = c(lm.gyrate.coef1.1, lm.gyrate.coef2.1, lm.gyrate.p1, 
                                                            lm.gyrate.rmse1, lm.gyrate.adr1, lm.gyrate.aic1)) 

#add row names and write data table as csv
rownames(gyrate_results.df1) <- c("SUPPRESSION", "MANAGED")
gyrate_results.df1
write.csv(gyrate_results.df, "./gyrate_regression_results.csv")

###########################################################################################################
###########################################################################################################
