###########################################################################################################
# 
# This script examines fire size as a function of high-severity class area
# Author: Megan Singleton 
# Updated: 11/24/2020
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
library(bbmle)
library(extrafont)
loadfonts()

# set working directory 
setwd("D:/Patch_Metrics/Patch_Metrics_Complete_Analysis/R/Fire_Size_Regression/")

# check to see if there are any unused packages 
funchir::stale_package_check('CA_Fire_Size_Regression.R')

# disable scientific notation
options(scipen=999)

# attach data
patch_metrics <- read.csv(file="GEE_1984_2017_final_patch_metrics_709.csv", header=TRUE, sep=",")
names(patch_metrics)

#create empty dataframe to input regression results in
ca_results.df <- data.frame(matrix(vector(), 1, 6,
                                   dimnames=list(c(), 
                                                 c("B0", "B1", "P", "RMSE", "R2", "AIC"))),
                            stringsAsFactors=F)

###########################################################################################################
# Analysis of CA by Fire Size - SUPPRESSION
###########################################################################################################

# filter by suppression fires only
filter.df <- filter(patch_metrics, FIRE_TYPE2 == "SUPPRESSION")
names(filter.df)

# set x and y 
x_ca <- filter.df$FORESTED_BURNED
y_ca <- filter.df$CA

# log transform
x_ca_log <- log(x_ca)
y_ca_log <- log(y_ca)
plot(x_ca_log, y_ca_log)

# create new dataframe and plot to get general idea of data
ca.df <- data.frame(x_ca,y_ca)
plot(x_ca, y_ca)

# create linear model
lm.ca <- lm(y_ca ~ x_ca, data = ca.df)
summary(lm.ca)

# test for appropriate transformation
# lambda close to 0 suggest log transformation
boxcox(lm.ca)

# create candidate models
lm.ca_log_x <- lm(y_ca ~ log(x_ca), data = ca.df)
lm.ca_log_y <- lm(log(y_ca) ~ x_ca, data = ca.df)
lm.ca_log_log <- lm(log(y_ca) ~ log(x_ca), data = ca.df)
lm.ca_sqrt_log <- lm(sqrt(y_ca) ~ log(x_ca), data = ca.df)
x2 <- x_ca^2
lm.ca_quad <- lm(y_ca ~ x_ca + x2, data = ca.df)
lm.ca_log_quad <- lm(log(y_ca) ~ x_ca + x2, data = ca.df)
lm.ca_log_log_quad <- lm(log(y_ca) ~ log(x_ca) + log(x2), data = ca.df)

# test models
AICctab(lm.ca_log_log, lm.ca_log_log_quad, base = T, delta = T, weights = T)
AICctab(lm.ca_log_log, lm.ca_log_y, base = T, delta = T, weights = T)

# plot model diagnostics
par(mfrow = c(2, 2))  # Split the plotting panel into a 2 x 2 grid
plot(lm.ca) 
plot(lm.ca_log_y)
plot(lm.ca_log_log)

# normality test 
shapiro.test(resid(lm.ca))
shapiro.test(resid(lm.ca_log_x))
shapiro.test(resid(lm.ca_log_log))
hist(fit_ca)

#testing homoscedasticity
ncvTest(lm.ca_log_x)
ncvTest(lm.ca_log_log)
ncvTest(lm.ca_sqrt_log)

# all diagnostics
gvlma(lm.ca_log_log)
gvlma(lm.ca_log_x)
gvlma(lm.ca)

# retrieve summary statistics and put into dataframe
lm.ca.coef1 <- lm.ca_log_log$coefficients[1]
lm.ca.coef2 <- lm.ca_log_log$coefficients[2]
lm.ca.p <- summary(lm.ca_log_log)$coefficients[2, 4]
lm.ca.rmse <- rmse(lm.ca_log_log)
lm.ca.adr <- summary(lm.ca_log_log)$r.squared
lm.ca.aic <- AIC(lm.ca_log_log)

ca_results.df <- rbind.data.frame(ca_results.df, "SUPP" = c(lm.ca.coef1, lm.ca.coef2, lm.ca.p, 
                                                            lm.ca.rmse, lm.ca.adr, lm.ca.aic)) 
ca_results.df <- ca_results.df[-c(1),]; row.names(ca_results.df) <- 1:nrow(ca_results.df)
ca_results.df

###########################################################################################################
# Analysis of CA by Fire Size - MANAGED
###########################################################################################################

# filter by managed fires only
filter.df1 <- filter(patch_metrics, FIRE_TYPE2 == "OTHER")
names(filter.df1)

# set x and y 
x_ca1 <- filter.df1$FORESTED_BURNED
y_ca1 <- filter.df1$CA

# create new dataframe and plot to get general idea of data
ca.df1 <- data.frame(x_ca1,y_ca1)
dev.off()
plot(ca.df1)

# create linear model
lm.ca1 <- lm(y_ca1 ~ x_ca1, data = ca.df1)
summary(lm.ca1)

# test for appropriate transformation
# lambda close to 0 suggest log transformation
boxcox(lm.ca1)

# create candidate models
lm.ca_log_x1 <- lm(y_ca1 ~ log(x_ca1), data = ca.df1)
lm.ca_log_y1 <- lm(log(y_ca1) ~ x_ca1, data = ca.df1)
lm.ca_log_log1 <- lm(log(y_ca1) ~ log(x_ca1), data = ca.df1)

# plot model diagnostics
par(mfrow = c(2, 2))  # Split the plotting panel into a 2 x 2 grid
plot(lm.ca1)
plot(lm.ca_log_x1)
plot(lm.ca_log_log1)

# normality test 
shapiro.test(resid(lm.ca1))
shapiro.test(resid(lm.ca_log_x1))
shapiro.test(resid(lm.ca_log_log1))
fit_ca1 <- studres(lm.ca_log_log1)
hist(fit_ca1)

# testing homoscedasticity
ncvTest(lm.ca_log_x1)
ncvTest(lm.ca_log_log1)

# all diagnositics
gvlma(lm.ca_log_log1)

#retrieve summary statistics and put into dataframe
lm.ca.f1 <- summary(lm.ca_log_log1)$fstatistic[1]
lm.ca.df1 <- summary(lm.ca_log_log1)$fstatistic[3]
lm.ca.p1 <- summary(lm.ca_log_log1)$coefficients[2, 4]
lm.ca.adr1 <- summary(lm.ca_log_log1)$r.squared

ca_results.df <- rbind.data.frame(ca_results.df, "MAN" = c(lm.ca.f1, lm.ca.df1,lm.ca.p1, lm.ca.adr1)) 
ca_results.df

#RMSE
summary(lm.ca_log_log1)
rmse(lm.ca_log_log1)
AIC(lm.ca_log_log1)

# retrieve summary statistics and put into dataframe
lm.ca.coef1.1 <- lm.ca_log_log1$coefficients[1]
lm.ca.coef2.1 <- lm.ca_log_log1$coefficients[2]
lm.ca.p1 <- summary(lm.ca_log_log1)$coefficients[2, 4]
lm.ca.rmse1 <- rmse(lm.ca_log_log1)
lm.ca.adr1 <- summary(lm.ca_log_log1)$r.squared
lm.ca.aic1 <- AIC(lm.ca_log_log1)

ca_results.df1 <- rbind.data.frame(ca_results.df, "MAN" = c(lm.ca.coef1.1, lm.ca.coef2.1, lm.ca.p1, 
                                                            lm.ca.rmse1, lm.ca.adr1, lm.ca.aic1)) 
ca_results.df1 <- ca_results.df[-c(1),]; row.names(ca_results.df) <- 1:nrow(ca_results.df)
ca_results.df1


#add row names and write data table as csv
rownames(ca_results.df) <- c("SUPPRESSION", "MANAGED")
ca_results.df
write.csv(ca_results.df, "./ca_regression_results_sup_man.csv")

###########################################################################################################
###########################################################################################################
