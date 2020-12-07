###########################################################################################################
# 
# This script examines proportion high-severity fire as a function of high-severity number of disjunnct
# core areas (> 200 m from patch edge)
# Author: Megan Singleton 
# Updated: 12/2/2020
# Notes: quadratic term mean centered 
#
###########################################################################################################

# load libraries 
library(ggplot2)
library(dplyr)
library(MASS)
library(bbmle)
library(sjstats)
library(countreg)
library(pscl) #zero-inflated

# set working directory 
setwd("D:/Patch_Metrics/Patch_Metrics_Complete_Analysis/R/Proportion_Regression/")

# check to see if there are any unused packages 
funchir::stale_package_check('NDCA_Proportion_HS_Regression.R')

# disable scientific notation
options(scipen=999)
options(digits=6)

# attach data
patch_metrics <- read.csv(file="GEE_1984_2017_final_patch_metrics_709.csv", header=TRUE, sep=",")
names(patch_metrics)

###########################################################################################################
# Analysis of NDCA by Proportion High-Severity - SUPPRESSION
###########################################################################################################
dev.off()
# filter by suppression fires only
filter.df <- filter(patch_metrics, FIRE_TYPE2 == "SUPPRESSION")

# set x and y 
x_ndca <- filter.df$PLAND
y_ndca <- filter.df$NDCA
y_log_ndca <- log(y_ndca + 1)
ndca.df <- data.frame(x_ndca,y_ndca)

# plot to test
plot(x_ndca, y_ndca)
plot(x_ndca, y_log_ndca)

# create models 
lm.ndca <- lm(y_ndca ~ x_ndca, data = ndca.df)
lm.ndca_log <- lm(y_log_ndca ~ x_ndca, data = ndca.df)
m1 <- glm(y_ndca ~ x_ndca, data=ndca.df, family = poisson)
m2 <- glm.nb(y_ndca ~ x_ndca, data=ndca.df, control=glm.control(maxit=100))
m3 <- zeroinfl(y_ndca ~ x_ndca | x_ndca,
               data = ndca.df, dist = "negbin")
# quadratic term
x_ndca_c <- scale(x_ndca, center = TRUE, scale = FALSE)
x_ndca_2 <- x_ndca_c^2
m4 <- zeroinfl(y_ndca ~ x_ndca + x_ndca_2 | x_ndca,
               data = ndca.df, dist = "negbin") #####******
summary(m4)

# rootogram for choosing best count model
root.pois <- countreg::rootogram(m1, style = "hanging", plot = FALSE)
root.nb   <- countreg::rootogram(m2, style = "hanging", plot = FALSE)
root.zinb <- countreg::rootogram(m3, style = "hanging", plot = FALSE)
ylims <- ylim(-5, 25)  # common scale for comparison
plot_grid(autoplot(root.pois) + ylims,
          autoplot(root.nb) + ylims, 
          autoplot(root.zinb) + ylims, ncol = 3, labels = "auto")

# test models with vuong
vuong(m2,m3)
vuong(m3,m4)

# select models with AIC
AICctab(m3, m4, base = T, delta = T, weights = T)

# histogram
fit_ndca <- studres(lm.ndca)
hist(fit_ndca)

# SUMMARY/RMSE/AIC/R2
summary(m4)
rmse(m4)
AIC(m4)
rcompanion::nagelkerke(m4)

# likelihood ratio test
lrtest(m4)


###########################################################################################################
# Analysis of NDCA by Proportion High-Severity - MANAGED
###########################################################################################################

dev.off()
# filter by managed fires only
filter.df1 <- filter(patch_metrics, FIRE_TYPE2 == "OTHER")

# set x and y
x_ndca1 <- filter.df1$PLAND
y_ndca1 <- filter.df1$NDCA
y_log_ndca1 <- log(y_ndca1 + 1)
ndca.df1 <- data.frame(x_ndca1,y_ndca1)
plot(x_ndca1, y_ndca1)

# create models 
lm.ndca1 <- lm(y_ndca1 ~ x_ndca1, data = ndca.df1)
lm.ndca_log1 <- lm(y_log_ndca1 ~ x_ndca1, data = ndca.df1)
m5 <- glm(y_ndca1 ~ x_ndca1, data=ndca.df1, family = poisson)
m6 <- glm.nb(y_ndca1 ~ x_ndca1, data=ndca.df1, control=glm.control(maxit=100))
m7 <- zeroinfl(y_ndca1 ~ x_ndca1 | x_ndca1,
               data = ndca.df1, dist = "negbin") #####******
# quadratic term
x_ndca1_2 <- x_ndca1^2
m8 <- zeroinfl(y_ndca1 ~ x_ndca1 + x_ndca1_2 | x_ndca1,
               data = ndca.df1, dist = "negbin") #####******

summary(m8)

#test models
vuong(m7,m8)

# select models with AIC
AICctab(m7, m8, base = T, delta = T, weights = T)

# SUMMARY/RMSE/AIC/R2
summary(m7)
rmse(m7)
AIC(m7)
rcompanion::nagelkerke(m7)

# likelihood ratio test
lrtest(m7)


