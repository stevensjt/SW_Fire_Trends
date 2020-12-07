###########################################################################################################
# 
# This script plots 6 different patch metrics as a function of time.
# Including: class area, patch radius of gyration, perimeter area ratio, number of core areas, percentage of 
# like adjacencies and the stand decay coefficient
# Author: Megan Singleton 
# Updated: 12/3/2020
#
###########################################################################################################

# load libraries
library(ggplot2)
library(dplyr)
library(cowplot)
library(bbmle)
library(betareg)
library(extrafont)
library(pscl)
font_import()
loadfonts()
loadfonts(device = "win")

# set working directory 
setwd("D:/Patch_Metrics/Patch_Metrics_Complete_Analysis/R/Time_Series/")

# check to see if there are any unused packages 
funchir::stale_package_check('ggplot_all_6_Time_Series.R')

# disable scientific notation
options(scipen=999)
options(digits = 6)

# attach data
patch_metrics <- read.csv(file="GEE_1984_2017_final_patch_metrics_709.csv", header=TRUE, sep=",")
names(patch_metrics)

# filter by suppression fires only
filter.df <- filter(patch_metrics, FIRE_TYPE2 == "SUPPRESSION")

# filter by managed fires only
filter.df1 <- filter(patch_metrics, FIRE_TYPE2 == "OTHER")

p <- theme(plot.title = element_text(family = "Calibri", color = "black", size =13, hjust = 0.5, lineheight = 0.8, face="bold"),
      #axis.text.x = element_text(family = "Calibri", color="black", size=8),
      axis.title.y = element_text(family = "Calibri", color="black", size=12),
      #axis.title.y = element_blank(),
      plot.margin = unit(c(0.1, 0.25, 0.1, 0.25), "cm"), #top #left #bottom #right
      panel.border = element_rect(colour = "black", fill=NA, size=0.5))

###########################################################################################################
# Plot Class Area (CA)
###########################################################################################################

# read in data and remove years with 1 fire 
mean_ca_sup <- read.csv("./mean_ca_sup.csv", header=TRUE, sep=",")
mean_ca_man <- read.csv("./mean_ca_man.csv", header=TRUE, sep=",")
mean_ca_sup <- mean_ca_sup[-c(1),]
mean_ca_man <- mean_ca_man[-c(1),]

# suppression fires
x_ca <- mean_ca_sup$YEAR
y_ca <- mean_ca_sup$CA
#y_log_ca <- log(y_ca)
ca.df <- data.frame(x_ca,y_ca)
myts.ca <- ts(mean_ca_sup$CA, start=1985, end= 2017)
#fit <- gls(myts.log ~ x_ca, correlation = corARMA(p=1), method = "ML") #******
fit <- gls(myts.ca ~ x_ca, correlation = corARMA(p=1), method = "ML") #******
ca_fitted.df <- fitted(fit)
ca.df <- data.frame(x = 1:33, fitted = ca_fitted.df)
ca.df <- cbind(ca.df, x_ca, y_ca)
names(ca.df) <- c("x_ca", 'predicted_ca', 'fire_ca', 'observed_ca')

# managed fires
x_ca1 <- mean_ca_man$YEAR
y_ca1 <- mean_ca_man$CA
#y_log_ca1 <- log(y_ca1)
ca.df1 <- data.frame(x_ca1,y_ca1)
plot(ca.df1)
lm.ca1 <- lm(y_ca1 ~ x_ca1, data=ca.df1)
ca_fitted.df1 <- fitted(lm.ca1)
ca.df1 <- data.frame(x = 1:21, fitted = ca_fitted.df1)
ca.df1 <- cbind(ca.df1, x_ca1, y_ca1)
names(ca.df1) <- c("x_ca", 'predicted_ca', 'fire_ca', 'observed_ca')

# plot
ca_combine <- ggplot(data = ca.df, aes(x = fire_ca, y = observed_ca)) + 
  geom_point(data = ca.df, aes(x = fire_ca, y = observed_ca), color = "indianred3", size = 2) + 
  geom_line(data = ca.df, aes(y = predicted_ca), color = "indianred3", size = 1) + 
  geom_line(data = ca.df, aes(x = fire_ca, y = observed_ca), color = "indianred3", size = 1) +
  geom_point(data = ca.df1, aes(x = fire_ca, y = observed_ca), color = "turquoise3", size = 2) + 
  geom_line(data = ca.df1, aes(y = predicted_ca), color = "turquoise3", size = 1) +
  geom_line(data = ca.df1, aes(x = fire_ca, y = observed_ca), color = "turquoise3", size = 1) +
  theme_bw() + 
  ggtitle("\nhigh-severity area\n(CA)*") + 
  labs(x="", y = "ha") +
  #scale_y_continuous(breaks=seq(0,10,2)) +
  scale_x_continuous(breaks=seq(1985,2015,5)) +
  panel_border() 

ca_combine <- ca_combine + p
ca_combine

###########################################################################################################
# Plot Patch Radius of Gyration (GYRATE_AM)
###########################################################################################################

# read in data and remove years with 1 fire 
mean_gyrate_sup <- read.csv("./mean_gyrate_sup.csv", head = TRUE, sep=",")
mean_gyrate_man <- read.csv("./mean_gyrate_man.csv", header=TRUE, sep=",")
mean_gyrate_sup <- mean_gyrate_sup[-c(1),]
mean_gyrate_man <- mean_gyrate_man[-c(1),]

# suppression fires
x_gyrate <- mean_gyrate_sup$YEAR
y_gyrate <- mean_gyrate_sup$GYRATE_AM
time_gyrate <- seq(from = 1, to = 33, by = 1)
time_gyrate <- time_gyrate - mean(time_gyrate)
time2_gyrate <- time_gyrate^2
gyrate.df <- data.frame(x_gyrate,y_gyrate)
plot(gyrate.df)
#x2 <- x_gyrate^2
#lm.gyrate <- lm(y_gyrate ~ x_gyrate + x2)
myts.gyrate <- ts(mean_gyrate_sup$GYRATE_AM, start=1985, end= 2017)
fit_gyrate <- gls(myts.gyrate ~ time_gyrate + time2_gyrate, correlation = corARMA(p=1), method = "ML") #*****
#lm.gyrate <- lm(y_log_gyrate ~ x_gyrate, data = gyrate.df)
summary(fit_gyrate)
gyrate_fitted.df <- fitted(fit_gyrate)
gyrate.df <- data.frame(x = 1:33, fitted = gyrate_fitted.df)
gyrate.df <- cbind(gyrate.df, x_gyrate, y_gyrate)
names(gyrate.df) <- c("x_gyrate", 'predicted_gyrate', 'fire_gyrate', 'observed_gyrate')
names(gyrate.df)

# managed fires
x_gyrate1 <- mean_gyrate_man$YEAR
y_gyrate1 <- mean_gyrate_man$GYRATE_AM
time_gyrate1 <- seq(from = 1, to = 21, by = 1)
time_gyrate1 <- time_gyrate1 - mean(time_gyrate1)
time2_gyrate1 <- time_gyrate1^2
#y_sqrt_gyrate1 <- sqrt(y_gyrate1)
gyrate.df1 <- data.frame(x_gyrate1,y_gyrate1)
plot(gyrate.df1)
#x22 <- x_gyrate1^2
#lm.gyrate1 <- lm(y_gyrate1 ~ x_gyrate1 + x22)
fit_gyrate1 <- lm(y_gyrate1 ~ time_gyrate1 + time2_gyrate1) #*****
#lm.gyrate1 <- lm(y_log_gyrate1 ~ x_gyrate1, data = gyrate.df1)
gyrate_fitted.df1 <- fitted(fit_gyrate1)
gyrate.df1 <- data.frame(x = 1:21, fitted = gyrate_fitted.df1)
gyrate.df1 <- cbind(gyrate.df1, x_gyrate1, y_gyrate1)
names(gyrate.df1) <- c("x_gyrate", 'predicted_gyrate', 'fire_gyrate', 'observed_gyrate')
names(gyrate.df1)

# plot
gyrate_combine <- ggplot(data = gyrate.df, aes(x = fire_gyrate, y = observed_gyrate)) + 
  geom_point(data = gyrate.df, aes(x = fire_gyrate, y = observed_gyrate), color = "indianred3", size = 2) + 
  geom_line(data = gyrate.df, aes(y = predicted_gyrate), color = "indianred3", size = 1) + 
  geom_line(data = gyrate.df, aes(x = fire_gyrate, y = observed_gyrate), color = "indianred3", size = 1) +
  geom_point(data = gyrate.df1, aes(x = fire_gyrate, y = observed_gyrate), color = "turquoise3", size = 2) + 
  geom_line(data = gyrate.df1, aes(y = predicted_gyrate), color = "turquoise3", size = 1) +
  geom_line(data = gyrate.df1, aes(x = fire_gyrate, y = observed_gyrate), color = "turquoise3", size = 1) +
  theme_bw() + 
  ggtitle("\npatch reach\n(GYRATE_AM)") +
  #ggtitle("\npatch radius of gyration\n(GYRATE_AM)") + 
  labs(x="", y = "meters") +
  #scale_y_continuous(breaks=seq(0,10,2)) +
  scale_x_continuous(breaks=seq(1985,2015,5)) +
  panel_border()

gyrate_combine <- gyrate_combine + p
gyrate_combine

###########################################################################################################
# Plot Number of Disjunct Core Areas (NDCA)
###########################################################################################################

# read in data and remove years with 1 fire 
sum_ndca_sup <- read.csv("./sum_ndca_sup.csv", header=TRUE, sep=",")
sum_ndca_man <- read.csv("./sum_ndca_man.csv", header=TRUE, sep=",")
sum_ndca_sup <- sum_ndca_sup[-c(1),]
sum_ndca_man <- sum_ndca_man[-c(1),]

# suppression fires
x_ndca <- sum_ndca_sup$YEAR
y_ndca <- sum_ndca_sup$NDCA
ndca.df <- data.frame(x_ndca,y_ndca)
plot(x_ndca, y_ndca)
m_ndca  <- glm.nb(y_ndca ~ x_ndca, data=ndca.df, control=glm.control(maxit=100))
ndca_fitted.df <- fitted(m_ndca)
ndca.df <- data.frame(x = 1:33, fitted = ndca_fitted.df)
ndca.df <- cbind(ndca.df, x_ndca, y_ndca)
names(ndca.df) <- c("x_ndca", 'predicted_ndca', 'fire_ndca', 'observed_ndca')
names(ndca.df)

# managed fires
x_ndca1 <- sum_ndca_man$YEAR
y_ndca1 <- sum_ndca_man$NDCA
ndca.df1 <- data.frame(x_ndca1,y_ndca1)
plot(x_ndca1, y_ndca1)
m1_ndca  <- glm.nb(y_ndca1 ~ x_ndca1, data=ndca.df1, control=glm.control(maxit=100))
ndca_fitted.df1 <- fitted(m1_ndca)
ndca.df1 <- data.frame(x = 1:21, fitted = ndca_fitted.df1)
ndca.df1 <- cbind(ndca.df1, x_ndca1, y_ndca1)
names(ndca.df1) <- c("x_ndca", 'predicted_ndca', 'fire_ndca', 'observed_ndca')
names(ndca.df1)

# plot
ndca_combine <- ggplot(data = ndca.df, aes(x = fire_ndca, y = observed_ndca)) + 
  geom_point(data = ndca.df, aes(x = fire_ndca, y = observed_ndca), color = "indianred3", size = 2) + 
  geom_line(data = ndca.df, aes(y = predicted_ndca), color = "indianred3", size = 1) + 
  geom_line(data = ndca.df, aes(x = fire_ndca, y = observed_ndca), color = "indianred3", size = 1) +
  geom_point(data = ndca.df1, aes(x = fire_ndca, y = observed_ndca), color = "turquoise3", size = 2) + 
  geom_line(data = ndca.df1, aes(y = predicted_ndca), color = "turquoise3", size = 1) +
  geom_line(data = ndca.df1, aes(x = fire_ndca, y = observed_ndca), color = "turquoise3", size = 1) +
  theme_bw() + 
  #ggtitle("\nnumber of disjunct\ncore areas (NDCA)") + 
  ggtitle("\nnumber of core areas\n(NDCA)") +
  labs(x="", y = "count") +
  #scale_y_continuous(breaks=seq(0,10,2)) +
  scale_y_continuous(breaks =seq(0,550,100)) +
  scale_x_continuous(breaks=seq(1985,2015,5)) +
  panel_border()

ndca_combine <- ndca_combine + p
ndca_combine

###########################################################################################################
# Plot Patch Perimeter:Ratio Ratio (PARA_AM)
###########################################################################################################

# read in data and remove years with 1 fire 
mean_para_sup <- read.csv("./mean_para_sup.csv", header=TRUE, sep=",")
mean_para_man <- read.csv("./mean_para_man.csv", header=TRUE, sep=",")
mean_para_sup <- mean_para_sup[-c(1),]
mean_para_man <- mean_para_man[-c(1),]

# suppression fires
x_para <- mean_para_sup$YEAR
y_para <- mean_para_sup$PARA_AM
time_para <- seq(from = 1, to = 33, by = 1)
time_para <- time_para - mean(time_para)
time2_para <- time_para^2
#x_log_para <- log(x_para)
#y_log_para <- log(y_para)
para.df <- data.frame(x_para,y_para)
plot(para.df)
myts_para <- ts(mean_para_sup$PARA_AM, start=1985, end= 2017)
fit_para <- gls(myts_para ~ time_para + time2_para, correlation = corARMA(q=1), method = "ML")  ##****
para_fitted.df <- fitted(fit_para)
para.df <- data.frame(x = 1:33, fitted = para_fitted.df)
para.df <- cbind(para.df, x_para, y_para)
names(para.df) <- c("x_para", 'predicted_para', 'fire_para', 'observed_para')
names(para.df)

# managed fires
x_para1 <- mean_para_man$YEAR
y_para1 <- mean_para_man$PARA_AM
y_per1 <- y_para1/100
#x_log_para1 <- log(x_para1)
#y_log_para1 <- log(y_para1)
para.df1 <- data.frame(x_para1,y_para1)
plot(para.df1)
lm.para1 <- lm(y_para1 ~ x_para1, data = para.df1)
para_fitted.df1 <- fitted(lm.para1)
para.df1 <- data.frame(x = 1:21, fitted = para_fitted.df1)
para.df1 <- cbind(para.df1, x_para1, y_para1)
names(para.df1) <- c("x_para", 'predicted_para', 'fire_para', 'observed_para')
names(para.df1)

# plot
para_combine <- ggplot(data = para.df, aes(x = fire_para, y = observed_para)) + 
  geom_point(data = para.df, aes(x = fire_para, y = observed_para), color = "indianred3", size = 2) + 
  geom_line(data = para.df, aes(y = predicted_para), color = "indianred3", size = 1) + 
  geom_line(data = para.df, aes(x = fire_para, y = observed_para), color = "indianred3", size = 1) +
  geom_point(data = para.df1, aes(x = fire_para, y = observed_para), color = "turquoise3", size = 2) + 
  geom_line(data = para.df1, aes(y = predicted_para), color = "turquoise3", size = 1) +
  geom_line(data = para.df1, aes(x = fire_para, y = observed_para), color = "turquoise3", size = 1) +
  theme_bw() + 
  ggtitle("\nshape complexity\n(PARA_AM)") +
  #ggtitle("\nperimeter:area\n (PARA_AM)") + 
  labs(x="", y = "perimeter:area") +
  #scale_y_continuous(breaks=seq(0,10,2)) +
  scale_x_continuous(breaks=seq(1985,2015,5)) +
  panel_border() 

para_combine <- para_combine + p
para_combine

###########################################################################################################
# Plot Patch Percentage of Like Adjacencies 
###########################################################################################################

# read in data and remove years with 1 fire 
mean_pladj_sup <- read.csv("./mean_pladj_sup.csv", header=TRUE, sep=",")
mean_pladj_man <- read.csv("./mean_pladj_man.csv", header=TRUE, sep=",")
mean_pladj_sup <- mean_pladj_sup[-c(1),]
mean_pladj_man <- mean_pladj_man[-c(1),]

# suppression fires
x_pladj <- mean_pladj_sup$YEAR
y_pladj <- mean_pladj_sup$PLADJ
time_pladj <- seq(from = 1, to = 33, by = 1)
time_pladj <- time_pladj - mean(time_pladj)
time2_pladj <- time_pladj^2
myts_pladj <- ts((mean_pladj_sup$PLADJ), start= 1985, end=2017)
#y_pladj_per <- y_pladj/100
pladj.df <- data.frame(x_pladj,y_pladj)
#x2 <- x_pladj^2
#lm.pladj <- lm(y_pladj ~ x_pladj + x2)
fit_pladj <- gls(myts_pladj ~ time_pladj + time2_pladj, correlation = corARMA(q=1), method = "ML")  ##**
pladj_fitted.df <- fitted(fit_pladj)
pladj.df <- data.frame(x = 1:33, fitted = pladj_fitted.df)
pladj.df <- cbind(pladj.df, x_pladj, y_pladj)
names(pladj.df) <- c("x_pladj", 'predicted_pladj', 'fire_pladj', 'observed_pladj')
names(pladj.df)

# managed fires
x_pladj1 <- mean_pladj_man$YEAR
y_pladj1 <- mean_pladj_man$PLADJ
yper1 <- y_pladj1/100
#y_pladj_per1 <- y_pladj1/100
#pladj.df1 <- data.frame(x_pladj1,y_pladj1)
#plot(x_pladj1, y_pladj1)
#lm.pladj1 <- lm(y_pladj1 ~ x_pladj1)
pladj_fit1 <- betareg(yper1 ~ x_pladj1, data = mean_pladj_man) #***************
pladj_fitted.df1 <- fitted(pladj_fit1)
pladj.df1 <- data.frame(x = 1:21, fitted = pladj_fitted.df1)
pladj.df1 <- cbind(pladj.df1, x_pladj1, y_pladj1)
names(pladj.df1) <- c("x_pladj", 'predicted_pladj', 'fire_pladj', 'observed_pladj')
names(pladj.df1)
pladj.df1$predicted_pladj <- pladj.df1$predicted_pladj * 100

# plot
pladj_combine <- ggplot(data = pladj.df, aes(x = fire_pladj, y = observed_pladj)) + 
  geom_point(data = pladj.df, aes(x = fire_pladj, y = observed_pladj), color = "indianred3", size = 2) + 
  geom_line(data = pladj.df, aes(y = predicted_pladj), color = "indianred3", size = 1) + 
  geom_line(data = pladj.df, aes(x = fire_pladj, y = observed_pladj), color = "indianred3", size = 1) +
  geom_point(data = pladj.df1, aes(x = fire_pladj, y = observed_pladj), color = "turquoise3", size = 2) + 
  geom_line(data = pladj.df1, aes(y = predicted_pladj), color = "turquoise3", size = 1) +
  geom_line(data = pladj.df1, aes(x = fire_pladj, y = observed_pladj), color = "turquoise3", size = 1) +
  theme_bw() + 
  ggtitle("\naggregation\n(PLADJ)") +
  #ggtitle("\npercentage like\nadjacencies (PLADJ)") + 
  labs(x="", y = "percent") +
  #scale_y_continuous(breaks=seq(0,10,2)) +
  scale_x_continuous(breaks=seq(1985,2015,5)) +
  panel_border()

pladj_combine <- pladj_combine + p
pladj_combine

###########################################################################################################
# Plot High-Severity Stand-Decay Coefficient 
###########################################################################################################

# read in data and remove years with 1 fire 
mean_sdc_sup <- read.csv("./mean_sdc_sup.csv", header=TRUE, sep=",")
mean_sdc_man <- read.csv("./mean_sdc_man.csv", header=TRUE, sep=",")
mean_sdc_sup <- mean_sdc_sup[-c(1,3,8),]
mean_sdc_man <- mean_sdc_man[-c(1:2),]

# suppression fires
x_sdc <- mean_sdc_sup$YEAR
y_sdc <- mean_sdc_sup$SDC
#y_log_sdc <- log(y_sdc)
sdc.df <- data.frame(x_sdc,y_sdc)
plot(x_sdc, y_sdc)
#myts_sdc <- ts(mean_sdc_sup$SDC, start=1985, end= 2017)
#time_sdc <- seq(from = 1, to = 33, by = 1)
#time_sdc <- time_sdc - mean(time_sdc)
x_sdc2 <- x_sdc^2
fit <- lm(y_sdc ~ x_sdc + x_sdc2) #*****
sdc_fitted.df <- fitted(fit)
sdc.df <- data.frame(x = 1:31, fitted = sdc_fitted.df)
sdc.df <- cbind(sdc.df, x_sdc, y_sdc)
names(sdc.df) <- c("x_sdc", 'predicted_sdc', 'fire_sdc', 'observed_sdc')
names(sdc.df)

# managed fires
x_sdc1 <- mean_sdc_man$YEAR
y_sdc1 <- mean_sdc_man$SDC
#time_sdc1 <- seq(from = 1, to = 20, by = 1)
#time_sdc1 <- time_sdc1 - mean(time_sdc1)
x_sdc12 <- x_sdc1^2
sdc.df1 <- data.frame(x_sdc1,y_sdc1)
plot(x_sdc1, y_sdc1)
fit_sdc1 <- lm(y_sdc1 ~ x_sdc1 + x_sdc12, data = sdc.df1) #***
#fit_sdc1 <- lm(ylog ~ time + time2) #***
sdc_fitted.df1 <- fitted(fit_sdc1)
sdc.df1 <- data.frame(x = 1:19, fitted = sdc_fitted.df1)
sdc.df1 <- cbind(sdc.df1, x_sdc1, y_sdc1)
names(sdc.df1) <- c("x_sdc", 'predicted_sdc', 'fire_sdc', 'observed_sdc')
names(sdc.df1)

# plot
sdc_combine <- ggplot(data = sdc.df, aes(x = fire_sdc, y = observed_sdc)) + 
  geom_point(data = sdc.df, aes(x = fire_sdc, y = observed_sdc), color = "indianred3", size = 2) + 
  geom_line(data = sdc.df, aes(y = predicted_sdc), color = "indianred3", size = 1) + 
  geom_line(data = sdc.df, aes(x = fire_sdc, y = observed_sdc), color = "indianred3", size = 1) +
  geom_point(data = sdc.df1, aes(x = fire_sdc, y = observed_sdc), color = "turquoise3", size = 2) + 
  geom_line(data = sdc.df1, aes(y = predicted_sdc), color = "turquoise3", size = 1) +
  geom_line(data = sdc.df1, aes(x = fire_sdc, y = observed_sdc), color = "turquoise3", size = 1) +
  theme_bw() + 
  ggtitle("\npatch decay\n(SDC)") +
  #ggtitle("\nstand-replacing decay\ncoefficient (SDC)") + 
  labs(x="", y = "sdc") +
  #scale_y_continuous(breaks=seq(0,10,2)) +
  scale_x_continuous(breaks=seq(1985,2015,5)) +
  panel_border() 

sdc_combine <- sdc_combine + p
sdc_combine

###########################################################################################################
# Plot Patch Percentage of High-Severity Fire 
###########################################################################################################

# read in data and remove years with 1 fire 
mean_pland_sup <- read.csv("./mean_pland_sup.csv", header=TRUE, sep=",")
mean_pland_man <- read.csv("./mean_pland_man.csv", header=TRUE, sep=",")
mean_pland_sup <- mean_pland_sup[-c(1),]
mean_pland_man <- mean_pland_man[-c(1),]

# suppression fires
x_pland <- mean_pland_sup$YEAR
y_pland <- mean_pland_sup$PLAND
time_pland <- seq(from = 1, to = 33, by = 1)
time_pland <- time_pland - mean(time_pland)
time2_pland <- time_pland^2
myts_pland <- ts((mean_pland_sup$PLAND), start= 1985, end=2017)
#y_pland_per <- y_pland/100
pland.df <- data.frame(x_pland,y_pland)
#x2 <- x_pland^2
#lm.pland <- lm(y_pland ~ x_pland + x2)
fit_pland <- gls(myts_pland ~ time_pland + time2_pland, correlation = corARMA(q=1), method = "ML")  ##**
pland_fitted.df <- fitted(fit_pland)
pland.df <- data.frame(x = 1:33, fitted = pland_fitted.df)
pland.df <- cbind(pland.df, x_pland, y_pland)
names(pland.df) <- c("x_pland", 'predicted_pland', 'fire_pland', 'observed_pland')
names(pland.df)

# managed fires
x_pland1 <- mean_pland_man$YEAR
y_pland1 <- mean_pland_man$PLAND
yper1 <- y_pland1/100
#y_pland_per1 <- y_pland1/100
#pland.df1 <- data.frame(x_pland1,y_pland1)
#plot(x_pland1, y_pland1)
#lm.pland1 <- lm(y_pland1 ~ x_pland1)
pland_fit1 <- betareg(yper1 ~ x_pland1, data = mean_pland_man) #***************
pland_fitted.df1 <- fitted(pland_fit1)
pland.df1 <- data.frame(x = 1:21, fitted = pland_fitted.df1)
pland.df1 <- cbind(pland.df1, x_pland1, y_pland1)
names(pland.df1) <- c("x_pland", 'predicted_pland', 'fire_pland', 'observed_pland')
names(pland.df1)
pland.df1$predicted_pland <- pland.df1$predicted_pland * 100

# plot
pland_combine <- ggplot(data = pland.df, aes(x = fire_pland, y = observed_pland)) + 
  geom_point(data = pland.df, aes(x = fire_pland, y = observed_pland), color = "indianred3", size = 2) + 
  geom_line(data = pland.df, aes(y = predicted_pland), color = "indianred3", size = 1) + 
  geom_line(data = pland.df, aes(x = fire_pland, y = observed_pland), color = "indianred3", size = 1) +
  geom_point(data = pland.df1, aes(x = fire_pland, y = observed_pland), color = "turquoise3", size = 2) + 
  geom_line(data = pland.df1, aes(y = predicted_pland), color = "turquoise3", size = 1) +
  geom_line(data = pland.df1, aes(x = fire_pland, y = observed_pland), color = "turquoise3", size = 1) +
  theme_bw() + 
  ggtitle("\nPercent High-Severity Fire\n(PLAND)") +
  #ggtitle("\npercentage like\nadjacencies (pland)") + 
  labs(x="", y = "percent") +
  #scale_y_continuous(breaks=seq(0,10,2)) +
  scale_x_continuous(breaks=seq(1985,2015,5)) +
  panel_border()

pland_combine <- pland_combine + p
pland_combine

###########################################################################################################
# Plot Combined
###########################################################################################################

dev.off()
#cowplot
p.combine <- plot_grid(ca_combine, gyrate_combine, ndca_combine, 
               para_combine, pladj_combine, sdc_combine,
               labels = "AUTO", 
               label_size = 12,
               align = "h")
p.combine

OutPath <- "./ggplot_all_6_Time_Series.pdf"
ggsave(OutPath, plot=p.combine, device = cairo_pdf)



