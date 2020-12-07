###########################################################################################################
# 
# This script plots 6 different patch metrics as a function of proportion high-severity fire.
# Including: class area, patch radius of gyration, perimeter area ratio, number of core areas, percentage of like
# adjacencies and the stand decay coefficient
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
setwd("D:/Patch_Metrics/Patch_Metrics_Complete_Analysis/R/Proportion_Regression/")

# check to see if there are any unused packages 
funchir::stale_package_check('ggplot_all_6_Proportion_HS_Regression.R')

# disable scientific notation
options(scipen=999)

# attach data
patch_metrics <- read.csv(file="GEE_1984_2017_final_patch_metrics_709.csv", header=TRUE, sep=",")
names(patch_metrics)

# filter by suppression fires only
filter.df <- filter(patch_metrics, FIRE_TYPE2 == "SUPPRESSION")

# filter by managed fires only
filter.df1 <- filter(patch_metrics, FIRE_TYPE2 == "OTHER")

# theme for all plots
p <- theme(plot.title = element_text(family = "Calibri", color = "black", size =13, hjust = 0.5, lineheight = 0.8, face="bold"),
           axis.title.x = element_text(family = "Calibri", color="black", size=12),
           axis.title.y = element_text(family = "Calibri", color="black", size=12),
           #axis.title.y = element_blank(),
           plot.margin = unit(c(0.1, 0.25, 0.1, 0.25), "cm"), #top #left #bottom #right
           panel.border = element_rect(colour = "black", fill=NA, size=0.5))

###########################################################################################################
# Plot Class Area (CA)
###########################################################################################################

options(scipen = 0)
options(digits = 6)

# set x and y for suppression fires
x_ca <- filter.df$PLAND
y_ca <- filter.df$CA

# transform data and put in dataframe
x_logit_ca <- car::logit(x_ca)
y_log_ca <- log(y_ca)
ca.df <- data.frame(x_ca,y_log_ca)
plot(ca.df)

# best model
ca_logit.model <- lm(y_log_ca ~ x_logit_ca, data = ca.df)             

# create new dataframe of predicted and observed values for suppression fires
ca_fitted.df <- fitted(ca_logit.model)
ca.df <- data.frame(x = 1:501, fitted = ca_fitted.df)
ca.df <- cbind(ca.df, x_logit_ca, y_log_ca)
names(ca.df) <- c("x_ca", 'predicted_ca', 'fire_ca', 'observed_ca')

# set x and y for managed fires
x_ca1 <- filter.df1$PLAND
y_ca1 <- filter.df1$CA

# transform data and put in dataframe
x_logit_ca1 <- car::logit(x_ca1)
y_log_ca1 <- log(y_ca1)
ca.df1 <- data.frame(x_ca1,y_ca1)
plot(ca.df1)

# best model
ca_logit.model1 <- lm(y_log_ca1 ~ x_logit_ca1, data = ca.df1) 

# create new dataframe of predicted and observed values for managed fires
ca_fitted.df1 <- fitted(ca_logit.model1)
ca.df1 <- data.frame(x = 1:208, fitted = ca_fitted.df1)
ca.df1 <- cbind(ca.df1, x_logit_ca1, y_log_ca1)
names(ca.df1) <- c("x_ca", 'predicted_ca', 'fire_ca', 'observed_ca')

# plot
ca_combine <- ggplot(data = ca.df, aes(x = fire_ca, y = observed_ca)) + 
  geom_point(data = ca.df, aes(x = fire_ca, y = observed_ca), color = "indianred3", size = 2, alpha = 0.5) + 
  geom_line(data = ca.df, aes(y = predicted_ca), color = "indianred3", size = 1) + 
  geom_point(data = ca.df1, aes(x = fire_ca, y = observed_ca), color = "turquoise3", size = 2, alpha = 0.5) + 
  geom_line(data = ca.df1, aes(y = predicted_ca), color = "turquoise3", size = 1) +
  theme_bw() + 
  ggtitle("\nhigh-severity area\n(CA)*") + 
  labs(x="logit(percent high-severity)", y = "ln(ha)") +
  #scale_y_continuous(breaks=seq(0,10,2)) +
  #scale_x_continuous(breaks=seq(2,12,2)) +
  panel_border() 

ca_combine <- ca_combine + p
ca_combine

###########################################################################################################
# Plot Patch Radius of Gyration (GYRATE_AM)
###########################################################################################################

# set x and y for suppression fires
x_gyrate <- filter.df$PLAND
y_gyrate <- filter.df$GYRATE_AM

# transform data and put in dataframe
x_logit_gyrate <- car::logit(x_gyrate)
y_log_gyrate <- log(y_gyrate)
gyrate.df <- data.frame(x_logit_gyrate,y_log_gyrate)
plot(gyrate.df)

# best model 
gyrate_logit.model <- lm(y_log_gyrate ~ x_logit_gyrate + I(x_logit_gyrate^2), data = gyrate.df) #####******

# create new dataframe of predicted and observed values for suppression fires
gyrate_fitted.df <- fitted(gyrate_logit.model)
gyrate.df <- data.frame(x = 1:501, fitted = gyrate_fitted.df)
gyrate.df <- cbind(gyrate.df, x_logit_gyrate, y_log_gyrate)
names(gyrate.df) <- c("x_gyrate", 'predicted_gyrate', 'fire_gyrate', 'observed_gyrate')

# set x and y for managed fires
x_gyrate1 <- filter.df1$PLAND
y_gyrate1 <- filter.df1$GYRATE_AM

# transform data and put in dataframe
x_logit_gyrate1 <- car::logit(x_gyrate1)
y_log_gyrate1 <- log(y_gyrate1)
gyrate.df1 <- data.frame(x_logit_gyrate1,y_log_gyrate1)
plot(gyrate.df1)

# best model
gyrate_logit.model1 <- lm(y_log_gyrate1 ~ x_logit_gyrate1 + I(x_logit_gyrate1^2), data = gyrate.df1) #####******

# create new dataframe of predicted and observed values for managed fires
gyrate_fitted.df1 <- fitted(gyrate_logit.model1)
gyrate.df1 <- data.frame(x = 1:208, fitted = gyrate_fitted.df1)
gyrate.df1 <- cbind(gyrate.df1, x_logit_gyrate1, y_log_gyrate1)
names(gyrate.df1) <- c("x_gyrate", 'predicted_gyrate', 'fire_gyrate', 'observed_gyrate')

# plot
gyrate_combine <- ggplot(data = gyrate.df, aes(x = fire_gyrate, y = observed_gyrate)) + 
  geom_point(data = gyrate.df, aes(x = fire_gyrate, y = observed_gyrate), color = "indianred3", size = 2, alpha = 0.5) + 
  geom_line(data = gyrate.df, aes(y = predicted_gyrate), color = "indianred3", size = 1) + 
  geom_point(data = gyrate.df1, aes(x = fire_gyrate, y = observed_gyrate), color = "turquoise3", size = 2, alpha = 0.5) + 
  geom_line(data = gyrate.df1, aes(y = predicted_gyrate), color = "turquoise3", size = 1) +
  theme_bw() + 
  #ggtitle("\npatch radius of gyration\n(GYRATE_AM)*") + 
  ggtitle("\npatch reach\n(GYRATE_AM)*") +
  labs(x="logit(percent high-severity)", y = "ln(m)") +
  #scale_y_continuous(breaks=seq(0,10,2)) +
  #scale_x_continuous(breaks=seq(2,12,2)) +
  panel_border()

gyrate_combine <- gyrate_combine + p
gyrate_combine

###########################################################################################################
# Plot Number of Disjunct Core Areas (NDCA)
###########################################################################################################

# set x and y for suppression fires and put into dataframe
x_ndca <- filter.df$PLAND
y_ndca <- filter.df$NDCA
ndca.df <- data.frame(x_ndca,y_ndca)
plot(x_ndca, y_ndca)

# best model 
x_ndca_c <- scale(x_ndca, center = TRUE, scale = FALSE)
x_ndca_2 <- x_ndca_c^2
m_ndca  <- zeroinfl(y_ndca ~ x_ndca + x_ndca_2 | x_ndca,
               data = ndca.df, dist = "negbin") #####******

# create new dataframe of predicted and observed values for suppression fires
ndca_fitted.df <- fitted(m_ndca)
ndca.df <- data.frame(x = 1:501, fitted = ndca_fitted.df)
ndca.df <- cbind(ndca.df, x_ndca, y_ndca)
names(ndca.df) <- c("x_ndca", 'predicted_ndca', 'fire_ndca', 'observed_ndca')

# set x and y and put into dataframe
x_ndca1 <- filter.df1$PLAND
y_ndca1 <- filter.df1$NDCA
ndca.df1 <- data.frame(x_ndca1,y_ndca1)
plot(x_ndca1, y_ndca1)

# best model
m1_ndca <- zeroinfl(y_ndca1 ~ x_ndca1 | x_ndca1,
                    data = ndca.df1, dist = "negbin") #####******

# create new dataframe of predicted and observed values for managed fires
ndca_fitted.df1 <- fitted(m1_ndca)
ndca.df1 <- data.frame(x = 1:208, fitted = ndca_fitted.df1)
ndca.df1 <- cbind(ndca.df1, x_ndca1, y_ndca1)
names(ndca.df1) <- c("x_ndca", 'predicted_ndca', 'fire_ndca', 'observed_ndca')

# plot
ndca_combine <- ggplot(data = ndca.df, aes(x = fire_ndca, y = observed_ndca)) + 
  geom_point(data = ndca.df, aes(x = fire_ndca, y = observed_ndca), color = "indianred3", size = 2, alpha = 0.5) + 
  geom_line(data = ndca.df, aes(y = predicted_ndca), color = "indianred3", size = 1) + 
  geom_point(data = ndca.df1, aes(x = fire_ndca, y = observed_ndca), color = "turquoise3", size = 2, alpha = 0.5) + 
  geom_line(data = ndca.df1, aes(y = predicted_ndca), color = "turquoise3", size = 1) +
  theme_bw() + 
  #ggtitle("\nnumber of disjunct\ncore areas (NDCA)*") + 
  ggtitle("\nnumber of core areas\n(NDCA)*") +
  labs(x="percent high-severity", y = "count") +
  #scale_y_continuous(breaks=seq(0,10,2)) +
  scale_y_continuous(limits =c(0,60)) +
  #scale_x_continuous(breaks=seq(2,12,2)) +
  panel_border()

ndca_combine <- ndca_combine + p
ndca_combine
# 6 missing values

###########################################################################################################
# Plot Patch Perimeter:Ratio Ratio (PARA_AM)
###########################################################################################################

# set x and y for suppression fires
x_para <- filter.df$PLAND
y_para <- filter.df$PARA_AM

# transform data and put in dataframe
x_logit_para <- car::logit(x_para)
y_log_para <- log(y_para)
para.df <- data.frame(x_logit_para,y_log_para)
plot(para.df)

# best model
para_logit.model <- lm(y_log_para ~ x_logit_para + I(x_logit_para^2), data = para.df) #####****** 

# create new dataframe of predicted and observed values for suppression fires
para_fitted.df <- fitted(para_logit.model)
fitted.df <- fitted(para_logit.model)
para.df <- data.frame(x = 1:501, fitted = para_fitted.df)
para.df <- cbind(para.df, x_logit_para, y_log_para)
names(para.df) <- c("x_para", 'predicted_para', 'fire_para', 'observed_para')

# set x and y for managed fires
x_para1 <- filter.df1$PLAND
y_para1 <- filter.df1$PARA_AM

# transform data and put in dataframe
x_logit_para1 <- car::logit(x_para1)
y_log_para1 <- log(y_para1)
para.df1 <- data.frame(x_logit_para1,y_log_para1)
plot(para.df1)

# best model
para_logit.model1 <- lm(y_log_para1 ~ x_logit_para1 + I(x_logit_para1^2), data = para.df1) #####****** 

# create new dataframe of predicted and observed values for managed fires
para_fitted.df1 <- fitted(para_logit.model1)
para.df1 <- data.frame(x = 1:208, fitted = para_fitted.df1)
para.df1 <- cbind(para.df1, x_logit_para1, y_log_para1)
names(para.df1) <- c("x_para", 'predicted_para', 'fire_para', 'observed_para')

# plot
para_combine <- ggplot(data = para.df, aes(x = fire_para, y = observed_para)) + 
  geom_point(data = para.df, aes(x = fire_para, y = observed_para), color = "indianred3", size = 2, alpha = 0.5) + 
  geom_line(data = para.df, aes(y = predicted_para), color = "indianred3", size = 1) + 
  geom_point(data = para.df1, aes(x = fire_para, y = observed_para), color = "turquoise3", size = 2, alpha = 0.5) + 
  geom_line(data = para.df1, aes(y = predicted_para), color = "turquoise3", size = 1) +
  theme_bw() + 
  #ggtitle("\nperimeter:area\n(PARA_AM)*") +
  ggtitle("\nshape complexity\n(PARA_AM)*") +
  labs(x="logit(percent high-severity)", y = "ln(perimeter:area)") +
  #scale_y_continuous(breaks=seq(0,10,2)) +
  #scale_x_continuous(breaks=seq(2,12,2)) +
  panel_border() 


para_combine <- para_combine + p
para_combine

###########################################################################################################
# Plot Patch Percentage of Like Adjacencies 
###########################################################################################################

# plot x and y 
x_pladj <- filter.df$PLAND
y_pladj <- filter.df$PLADJ
y_pladj_per <- y_pladj/100
pladj.df <- data.frame(x_pladj,y_pladj_per)

# best model 
beta.pladj <- betareg(y_pladj_per ~ x_pladj, data = pladj.df, link = "loglog") #####******

# create new dataframe of predicted and observed values for managed fires
pladj_fitted.df <- fitted(beta.pladj)
pladj.df <- data.frame(x = 1:501, fitted = pladj_fitted.df)
pladj.df <- cbind(pladj.df, x_pladj, y_pladj_per)
names(pladj.df) <- c("x_pladj", 'predicted_pladj', 'fire_pladj', 'observed_pladj')

# plot x and y 
x_pladj1 <- filter.df1$PLAND
y_pladj1 <- filter.df1$PLADJ
y_pladj_per1 <- y_pladj1/100
pladj.df1 <- data.frame(x_pladj1,y_pladj1)
plot(x_pladj1, y_pladj_per1)

# best model
beta.pladj1 <- betareg(y_pladj_per1 ~ x_pladj1, data = pladj.df1, link = "loglog") #####******

# create new dataframe of predicted and observed values for managed fires
pladj_fitted.df1 <- fitted(beta.pladj1)
pladj.df1 <- data.frame(x = 1:208, fitted = pladj_fitted.df1)
pladj.df1 <- cbind(pladj.df1, x_pladj1, y_pladj_per1)
names(pladj.df1) <- c("x_pladj", 'predicted_pladj', 'fire_pladj', 'observed_pladj')

# plot
pladj_combine <- ggplot(data = pladj.df, aes(x = fire_pladj, y = observed_pladj)) + 
  geom_point(data = pladj.df, aes(x = fire_pladj, y = observed_pladj), color = "indianred3", size = 2, alpha = 0.5) + 
  geom_line(data = pladj.df, aes(y = predicted_pladj), color = "indianred3", size = 1) + 
  geom_point(data = pladj.df1, aes(x = fire_pladj, y = observed_pladj), color = "turquoise3", size = 2, alpha = 0.5) + 
  geom_line(data = pladj.df1, aes(y = predicted_pladj), color = "turquoise3", size = 1) +
  theme_bw() + 
  #ggtitle("\npercentage like\nadjacencies (PLADJ)*") + 
  ggtitle("\naggregation\n(PLADJ)*") + 
  labs(x="percent high-severity", y = "proportion") +
  #scale_y_continuous(breaks=seq(0,10,2)) +
  #scale_x_continuous(breaks=seq(2,12,2)) +
  panel_border()

pladj_combine <- pladj_combine + p


pladj_combine

###########################################################################################################
# Plot High-Severity Stand-Decay Coefficient 
###########################################################################################################

# filter by suppression fires only with SDC value
filter.df2 <- filter(patch_metrics, FIRE_TYPE2 == "SUPPRESSION", SDC > 0)

# filter by managed fires only with SDC value
filter.df3 <- filter(patch_metrics, FIRE_TYPE2 == "OTHER", SDC > 0)

# set x and y for suppression fires
x_sdc <- filter.df2$PLAND
y_sdc <- filter.df2$SDC

# transform data and put into dataframe
x_logit_sdc <- car::logit(x_sdc)
y_log_sdc <- log(y_sdc)
sdc.df <- data.frame(x_logit_sdc,y_log_sdc)
plot(x_logit_sdc, y_log_sdc)

# best model
sdc_logit.model <- lm(y_log_sdc ~ x_logit_sdc + I(x_logit_sdc^2), data = sdc.df) #####******

# create new dataframe of predicted and observed values for suppression fires
sdc_fitted.df <- fitted(sdc_logit.model)
sdc.df <- data.frame(x = 1:416, fitted = sdc_fitted.df)
sdc.df <- cbind(sdc.df, x_logit_sdc, y_log_sdc)
names(sdc.df) <- c("x_sdc", 'predicted_sdc', 'fire_sdc', 'observed_sdc')

# set x and y for managed fires
x_sdc1 <- filter.df3$PLAND
y_sdc1 <- filter.df3$SDC

# transform data and put into dataframe
x_logit_sdc1 <- car::logit(x_sdc1)
y_log_sdc1 <- log(y_sdc1)
sdc.df1 <- data.frame(x_sdc1,y_sdc1)
plot(x_logit_sdc1, y_log_sdc1)

# best model
sdc_logit.model1 <- lm(y_log_sdc1 ~ x_logit_sdc1 + I(x_logit_sdc1^2), data = sdc.df1) #####******

# create new dataframe of predicted and observed values for managed fires
sdc_fitted.df1 <- fitted(sdc_logit.model1)
sdc.df1 <- data.frame(x = 1:139, fitted = sdc_fitted.df1)
sdc.df1 <- cbind(sdc.df1, x_logit_sdc1, y_log_sdc1)
names(sdc.df1) <- c("x_sdc", 'predicted_sdc', 'fire_sdc', 'observed_sdc')

# plot
sdc_combine <- ggplot(data = sdc.df, aes(x = fire_sdc, y = observed_sdc)) + 
  geom_point(data = sdc.df, aes(x = fire_sdc, y = observed_sdc), color = "indianred3", size = 2, alpha = 0.5) + 
  geom_line(data = sdc.df, aes(y = predicted_sdc), color = "indianred3", size = 1) + 
  geom_point(data = sdc.df1, aes(x = fire_sdc, y = observed_sdc), color = "turquoise3", size = 2, alpha = 0.5) + 
  geom_line(data = sdc.df1, aes(y = predicted_sdc), color = "turquoise3", size = 1) +
  theme_bw() + 
  #ggtitle("\nstand-replacing decay\ncoefficient (SDC)*") +
  ggtitle("\npatch decay\n(SDC)*") +
  labs(x="logit(percent high-severity)", y = "ln(sdc)") +
  #scale_y_continuous(breaks=seq(0,10,2)) +
  #scale_x_continuous(breaks=seq(2,12,2)) +
  panel_border() 

sdc_combine <- sdc_combine + p
sdc_combine

###########################################################################################################
# Plot Combined
###########################################################################################################

dev.off()
# cowplot
p.combine <- plot_grid(ca_combine, gyrate_combine, ndca_combine, 
               para_combine, pladj_combine, sdc_combine,
               labels = "AUTO", label_size = 12,
               align = "h")
p.combine

OutPath <- "./ggplot_all_6_Proportion_HS_Regression.pdf"
ggsave(OutPath, plot=p, device = cairo_pdf)
# need to export as 6x8

