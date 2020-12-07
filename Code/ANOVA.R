###########################################################################################################
# 
# This script performs Kruskal-Wallis one way ANOVA test on all 6 patch metrics
# It also creates boxplots for all 6 metrics for data visualization of differences between managed and 
# and suppression fires
# Author: Megan Singleton 
# Updated: 11/23/2020
#
###########################################################################################################

# load libraries
library(bbmle)
library(car)
library(dplyr)
library(cowplot)
library(ggplot2)
library(extrafont)
# import fonts for ggplot
font_import()
loadfonts()
loadfonts(device = "win")

# set working directory 
setwd("D:/Patch_Metrics/Patch_Metrics_Complete_Analysis/R/ANOVA/")

# check to see if there are any unused packages 
funchir::stale_package_check('ANOVA.R')

# disable scientific notation
options(scipen=999)

# load file
patch_metrics <- read.csv(file="GEE_1984_2017_final_patch_metrics_709.csv", header=TRUE, sep=",")
names(patch_metrics)

# kruskal wallis test for all patch metrics
kruskal.test(CA ~ FIRE_TYPE2, data = patch_metrics)
kruskal.test(GYRATE_AM ~ FIRE_TYPE2, data = patch_metrics)
kruskal.test(NDCA ~ FIRE_TYPE2, data = patch_metrics)
kruskal.test(PARA_AM ~ FIRE_TYPE2, data = patch_metrics)
kruskal.test(PLADJ ~ FIRE_TYPE2, data = patch_metrics)
filter.df <- filter(patch_metrics, SDC > 0)
kruskal.test(SDC ~ FIRE_TYPE2, data = filter.df)

# set theme
p <- theme(plot.title = element_text(family = "Calibri", color = "black", size =13, hjust = 0.5, lineheight = 0.8, face="bold"),
      axis.title.x = element_text(family = "Calibri", color="black", size=12),
      axis.text.x = element_text(family = "Calibri", color="black", size=12),
      axis.text.y = element_text(family = "Calibri", color="black", size=12),
      plot.margin = unit(c(0.1, 0.25, 0.1, 0.25), "cm"), #top #left #bottom #right
      panel.border = element_rect(colour = "black", fill=NA, size=0.5))


######################################################################################################
#
# Box Plots for Class Area
#
######################################################################################################

filter.df <- filter(patch_metrics, FIRE_TYPE3 == "SUPPRESSION")
filter.df1 <- filter(patch_metrics, FIRE_TYPE3 == "MANAGED")

log_ca <- log(filter.df$CA)
d <- cbind(filter.df, log_ca)
log_ca1 <- log(filter.df1$CA)
d1 <- cbind(filter.df1, log_ca1)

ca_combine <- ggplot(data = d, aes(x = FIRE_TYPE3, y = log_ca)) + 
  geom_boxplot(data = d, aes(x = FIRE_TYPE3, y = log_ca), fill = "indianred3", size = 1) + 
  geom_boxplot(data = d1, aes(x = FIRE_TYPE3, y = log_ca1), fill = "turquoise3", size = 1) +
  theme_minimal() + 
  theme(legend.position="none") + 
  ggtitle("\nhigh-severity area\n(CA)*") + 
  labs(x="", y = "ln(ha)") +
  panel_border() 

ca_combine <- ca_combine + p
ca_combine

######################################################################################################
#
# Box Plots for Patch Reach
#
######################################################################################################

log_gyrate <- log(filter.df$GYRATE_AM)
d <- cbind(filter.df, log_gyrate)
log_gyrate1 <- log(filter.df1$GYRATE_AM)
d1 <- cbind(filter.df1, log_gyrate1)

gyrate_combine <- ggplot(data = d, aes(x = FIRE_TYPE3, y = log_gyrate)) + 
  geom_boxplot(data = d, aes(x = FIRE_TYPE3, y = log_gyrate), fill = "indianred3", size = 1) + 
  geom_boxplot(data = d1, aes(x = FIRE_TYPE3, y = log_gyrate1), fill = "turquoise3", size = 1) +
  theme_minimal() + 
  theme(legend.position="none") + 
  ggtitle("\npatch reach\n(GYRATE_AM)*") +
  #ggtitle("\npatch radius of gyration\n(GYRATE_AM)*") + 
  labs(x="", y = "ln(m)") +
  panel_border() 

gyrate_combine <- gyrate_combine + p
gyrate_combine

######################################################################################################
#
# Box Plots for Number of Disjunct Core Areas 
#
######################################################################################################

log_ndca <- filter.df$NDCA + 1
log_ndca <- log(log_ndca)
d <- cbind(filter.df, log_ndca)
log_ndca1 <- filter.df1$NDCA + 1
log_ndca1 <- log(log_ndca1)
d1 <- cbind(filter.df1, log_ndca1)

filter.df4 <- filter(patch_metrics, FIRE_TYPE2 == "SUPPRESSION", NDCA > 0)
filter.df5 <- filter(patch_metrics, FIRE_TYPE2 == "OTHER", NDCA > 0)
log_ndca <- log(filter.df4$NDCA)
d4 <- cbind(filter.df4, log_ndca)
log_ndca1 <- log(filter.df5$NDCA)
d5 <- cbind(filter.df5, log_ndca1)

ndca_combine <- ggplot(data = d4, aes(x = FIRE_TYPE3, y = log_ndca)) + 
  geom_boxplot(data = d4, aes(x = FIRE_TYPE3, y = log_ndca), fill = "indianred3", size = 1) + 
  geom_boxplot(data = d5, aes(x = FIRE_TYPE3, y = log_ndca1), fill = "turquoise3", size = 1) +
  theme_minimal() + 
  theme(legend.position="none") +
  ggtitle("\nnumber of core areas\n(NDCA)*") + 
  #ggtitle("\nnumber of disjunct\ncore areas (NDCA)*") + 
  labs(x="", y = "ln(count)") +
  panel_border() 

ndca_combine <- ndca_combine + p
ndca_combine

######################################################################################################
#
# Box Plots for Patch Shape
#
######################################################################################################

log_para <- log(filter.df$PARA_AM)
d <- cbind(filter.df, log_para)
log_para1 <- log(filter.df1$PARA_AM)
d1 <- cbind(filter.df1, log_para1)

para_combine <- ggplot(data = d, aes(x = FIRE_TYPE3, y = log_para)) + 
  geom_boxplot(data = d, aes(x = FIRE_TYPE3, y = log_para), fill = "indianred3", size = 1) + 
  geom_boxplot(data = d1, aes(x = FIRE_TYPE3, y = log_para1), fill = "turquoise3", size = 1) +
  theme_minimal() + 
  theme(legend.position="none") + 
  #ggtitle("\nperimeter:area\n(PARA_AM)*") + 
  ggtitle("\nshape complexity\n(PARA_AM)*") +
  labs(x="", y = "ln(perimeter:area)") +
  panel_border() 

para_combine <- para_combine + p
para_combine

######################################################################################################
#
# Box Plots for Patch Aggregation
#
######################################################################################################

logit_pladj <- logit(filter.df$PLADJ)
d <- cbind(filter.df, logit_pladj)
logit_pladj1 <- logit(filter.df1$PLADJ)
d1 <- cbind(filter.df1, logit_pladj1)

pladj_combine <- ggplot(data = d, aes(x = FIRE_TYPE3, y = logit_pladj)) + 
  geom_boxplot(data = d, aes(x = FIRE_TYPE3, y = logit_pladj), fill = "indianred3", size = 1) + 
  geom_boxplot(data = d1, aes(x = FIRE_TYPE3, y = logit_pladj1), fill = "turquoise3", size = 1) +
  theme_minimal() + 
  theme(legend.position="none") + 
  #ggtitle("\npercentage like\nadjacencies (PLADJ)*") + 
  ggtitle("\naggregation\n(PLADJ)*") +
  labs(x="", y = "logit(percent)") +
  panel_border() 

pladj_combine <- pladj_combine + p
pladj_combine

######################################################################################################
#
# Box Plots for Patch Decay
#
######################################################################################################

filter.df2 <- filter(patch_metrics, FIRE_TYPE2 == "SUPPRESSION", SDC > 0)
filter.df3 <- filter(patch_metrics, FIRE_TYPE2 == "OTHER", SDC > 0)

log_sdc <- log(filter.df2$SDC)
d2 <- cbind(filter.df2, log_sdc)
log_sdc1 <- log(filter.df3$SDC)
d3 <- cbind(filter.df3, log_sdc1)

sdc_combine <- ggplot(data = d2, aes(x = FIRE_TYPE3, y = log_sdc)) + 
  geom_boxplot(data = d2, aes(x = FIRE_TYPE3, y = log_sdc), fill = "indianred3", size = 1) + 
  geom_boxplot(data = d3, aes(x = FIRE_TYPE3, y = log_sdc1), fill = "turquoise3", size = 1) +
  theme_minimal() + 
  theme(legend.position="none") + 
  #ggtitle("\nstand-replacing decay\ncoefficient(SDC)*") + 
  ggtitle("\npatch decay\n(SDC)*") +
  labs(x="", y = "ln(sdc)") +
  panel_border() 

sdc_combine <- sdc_combine + p
sdc_combine

######################################################################################################
#
# Combined Box Plots for all 6 Metrics
#
######################################################################################################

dev.off()

# combine all plots with cowplot
p.combine <- plot_grid(ca_combine, gyrate_combine, ndca_combine, 
               para_combine, pladj_combine, sdc_combine,
               labels = "AUTO", label_size = 12,
               align = "h")
p.combine

OutPath <- "./anova_ggplot.pdf"
ggsave(OutPath, plot=p.combine, device = cairo_pdf)




