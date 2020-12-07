###########################################################################################################
# 
# This script examines and plots fire size as a function of high-severity number of disjunct core areas for
# suppression and managed fires
# Final plotted models are zero-inflated negative binomial models
# Author: Megan Singleton 
# Updated: 12/2/2020
#
###########################################################################################################

library(ggplot2)
library(cowplot)
library(dplyr)
library(MASS)
library(bbmle)
library(EnvStats)
library(betareg)
library(pscl)
library(lmtest)
library(sjstats) #rmse
library(countreg)
library(ggeffects)
install.packages("countreg", repos="http://R-Forge.R-project.org")

# check to see if there are any unused packages 
funchir::stale_package_check('NDCA_Fire_Size_Regression.R')

# set working directory 
setwd("D:/Patch_Metrics/Patch_Metrics_Complete_Analysis/R/Fire_Size_Regression/")

# disable scientific notation
options(scipen=999)

# attach data
patch_metrics <- read.csv(file="GEE_1984_2017_final_patch_metrics_709.csv", header=TRUE, sep=",")
names(patch_metrics)

###########################################################################################################
# Analysis of NDCA by Fire Size - SUPPRESSION
###########################################################################################################

# filter dataset for all suppression fires
filter.df <- filter(patch_metrics, FIRE_TYPE2 == "SUPPRESSION")
names(filter.df)

# set x and y 
x_ndca <- filter.df$FORESTED_BURNED
y_ndca <- filter.df$NDCA
ndca.df <- data.frame(x_ndca,y_ndca)
plot(x_ndca, y_ndca)

# histogram of data
p <-ggplot(ndca.df, aes(x=x_ndca)) + 
  geom_histogram(binwidth = 1, color="black", fill="white")
p

# poisson model
glm.pois <- glm(y_ndca ~ x_ndca, data=ndca.df, family = poisson)
glm.pois.fit <- goodfit(y_ndca, type = "poisson")
summary(glm.pois)
pchisq(glm.pois$deviance, glm.pois$df.residual, lower.tail=F) # should get about 1 if data is not overdispersed
# data is overdispersed need better model

# negative binomial model 
glm.nb <- glm.nb(y_ndca ~ x_ndca, data=ndca.df, control=glm.control(maxit=100))
summary(glm.nb)

# zero-inflated poisson model 
glm.zinpois <- zeroinfl(y_ndca ~ x_ndca | x_ndca,
                        data = ndca.df,
                        dist = "poisson")
summary(glm.zinpois)

# zero-inflated negative binomial model 
glm.zinb  <- zeroinfl(y_ndca ~ x_ndca | x_ndca,
               data = ndca.df, dist = "negbin")
summary(glm.zinb)

# rootogram for choosing best model
root.pois <- countreg::rootogram(glm.pois, style = "hanging", plot = FALSE)
root.nb   <- countreg::rootogram(glm.nb, style = "hanging", plot = FALSE)
root.zinb <- countreg::rootogram(glm.zinb, style = "hanging", plot = FALSE)
root.zinpois <- countreg::rootogram(glm.zinpois, style = "hanging", plot = FALSE)
ylims <- ylim(-5, 25)  # common scale for comparison
plot_grid(autoplot(root.pois) + ylims, autoplot(root.nb) + ylims, 
          autoplot(root.zinpois) + ylims,
          autoplot(root.zinb) + ylims, ncol = 3, labels = "auto")

# test models - model selection
vuong(glm.nb,glm.zinb)

# likelihood ratio
lrtest(glm.zinb)
summary(glm.zinb)
rmse(glm.zinb)
AIC(glm.zinb)

# Pseudo R2 - Use McFadden
rcompanion::nagelkerke(glm.zinb)

###########################################################################################################
# Analysis of NDCA by Fire Size - SUPPRESSION
###########################################################################################################

# generate predicted means for suppression fires
plot(ggpredict(glm.zinb, terms="x_ndca"))

pr.zinb <- predict(glm.zinb)
pr.zinb <- round(pr.zinb, digits = 0)
pr.df <- data.frame(x = 1:501, Zinb = pr.zinb)
pr.df <- cbind(pr.df, x_ndca, y_ndca)
names(pr.df) <- c("x", 'predicted', 'fire_size', 'observed')

yfit <- fitted(glm.zinb)
cols <- c("SUP" = "indianred3")

ndca_predict <- ggplot(data = pr.df, aes(x = fire_size, y = observed)) +
  geom_point(data = pr.df, aes(x = fire_size, y = observed, colour = "SUP"), size = 2) +
  geom_line(data = pr.df, aes(y = yfit), color = "indianred", linetype = 1, size = 1) +
  theme_bw() +  
  ggtitle("\nNumber of Disjunct Core Areas > 200m from \nForest Edge in all Suppression Fires") + 
  labs(x="Foresed Burned (ha)", y = "\nCount") +
  scale_x_continuous(limits=c(0,50000)) +
  scale_y_continuous(limits=c(0,60)) +
  scale_colour_manual(name="Legend",values=cols, 
                      guide = guide_legend(override.aes=aes(fill=NA))) +
  panel_border()

p <- theme(plot.title = element_text(family = "Calibri", color = "black", size =12, hjust = 0.5, lineheight = 0.8, face="bold"),
      axis.title.x = element_text(family = "Calibri", color="black", size=12, face = "bold"),
      axis.text.x = element_text(family = "Calibri", color="black", size=12),
      axis.title.y = element_text(family = "Calibri", color="black", size=12, face = "bold"),
      axis.text.y = element_text(family = "Calibri", color="black", size=12),
      plot.margin = unit(c(1, 1.5, 1, 1.5), "cm"), #top, left, bottom, right
      panel.border = element_rect(colour = "black", fill=NA, size=0.5))

ndca_predict <- ndca_predict + p
ndca_predict

OutPath <- "./NDCA_supp_zinb.pdf"
ggsave(OutPath, plot=ndca_predict, device = cairo_pdf)


###########################################################################################################
# Analysis of NDCA by Fire Size - MANAGED
###########################################################################################################

# filter dataset for all managed fires
filter.df1 <- filter(patch_metrics, FIRE_TYPE2 == "OTHER")
names(filter.df1)

# set x and y
x_ndca1 <- filter.df1$FORESTED_BURNED
y_ndca1 <- filter.df1$NDCA
ndca.df1 <- data.frame(x_ndca1,y_ndca1)

# rescale FORESTED BURNED to resolve NaNs issues in zero inflated binomial regression model below
# rescale since Hessian matrix is not positive definite
x_ndca1 <- x_ndca1/1000
ndca.df1 <- data.frame(x_ndca1,y_ndca1)
plot(ndca.df1)

# poisson model
glm.pois1 <- glm(y_ndca1 ~ x_ndca1, data=ndca.df1, family = poisson) 
summary(glm.pois1)
pchisq(glm.pois1$deviance, glm.pois1$df.residual, lower.tail=F) # should get about 1 if data is not overdispersed

# negative binomial model 
glm.nb1 <- glm.nb(y_ndca1 ~ x_ndca1, data=ndca.df1, control=glm.control(maxit=100))
summary(glm.nb1)

# zero-inflated negative binomial model 
glm.zinb1 <- zeroinfl(y_ndca1  ~ x_ndca1 | x_ndca1,                     
               data = ndca.df1, dist = "negbin", EM = TRUE,   
               control = zeroinfl.control(maxit=1000))
summary(glm.zinb1)
# assess variance-covariance matrix
diag(glm.zinb$vcov)

# zero-inflated poisson model 
glm.zinpois1 <- zeroinfl(y_ndca1 ~ x_ndca1 | x_ndca1,
               data = ndca.df1,
               dist = "poisson", 
               EM = TRUE)
summary(glm.zinpois)

# rootogram analysis for model selection 
root.pois1 <- countreg::rootogram(glm.pois1, style = "hanging", plot = FALSE)
root.nb1   <- countreg::rootogram(glm.nb1, style = "hanging", plot = FALSE)
root.zinb1 <- countreg::rootogram(glm.zinb1, style = "hanging", plot = FALSE)
root.zinpois1 <- countreg::rootogram(glm.zinpois1, style = "hanging", plot = FALSE)
ylims <- ylim(-5, 25)  # common scale for comparison
plot_grid(autoplot(root.pois) + ylims, autoplot(root.nb) + ylims, 
          autoplot(root.zinb) + ylims, autoplot(root.zinpois) + ylims,
          ncol = 3, labels = "auto")

# test competing models
vuong(glm.zinb1, glm.zinpois)

# AIC/RMSE
summary(glm.zinb1)
rmse(glm.zinb1)
AIC(glm.zinb1, k=2)

# likelihood ratio against null model for p-value
lrtest(glm.zinb1)

# Pseudo R2 - use McFadden
rcompanion::nagelkerke(glm.zinb1)

###########################################################################################################
# Plotting NDCA by Fire Size - MANAGED
###########################################################################################################

# generate predicted means for managed fires
plot(ggpredict(glm.zinb1, terms="x_ndca1"))

pr.zinb1 <- predict(glm.zinb1)
pr.zinb1 <- round(pr.zinb1, digits = 0)
pr.df1 <- data.frame(x = 1:208, Zinb = pr.zinb1)
pr.df1 <- cbind(pr.df1, x_ndca1, y_ndca1)
names(pr.df1) <- c("x", 'predicted', 'fire_size', 'observed')

yfit1 <- fitted(glm.zinb1)
cols <- c("MAN" = "turquoise3")

ndca_predict2 <- ggplot(data = pr.df1, aes(x = fire_size, y = observed)) + 
  geom_point(data = pr.df1, aes(x = fire_size, y = observed, colour = "MAN"), size = 2) +
  theme_minimal() + 
  ggtitle("\nNumber of Disjunct Core Areas > 200m from \nForest Edge in all Managed Fires") + 
  geom_line(data = pr.df1, aes(y = yfit1), color = "turquoise3", linetype = 1, size = 1.3) +
  labs(x="Forested Burned", y = "Count") +
  #scale_y_continuous(breaks=seq(0,10,1)) +
  #scale_x_continuous(breaks=seq(0,25000,5000)) +
  panel_border() + 
  scale_color_manual(name = "class",
                     labels = c("MAN"),
                     values = cols)

ndca_predict2 <- ndca_predict2 + p

ndca_predict2

OutPath3 <- "D:/Patch_Metrics/GEE_Severity_Dataset/R/GEE_Regression_Analysis/NDCA/PDFs/NDCA_man_zinb.pdf"
ggsave(OutPath3, plot=ndca_predict2, device = cairo_pdf)


###########################################################################################################
# Plotting NDCA by Fire Size - COMBINED
###########################################################################################################

cols <- c("SUP" = "indianred3", "MAN" = "turquoise3")

ndca_combine <- ggplot(data = pr.df, aes(x = fire_size, y = observed)) + 
  geom_point(data = pr.df, aes(x = fire_size, y = observed, colour = "SUP"), size = 2, alpha = 0.5) + 
  geom_point(data = pr.df1, aes(x = fire_size, y = observed, colour = "MAN"), size = 2, alpha = 0.5) + 
  theme_bw() + 
  ggtitle("\nNumber of Disjunct Core Areas\nin all Suppression and Managed Fires") + 
  geom_line(data = pr.df, aes(y = yfit), color = "indianred3", linetype = 1, size = 1.3) +
  geom_line(data = pr.df1, aes(y = yfit1), color = "turquoise3", linetype = 1, size = 1.3) +
  labs(x="\nForested Burned Area", y = "\nNDCA") +
  scale_x_continuous(limits=c(0,50000)) +
  scale_y_continuous(limits=c(0,60)) +
  panel_border() + 
  scale_colour_manual(name="Legend",
                      values=cols, 
                      guide = guide_legend(override.aes=aes(fill=NA))) 

ndca_combine <- ndca_combine + p
ndca_combine

OutPath <- "./NDCA_sup_man_beta.pdf"
ggsave(OutPath, plot=ndca_combine, device = cairo_pdf)
