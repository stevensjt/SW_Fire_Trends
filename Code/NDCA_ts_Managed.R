###########################################################################################################
# 
# This script examines the number of disjunct core areas as a function of time for managed fires
# Author: Megan Singleton 
# Updated: 12/04/2020
#
###########################################################################################################

# load libraries
library(forecast)
library(tseries)
library(dplyr) #formating
library(bbmle) #AICctab 
library(car) #durbinwatson test
library(ggplot2) #graphics
library(cowplot) #graphics
library(extrafont) #fonts
library(pscl) #vuong
library(MASS) #neg binomial
library(trend) #mk test
library(lmtest) #lrtest
font_import() #register fonts, must do once per session
loadfonts()
loadfonts(device = "win")

# set working directory 
setwd("D:/Patch_Metrics/Patch_Metrics_Complete_Analysis/R/Time_Series/")

# check to see if there are any unused packages 
funchir::stale_package_check('NDCA_ts_Managed.R')

# disable scientific notation
options(scipen=999)
options(digits=6)

# attach data
ts_patch_metrics <- read.csv(file="GEE_1984_2017_final_patch_metrics_709.csv", header=TRUE, sep=",")
names(ts_patch_metrics)

###########################################################################################################
# SUM NDCA by YEAR - MANAGED
###########################################################################################################

# filter dataset by managed fires 
filter.df <- filter(ts_patch_metrics, FIRE_TYPE2 == "OTHER")

##aggregate data by NDCA by year then calculate mean 
# mean_ndca.df <- aggregate(filter.df[, 7], list(filter.df[,2]), mean)
sum_ndca.df <- aggregate(filter.df[, 32], list(filter.df[,2]), sum)
## rename columns
colnames(sum_ndca.df)[colnames(sum_ndca.df)=="Group.1"] <- "YEAR"
colnames(sum_ndca.df)[colnames(sum_ndca.df)=="x"] <- "NDCA"
## convert character data to numeric 
sum_ndca.df <- mutate_all(sum_ndca.df, function(x) as.numeric(as.character(x)))
## format to number of digits
sum_ndca.df <- sum_ndca.df %>% 
  mutate_if(is.numeric, round, digits = 0)
sum_ndca.df

# calculate number of fires per year
for(r in 1:nrow(sum_ndca.df)){
  y <- sum_ndca.df$YEAR[r]
  sum_ndca.df$n[r] <- 
    length(which(filter.df$YEAR==y))
}

# save
write.csv(sum_ndca.df, "./sum_ndca_man.csv")

# remove years with only one fire per year
sum_ndca.df <- sum_ndca.df[-c(1),]

# time series analysis 
myts <- ts((sum_ndca.df$NDCA), start=c(1995:2001, 2003:2014, 2016:2017))
myts.log <- ts(log(sum_ndca.df$NDCA + 1), start=c(1995:2001, 2003:2014, 2016:2017))
mk.test(myts, alternative = c("two.sided"),
        continuity = TRUE)

plot(myts)
list(myts)
plot(myts.log)
list(myts.log)

#test if autocorrelated using Durbin Watson test 
x <- sum_ndca.df$YEAR
y <- sum_ndca.df$NDCA
set.seed(123)
durbinWatsonTest(mylm)

#examine ACF and PACF plots for autocorrelation as well
acf(myts)
pacf(myts)
acf(myts.log)
pacf(myts.log)
acf(diff(myts.log))
pacf(diff(myts.log))
auto.arima(diff(myts.log))
auto.arima(myts)

###########################################################################################################
# Test OLS Assumptions
###########################################################################################################

# set x and y 
x <- sum_ndca.df$YEAR
y <- sum_ndca.df$NDCA
ylog <- log(y + 1)
mylm <- lm(y~x)
myloglm <- lm(ylog~x)

# fit candidate models
fit <- lm(y ~ x)
summary(fit)
plot(x, y)

# normality test
shapiro.test(resid(fit))
# variance test 
ncvTest(fit)
# all diagnositics
gvlma(fit)
autoplot(fit)

###########################################################################################################
# Binomial and Poisson Distribution and Zero Inflated Models
###########################################################################################################

#Poisson model 
summary(m1 <- glm(NDCA ~ YEAR, family = "poisson", data = sum_ndca.df))

#negative binomial regression 
summary(m2 <- glm.nb(NDCA ~ YEAR, data=sum_ndca.df, control=glm.control(maxit=100))) #*****
exp(0.112)

# quadratic term
x_ndca <- sum_ndca.df$YEAR
y_ndca <- sum_ndca.df$NDCA
x_ndca_c <- scale(x_ndca, center = TRUE, scale = FALSE)
x_ndca_2 <- x_ndca_c^2

m3 <- glm.nb(y_ndca ~ x_ndca + x_ndca_2, data=sum_ndca.df, control=glm.control(maxit=100))

#test models - model selection
vuong(m1,m2)
vuong(m2,m3)

#likelihood ratio
lrtest(m2)

#Pseudo R2
rcompanion::nagelkerke(m2)

#rmse
rmse(m2)

#AIC
AIC(m2, k=2)

###########################################################################################################
# Plot Best Model 
###########################################################################################################

m2 <- glm.nb(NDCA ~ YEAR, data=sum_ndca.df, control=glm.control(maxit=100))
yfit <- fitted(m2)
cols1 <- c("MAN" = "turquoise3")

ndca_ts <- ggplot(data = sum_ndca.df, aes(x = YEAR, y = NDCA)) + 
  geom_point(data = sum_ndca.df, aes(x = YEAR, y = NDCA, colour = "MAN"), size = 2) +
  #geom_point(color= 'turquoise3', size = 3) + 
  geom_line(color = "turquoise3", size = 0.7) +
  geom_line(aes(y = yfit), color = "turquoise3", linetype = 1, size = 1) +
  theme_bw() + 
  ggtitle("\nNumber of Disjunct Core Areas in All Managed Fires") + 
  labs(x="\nYEAR", y = "\nNDCA (count)") +
  scale_y_continuous(breaks=seq(0, 10, 2)) +
  scale_x_continuous(breaks=seq(1995,2017,5)) +
  #ggpubr::stat_cor(label.x = 1990, label.y = 0.20) + 
  panel_border() +
  scale_color_manual(name = "class",
                     labels = c("MAN"), 
                     values = cols1)

ndca_ts <- ndca_ts + theme(plot.title = element_text(family = "Calibri", color = "black", size =16, hjust = 0.5, lineheight = 0.8, face="bold"),
                           axis.title.x = element_text(family = "Calibri", color="black", size=12, face = "bold"),
                           axis.text.x = element_text(family = "Calibri", color="black", size=12),
                           axis.title.y = element_text(family = "Calibri", color="black", size=12, face = "bold"),
                           axis.text.y = element_text(family = "Calibri", color="black", size=12),
                           #plot.margin = unit(c(1, 1.5, 1, 1.5), "cm"), #top, left, bottom, right
                           panel.border = element_rect(colour = "black", fill=NA, size=0.5))
ndca_ts

# save
OutPath <- "D./NDCA_ts_man.pdf"
ggsave(OutPath, plot=ndca_ts, device = cairo_pdf)
