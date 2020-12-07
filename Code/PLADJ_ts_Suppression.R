###########################################################################################################
# 
# This script examines high-severity percentage of like adjacencies as a function of time for suppression fires
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
library(ggfortify) #autoplot
library(trend) #mktest
library(sjstats) #rmse
library(gvlma) #model diagnostics
library(nlme) #gls models
library(lmtest) #lr test
font_import() #register fonts, must do once per session
loadfonts()
loadfonts(device = "win")

# set working directory 
setwd("D:/Patch_Metrics/Patch_Metrics_Complete_Analysis/R/Time_Series/")

# check to see if there are any unused packages 
funchir::stale_package_check('PLADJ_ts_Suppression.R')

# disable scientific notation
options(scipen=999)
options(digits=6)

# attach data
ts_patch_metrics <- read.csv(file="GEE_1984_2017_final_patch_metrics_709.csv", header=TRUE, sep=",")
names(ts_patch_metrics)

###########################################################################################################
# MEAN PLADJ by YEAR - SUPPRESSION
###########################################################################################################

# filter dataset by all suppression fires
filter.df <- filter(ts_patch_metrics, FIRE_TYPE2 == "SUPPRESSION")

##aggregate data by PLADJ by year then calculate mean 
mean_pladj.df <- aggregate(filter.df[, 44], list(filter.df[,2]), mean)
##rename columns
colnames(mean_pladj.df)[colnames(mean_pladj.df)=="Group.1"] <- "YEAR"
colnames(mean_pladj.df)[colnames(mean_pladj.df)=="x"] <- "PLADJ"
##convert character data to numeric 
mean_pladj.df <- mutate_all(mean_pladj.df, function(x) as.numeric(as.character(x)))
##format to number of digits
mean_pladj.df <- mean_pladj.df %>% 
  mutate_if(is.numeric, round, digits = 5)
mean_pladj.df

# calculate the number of fires per year
for(r in 1:nrow(mean_pladj.df)){
  y <- mean_pladj.df$YEAR[r]
  mean_pladj.df$n[r] <- 
    length(which(filter.df$YEAR==y))
}

# save
write.csv(mean_pladj.df, "./mean_pladj_sup.csv")

# remove years with only 1 fire 
mean_pladj.df <- mean_pladj.df[-c(1),]

#time series analysis 
myts <- ts((mean_pladj.df$PLADJ), start= 1985, end=2017)
myts.assqt <- ts(asin(sqrt(mean_pladj.df$PLADJ/100)), start= 1985, end=2017)
myts.logit <- ts(logit(mean_pladj.df$PLADJ/100), start= 1985, end=2017)
plot(myts)
plot(myts.assqt)
plot(myts.logit)
mk.test(myts, alternative = c("two.sided"),
        continuity = TRUE)

# test for autocorrelation in first-lag using Durbin Watson test 
x <- mean_pladj.df$YEAR
y <- mean_pladj.df$PLADJ
yarc <- asin(sqrt(y/100))
mylm <- lm(y~x)
myarclm <- lm(yarc~x)
set.seed(123)
durbinWatsonTest(mylm)

# examine ACF and PACF plots for autocorrelation as well
acf(myts)
pacf(myts)
acf(myts.assqt)
pacf(myts.assqt)
auto.arima(myts)
auto.arima(myts.assqt)

###########################################################################################################
# If Autocorrelated - Use GLS Model
###########################################################################################################

# set x and y 
x <- mean_pladj.df$YEAR
y <- mean_pladj.df$PLADJ

# mean center to mitigate collinearity
time.lm <- mean_pladj.df$YEAR
time <- seq(from = 1, to = 33, by = 1)
time <- time - mean(time)
time2 <- time^2

# create candidate models
fit1 <- gls(myts ~ time.lm, correlation = corARMA(q=1), method = "ML")
fit2 <- gls(myts ~ time.lm, correlation = corARMA(p=1), method = "ML")
fit3 <- gls(myts ~ time + time2, correlation = corARMA(p=1), method = "ML")
fit4 <- gls(myts ~ time + time2, correlation = corARMA(q=1), method = "ML") 
fit5 <- gls(myts ~ time + time2, correlation = corARMA(q=1))
fit6 <- gls(myts.assqt ~ time + time2, correlation = corARMA(p=1), method = "ML")
fit7 <- gls(myts.logit ~ time + time2, correlation = corARMA(p=1), method = "ML") 
fit8 <- gls(myts.logit ~ time + time2, correlation = corARMA(q=1), method = "ML") ##****

# compare competing models 
bbmle::AICctab(fit7,fit8,base=TRUE,logLik=TRUE, weights=TRUE)

# summary statistics 
summary(fit8) #**

#lr test 
lrtest(fit8) #**

# pseudo R squared
rcompanion::nagelkerke(fit3)
rcompanion::nagelkerke(fit4) #**
rcompanion::nagelkerke(fit5)

# plot model residuals
plot(fit8)

# normality test 
shapiro.test(resid(fit6))
shapiro.test(resid(fit7))
shapiro.test(resid(fit8))

# goodness of fit
rmse(fit8)
AIC(fit8, k=2)

###########################################################################################################
# Plot Untransformed Model 
###########################################################################################################

y_logit <- car::logit(y)
fit4 <- gls(myts ~ time + time2, correlation = corARMA(q=1), method = "ML") 
yfit <- fitted(fit4)

cols <- c("SUP" = "indianred3")

pladj_ts <- ggplot(data = mean_pladj.df, aes(x = x, y = y)) + 
  geom_point(data = mean_pladj.df, aes(x = x, y = y, colour = "SUP"), size = 2) +
  geom_line(color = "indianred3", size = 0.7) +
  geom_line(aes(y = yfit), color = "indianred3", linetype = 1, size = 1) +
  theme_bw() + 
  ggtitle("\nPercent Like Adjacencies (PLADJ) \nin all Suppression Fires") + 
  labs(x="\nYEAR", y = "\nPercent") +
  #scale_y_continuous(breaks=seq(0, 90, 10)) +
  scale_x_continuous(breaks=seq(1984,2017,5)) +
  panel_border() + 
  scale_color_manual(name = "class",
                     labels = c("SUP"),
                     values = cols)

pladj_ts <- pladj_ts + theme(plot.title = element_text(family = "Calibri", color = "black", size =16, hjust = 0.5, lineheight = 0.8, face="bold"),
                             axis.title.x = element_text(family = "Calibri", color="black", size=12, face = "bold"),
                             axis.text.x = element_text(family = "Calibri", color="black", size=12),
                             axis.title.y = element_text(family = "Calibri", color="black", size=12, face = "bold"),
                             axis.text.y = element_text(family = "Calibri", color="black", size=12),
                             plot.margin = unit(c(1, 1.5, 1, 1.5), "cm"), #top, left, bottom, right
                             panel.border = element_rect(colour = "black", fill=NA, size=0.5))
pladj_ts

# save
OutPath <- "./PLADJ_ts_suppression.pdf"
ggsave(OutPath, plot=pladj_ts, device = cairo_pdf)


