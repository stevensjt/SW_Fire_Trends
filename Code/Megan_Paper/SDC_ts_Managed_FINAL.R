####################################################################################################
###    CODE FOR MEAN SDC - MANAGED FIRES
####################################################################################################
library(forecast) #time series
library(tseries) #time series
library(dplyr) #formating
library(bbmle) #AICctab 
library(car) #durbinwatson test, non-constant variance
library(ggplot2) #graphics
library(cowplot) #graphics
library(extrafont) #fonts
library(gvlma) #test linear model assumptions
library(nlme) #generalized least squares
library(lmtest) #lr test
library(ggfortify) #autoplot
library(rcompanion) #rsquared
library(trend) #mk test 
font_import() #register fonts, must do once per session
loadfonts()
loadfonts(device = "win")

setwd("D:/Patch_Metrics/GEE_Severity_Dataset/R/Time_Series/")
ts_patch_metrics <- read.csv(file="GEE_1984_2017_final_patch_metrics_709.csv", header=TRUE, sep=",")
names(ts_patch_metrics)
attach(ts_patch_metrics)

####################################################################################################
###    MEAN SDC by MANAGED FIRES
####################################################################################################
filter.df <- filter(ts_patch_metrics, FIRE_TYPE2 == "OTHER", SDC > 0)
names(filter.df)
attach(filter.df)

##aggregate data by SDC by year then calculate mean 
mean_sdc.df <- aggregate(filter.df[, 5], list(filter.df[,2]), mean)
##rename columns
colnames(mean_sdc.df)[colnames(mean_sdc.df)=="Group.1"] <- "YEAR"
colnames(mean_sdc.df)[colnames(mean_sdc.df)=="x"] <- "SDC"
##convert character data to numeric 
mean_sdc.df <- mutate_all(mean_sdc.df, function(x) as.numeric(as.character(x)))
##format number of digits
mean_sdc.df <- mean_sdc.df %>% 
  mutate_if(is.numeric, round, digits = 5)
mean_sdc.df
mean_sdc.df <- mean_sdc.df[-c(1),]

#time series analysis 
myts <- ts((mean_sdc.df$SDC), start=c(1994:2001, 2003:2014, 2016:2017))
myts.log <- ts(log(mean_sdc.df$SDC), start=c(1994:2001, 2003:2014, 2016:2017))
plot(myts)
plot(myts.log)

#test for autocorrelation at first lag using Durbin Watson test 
x <- mean_sdc.df$YEAR
y <- mean_sdc.df$SDC
mylm <- lm(y~x)
set.seed(123)
durbinWatsonTest(mylm)

# log transform to stabilize variance 
ylog <- log(y + 1)
myloglm <- lm(ylog~x)
set.seed(123)
durbinWatsonTest(myloglm)

#man-kendall test for trends - if not autocorrelated
mk.test(myts, alternative = c("two.sided"),
        continuity = TRUE)

#examine ACF and PACF plots for autocorrelation beyond first lag
acf(myts)
pacf(myts)
acf(myts.log)
pacf(myts.log)
auto.arima(myts)
auto.arima(myts.log)

##### Testing Linear Model Assumptions #######################################
x <- mean_sdc.df$YEAR
y <- mean_sdc.df$SDC
#mean-center for collinearity of time
time <- seq(from = 1, to = 20, by = 1)
time <- time - mean(time)
time2 <- time^2
ylog <- log(y)
fit1 <- lm(y ~ x)
fit2 <- lm(ylog ~ x)
fit3 <- lm(y ~ time + time2)
fit4 <- lm(ylog ~ time + time2) #****

#test competing models
bbmle::AICctab(fit2,fit4,base=TRUE,logLik=TRUE, weights=TRUE)

#normality test
shapiro.test(resid(fit4))
#all diagnositics
autoplot(fit4)
gvlma(fit4)
#summary
summary(fit4)

########## PLotting Best Model ###############################################
cols <- c("MAN" = "turquoise3")
yfit <- fitted(fit4)

sdc_ts <- ggplot(data = mean_sdc.df, aes(x = x, y = ylog)) +
  geom_point(data = mean_sdc.df, aes(x = x, y = ylog, colour = "MAN"), size = 2) +
  geom_line(data = mean_sdc.df, aes(x = x, y = ylog), color = "turquoise3", size = 1) + 
  geom_line(data = mean_sdc.df, aes(y = yfit), color = "turquoise3", size = 1) +
  theme_bw() + 
  ggtitle("\nSDC in all Managed Fires") + 
  labs(x="\nYEAR", y = "\nln(SDC)") +
  scale_x_continuous(breaks=seq(1985,2017,5)) +
  panel_border() + 
  scale_color_manual(name = "class",
                     labels = c("MAN"),
                     values = cols)


sdc_ts <- sdc_ts + theme(plot.title = element_text(family = "Calibri", color = "black", size =16, hjust = 0.5, lineheight = 0.8, face="bold"),
                         axis.title.x = element_text(family = "Calibri", color="black", size=12, face = "bold"),
                         axis.text.x = element_text(family = "Calibri", color="black", size=12),
                         axis.title.y = element_text(family = "Calibri", color="black", size=12, face = "bold"),
                         axis.text.y = element_text(family = "Calibri", color="black", size=12),
                         plot.margin = unit(c(1, 1.5, 1, 1.5), "cm"), #top, left, bottom, right
                         panel.border = element_rect(colour = "black", fill=NA, size=0.5))
sdc_ts

OutPath <- "D:/Patch_Metrics/GEE_Severity_Dataset/R/Time_Series/SDC_Time_Series/PDFs/SDC_ts_MAN.pdf"
ggsave(OutPath, plot=sdc_ts, device = cairo_pdf)
