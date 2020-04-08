####################################################################################################
###    CODE FOR MEAN SDC - SUPPRESSION FIRES
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
###    MEAN SDC by Suppression Fires
####################################################################################################
filter.df <- filter(ts_patch_metrics, FIRE_TYPE2 == "SUPPRESSION", SDC > 0)
names(filter.df)
attach(filter.df)

##aggregate data by SDC by year then calculate mean 
mean_sdc.df <- aggregate(filter.df[, 5], list(filter.df[,2]), mean)
##rename columns
colnames(mean_sdc.df)[colnames(mean_sdc.df)=="Group.1"] <- "YEAR"
colnames(mean_sdc.df)[colnames(mean_sdc.df)=="x"] <- "SDC"
##convert character data to numeric 
mean_sdc.df <- mutate_all(mean_sdc.df, function(x) as.numeric(as.character(x)))
##format to number of digits
mean_sdc.df <- mean_sdc.df %>% 
  mutate_if(is.numeric, round, digits = 5)
mean_sdc.df
mean_sdc.df <- mean_sdc.df[-c(1),]

#time series analysis 
myts <- ts((mean_sdc.df$SDC), start=1985, end= 2017)
myts.log <- ts(log(mean_sdc.df$SDC), start=1985, end= 2017)
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
ylog <- log(y)
fit1 <- lm(y ~ x)
fit2 <- lm(ylog ~ x)
x2 <- x^2
fit3 <- lm(y ~ x + x2)
fit4 <- lm(ylog ~ x + x2)

#normality test
shapiro.test(resid(fit4))
#all diagnositics
autoplot(fit4)
gvlma(fit4)

bbmle::AICctab(fit2,fit5,base=TRUE,logLik=TRUE, weights=TRUE)

####################### GLS MODEL #############################################################
x <- mean_sdc.df$YEAR
y <- mean_sdc.df$SDC
#mean-center for collinearity of time
time <- seq(from = 1, to = 33, by = 1)
time <- time - mean(time)
time2 <- time^2

#gls models
gls.fit1 <- gls(myts.log ~ x, correlation = corARMA(p=1), method = "ML")
gls.fit2 <- gls(myts.log ~ x, correlation = corARMA(q=1), method = "ML")
gls.fit3 <- gls(myts.log ~ time + time2, correlation = corARMA(p=1), method = "ML") #*****
gls.fit4 <- gls(myts.log ~ time + time2, correlation = corARMA(q=1), method = "ML")

bbmle::AICctab(gls.fit3,gls.fit4,base=TRUE,logLik=TRUE, weights=TRUE)

#summary
summary(gls.fit3)
summary(gls.fit4) 

#liklihood ratio test 
lrtest(gls.fit3)
lrtest(gls.fit4)

#pseudo r squared
rcompanion::nagelkerke(gls.fit3)
rcompanion::nagelkerke(gls.fit4)

#plot residuals
plot(gls.fit3) #**
plot(gls.fit4) 

#goodness of fit
rmse(gls.fit3)
AIC(gls.fit3, k=2)


########## PLotting Best Model ###############################################
y <- mean_sdc.df$SDC
ylog <- log(y)
gls.fit3 <- gls(myts.log ~ time + time2, correlation = corARMA(p=1), method = "ML") #*****
yfit <- fitted(fit3)

cols <- c("SUP" = "indianred3")

sdc_ts <- ggplot(data = mean_sdc.df, aes(x = x, y = y)) +
  geom_point(data = mean_sdc.df, aes(x = x, y = y, colour = "SUP"), size = 2) +
  geom_line(data = mean_sdc.df, aes(x = x, y = y), color = "indianred3", size = 1) + 
  geom_line(aes(y = yfit), color = "indianred3", linetype = 1, size = 1) +
  theme_bw() + 
  ggtitle("\nSDC in all Suppression Fires") + 
  labs(x="\nYEAR", y = "\nln(SDC)") +
  #scale_y_continuous(breaks=seq(-7, -3, 0.25)) +
  scale_x_continuous(breaks=seq(1985,2017,5)) +
  scale_colour_manual(name="Legend",values=cols, 
                      guide = guide_legend(override.aes=aes(fill=NA))) +
  panel_border()


sdc_ts <- sdc_ts + theme(plot.title = element_text(family = "Calibri", color = "black", size =16, hjust = 0.5, lineheight = 0.8, face="bold"),
                         axis.title.x = element_text(family = "Calibri", color="black", size=12, face = "bold"),
                         axis.text.x = element_text(family = "Calibri", color="black", size=12),
                         axis.title.y = element_text(family = "Calibri", color="black", size=12, face = "bold"),
                         axis.text.y = element_text(family = "Calibri", color="black", size=12),
                         plot.margin = unit(c(1, 1.5, 1, 1.5), "cm"), #top, left, bottom, right
                         panel.border = element_rect(colour = "black", fill=NA, size=0.5))
sdc_ts

OutPath <- "D:/Patch_Metrics/GEE_Severity_Dataset/R/Time_Series/SDC_Time_Series/PDFs/SDC_ts_supp_709.pdf"
ggsave(OutPath, plot=sdc_ts, device = cairo_pdf)
