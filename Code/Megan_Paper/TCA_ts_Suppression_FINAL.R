####################################################################################################
###    CODE FOR MEAN TCA - SUPPRESSION FIRES
####################################################################################################
library(forecast) #time series
library(tseries) #time series
library(dplyr) #formating
library(bbmle) #AICctab 
library(car) #durbinwatson test
library(ggplot2) #graphics
library(cowplot) #graphics
library(extrafont) #fonts
library(ggfortify) #autoplot
library(trend) #mktest
library(sjstats) #rmse
library(gvlma) #diagnostics
library(lmtest) #lrtest
library(nlme) #gls
font_import() #register fonts, must do once per session
loadfonts()
loadfonts(device = "win")


setwd("D:/Patch_Metrics/GEE_Severity_Dataset/R/Time_Series/")
ts_patch_metrics <- read.csv(file="GEE_1984_2017_final_patch_metrics_709.csv", header=TRUE, sep=",")
names(ts_patch_metrics)
attach(ts_patch_metrics)

####################################################################################################
###    MEAN TCA - Suppression Fires
####################################################################################################
filter.df <- filter(ts_patch_metrics, FIRE_TYPE2 == "SUPPRESSION")
names(filter.df)
attach(filter.df)

##aggregate data by TCA by year then calculate mean 
mean_tca.df <- aggregate(filter.df[, 30], list(filter.df[,2]), mean)
##rename columns
colnames(mean_tca.df)[colnames(mean_tca.df)=="Group.1"] <- "YEAR"
colnames(mean_tca.df)[colnames(mean_tca.df)=="x"] <- "TCA"
##convert character data to numeric 
mean_tca.df <- mutate_all(mean_tca.df, function(x) as.numeric(as.character(x)))
##format number of digits
mean_tca.df <- mean_tca.df %>% 
  mutate_if(is.numeric, round, digits = 5)
mean_tca.df
#remove 1984 - only 1 fire
mean_tca.df <- mean_tca.df[-c(1),]

#time series analysis 
myts <- ts((mean_tca.df$TCA), start=1985, end=2017)
myts.log <- ts(log(mean_tca.df$TCA + 1), start=1985, end=2017)
plot(myts)
plot(myts.log)

#test for autocorrelation at first lag using Durbin Watson test 
x <- mean_tca.df$YEAR
y <- mean_tca.df$TCA
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

########## Testing Linear Model Assumptions ####################################
x <- mean_tca.df$YEAR
y <- mean_tca.df$TCA
ylog <- log(y + 1)
fit1 <- lm(y ~ x)
fit2 <- lm(ylog ~ x)

#normality test
shapiro.test(resid(fit2))
#all diagnositics
autoplot(fit2)
gvlma(fit2)
summary(fit2)

##### GLS MODEL #################################################################
time.lm <- mean_tca.df$YEAR
time <- seq(from = 1, to = 33, by = 1)
time <- time - mean(time)
time2 <- time^2
gls.fit1 <- gls(myts ~ time.lm, correlation = corARMA(p=1), method = "ML")
gls.fit2 <- gls(myts.log ~ time.lm, correlation = corARMA(p=1), method = "ML")  ##****
gls.fit3 <- gls(myts.log ~ time.lm, correlation = corARMA(q=1), method = "ML")  

#Test competing models 
AICctab(gls.fit2, gls.fit3, weights = TRUE, delta = TRUE, base = TRUE)

#summary
summary(gls.fit1)
summary(gls.fit2) 

#liklihood ratio test 
lrtest(gls.fit1)
lrtest(gls.fit2)

#pseudo r squared
rcompanion::nagelkerke(gls.fit1)
rcompanion::nagelkerke(gls.fit2)

#plot residuals
plot(gls.fit1)
plot(gls.fit2) #**

#goodness of fit
rmse(gls.fit2)
AIC(gls.fit2, k=2)

########## PLotting Best Model ###############################################
x <- mean_tca.df$YEAR
y <- mean_tca.df$TCA
ylog <- log(y + 1)
yfit <- fitted(gls.fit2)
cols <- c("SUP" = "indianred3")

tca_ts <- ggplot(data = mean_tca.df, aes(x = YEAR, y = ylog)) + 
  geom_point(color= 'indianred3', size = 3) + 
  geom_line(color = "indianred3", size = 0.7) +
  geom_line(aes(y = yfit), color = "indianred3", linetype = 1, size = 1) +
  theme_bw() + 
  ggtitle("\nTotal Core Area in all Suppression Fires") + 
  labs(x="\nYEAR", y = "\nln(ha)") +
  #scale_y_continuous(breaks=seq(200, 800, 150)) +
  scale_x_continuous(breaks=seq(1985,2017,5)) +
  panel_border()

tca_ts <- tca_ts + theme(plot.title = element_text(family = "Calibri", color = "black", size =16, hjust = 0.5, lineheight = 0.8, face="bold"),
                           axis.title.x = element_text(family = "Calibri", color="black", size=12, face = "bold"),
                           axis.text.x = element_text(family = "Calibri", color="black", size=12),
                           axis.title.y = element_text(family = "Calibri", color="black", size=12, face = "bold"),
                           axis.text.y = element_text(family = "Calibri", color="black", size=12),
                           plot.margin = unit(c(1, 1.5, 1, 1.5), "cm"), #top, left, bottom, right
                           panel.border = element_rect(colour = "black", fill=NA, size=0.5))
tca_ts

OutPath <- "D:/Patch_Metrics/R/Time_Series/TCA_Time_Series/PDFs/TCA_ts_USFS_supp.pdf"
ggsave(OutPath, plot=tca_ts, device = cairo_pdf)

##################################################################################
##################################################################################