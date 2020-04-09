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

#setwd("D:/Patch_Metrics/GEE_Severity_Dataset/R/Time_Series/")
ts_patch_metrics <- 
  read.csv(file="./Data/GEE_1984_2017_final_patch_metrics_709.csv", header=TRUE, sep=",")
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
rmse(gls.fit3) #JTS note: requires modelr or another package not called above?
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


sdc_ts2 <- sdc_ts + 
  theme(plot.title = element_text(family = "Calibri", color = "black", size =16, hjust = 0.5, 
                                  lineheight = 0.8, face="bold"),
        axis.title.x = element_text(family = "Calibri", color="black", size=12, face = "bold"),
        axis.text.x = element_text(family = "Calibri", color="black", size=12),
        axis.title.y = element_text(family = "Calibri", color="black", size=12, face = "bold"),
        axis.text.y = element_text(family = "Calibri", color="black", size=12),
        plot.margin = unit(c(1, 1.5, 1, 1.5), "cm"), #top, left, bottom, right
        panel.border = element_rect(colour = "black", fill=NA, size=0.5))
sdc_ts


ggsave("./Figures/EDA/sdc_ts_baseline.pdf", plot=sdc_ts, device = cairo_pdf)

####JTS code starts here####
#Task 1: Figure out how many fires represented in each year for SDC, CA
for(r in 1:nrow(mean_sdc.df)){
  y <- mean_sdc.df$YEAR[r]
  mean_sdc.df$n[r] <- 
    length(which(filter.df$YEAR==y))
}
pdf("./Figures/JTS_troubleshoot/1.SDC~sample_size2.pdf")
ggplot(mean_sdc.df)+
  geom_point(aes(x = n, y = SDC))
dev.off()
pdf("./Figures/JTS_troubleshoot/2.sample_size~year2.pdf")
ggplot(mean_sdc.df)+
  geom_point(aes(x = YEAR, y = n))
dev.off()

#Observation: years with a lot (n>20 fires) tend to have middling SDC values and get pulled towards the central tendency. With more fires, the impact of the "bad ones" gets reduced when taking the mean, especially for a value like SDC which doesn't have a lot of range (compare to range of CA for example, where, the value for a "bad" fire is WAY different than for an average fire). In 2011, the famous Las Conchas had the second-lowest SDC (0.0039), but there were a TON of fires that year (44), and the median SDC is 0.015 and there are many values in the 0.02-0.03 range. And there is a trend towards more fires over time.

#Task 2: Look for anomalies in CA vs SDC (could replace CA with NDCA, and SDC with any of the other metrics, if interested, but need updated data frame)
filter.df$CA_pct <-percent_rank(filter.df$CA)
filter.df$SDC_pct <-percent_rank(filter.df$SDC)
pdf("./Figures/JTS_troubleshoot/3.anomalies_SDCvsCA2.pdf")
ggplot(filter.df)+
  geom_point(aes(x = CA_pct, y = SDC_pct,col = YEAR))
dev.off()
filter.df$SDC_CA_anomalies <-
  resid(lm(SDC_pct~CA_pct,data=filter.df))
ggplot(filter.df)+
  geom_point(aes(x = CA_pct, y = SDC_pct,col = SDC_CA_anomalies))
pdf("./Figures/JTS_troubleshoot/4.anomalies_SDC~year2.pdf")
ggplot(filter.df,aes(x = YEAR, y = SDC_CA_anomalies))+
  geom_point()+
  geom_smooth()
dev.off()
Low_SDC <- #These are the ten fires with suspiciously small (coarse-grained) SDC values for how much HS area they have. They are mostly smaller (bottom half of distribution) CA value fires. Perhaps they had very round HS patches, that can give smaller SDC values. But might be worth taking a look at them in a GIS.
  filter.df[order(filter.df$SDC_CA_anomalies),"FIRE_NAME"][c(1:10)]
High_SDC <- #These are the ten fires with suspiciously high (fine-grained) SDC values for how much HS area they have. They are mostly larger (top half of distribution) CA value fires. Might be worth taking a look at them in a GIS. Notable that *7 of these 10 fires* ocurred in 2011 and one ocurred in 2012.
  filter.df[order(filter.df$SDC_CA_anomalies, decreasing = T),"FIRE_NAME"][c(1:10)]

#Observation: It seems clear that the high number of fires in 2011 is pulling the SDC curve back "up" relative to the CA curve and giving it the quadratic shape. The 2011 and 2012 CA values in the manuscript are the reason that a linear model was the best fit and the slope is positive. Same for NDCA although I don't have those data in my older version of the database.

#Task 3: Look at temporal relationship in SDC but rather than annual means, use all data, median, and minimum values.
ggplot(mean_sdc.df,aes(x = YEAR, y = SDC))+ #Basically the original plot
  geom_point()+
  geom_smooth()+
  labs(y = "mean SDC")
ggplot(filter.df,aes(x = YEAR, y = SDC))+ #Plot with all data resembles the original plot
  geom_point()+
  geom_smooth()+
  labs(y = "all SDC")

median_sdc.df <- aggregate(filter.df[, 5], list(filter.df[,2]), median)
##rename columns
colnames(median_sdc.df)[colnames(median_sdc.df)=="Group.1"] <- "YEAR"
colnames(median_sdc.df)[colnames(median_sdc.df)=="x"] <- "SDC"
##convert character data to numeric 
median_sdc.df <- mutate_all(median_sdc.df, function(x) as.numeric(as.character(x)))
ggplot(median_sdc.df,aes(x = YEAR, y = SDC))+
  geom_point()+
  geom_smooth()+
  labs(y = "median SDC")

min_sdc.df <- aggregate(filter.df[, 5], list(filter.df[,2]), min)
##rename columns
colnames(min_sdc.df)[colnames(min_sdc.df)=="Group.1"] <- "YEAR"
colnames(min_sdc.df)[colnames(min_sdc.df)=="x"] <- "SDC"
##convert character data to numeric 
min_sdc.df <- mutate_all(min_sdc.df, function(x) as.numeric(as.character(x)))
pdf("./Figures/JTS_troubleshoot/5.min_SDC~year2.pdf")
ggplot(min_sdc.df,aes(x = YEAR, y = SDC))+
  geom_point()+
  geom_smooth()+
  labs(y = "minimum SDC")
dev.off()

#Observation: minimum SDC would give a negative linear trend because fires in 1980's and early 90's were mild. There does still seem to be an uptick in minimum SDC since 2013 as there haven't been any real megafires since then (at least in NM, the last really intense season was 2013 with Thompson Ridge, etc)