library(tidyverse)
library(glmulti) #For glmulti(); version 1.0.7
library(rpart) #For rpart(); version 4.1-11
library(rpart.plot) #for prp(); version 2.1.2
library(RColorBrewer) #for brewer.pal(); version 1.1-2
library(grid) #for viewport(); version 3.4.3
library(dplyr)
library(MuMIn)
library(ggplot2)
library(MASS)
Sys.setenv(JAVA_HOME='C:\\Program Files\\Java\\jdk1.8.0_231')
install.packages("rJava","http://rforge.net")
library(rJava)
library(extrafont)
font_import()
loadfonts()
loadfonts(device = "win")

options(scipen=999)

#### 1a. Read and process data ####
setwd("D:/Patch_Metrics/GEE_Severity_Dataset/")

d <- read.csv("./R/Climate_Analysis/GEE_1984_2017_final_patch_metrics.csv") 
d_climate <- read.csv("./R/Climate_Analysis/sdc_gee_1984_2017_climate_1_593.csv", header=TRUE, sep=",")
d_climate <- d_climate[,-(1)] #remove first column 

#Process master data file for merging
d_climate$FIRE_NAME <- as.character(d_climate$FIRE_NAME)
d$FIRE_NAME<-as.character(d$FIRE_NAME)

d2 <-merge.data.frame(d,d_climate,by = "FIRE_NAME")
names(d2)
#write.csv(d2, "./R/Climate_Analysis/gee_master_list_version1.csv")


#### 1b. clean up data frame by removing unneccessary columns for analysis ####
d3 <- d2 %>%
  dplyr::select(-(13:15), -(17:20))
d4 <- d3[, -(11)]

#### 1c. Transform response variable (SDC), add to data frame, and 
####    visualize distribution of the patch metric with histogram   
####     NOTE: the transformed/untransformed SDC returns the same best model in part 2a. 
####     But the transformations do not make SDC distribution normal  ####
sdc_log <- log(d4[,5])
sdc_sq <- sqrt(d4[,5])
d5 <- cbind(d4, sdc_log)
d5 <- cbind(d4, sdc_sq)
names(d5)

ggplot(d5, aes(x=SDC.x)) +
  geom_histogram(color = "black", fill = "white", binwidth = 0.0035)

ggplot(d5, aes(x=sdc_log)) + 
  geom_histogram(color = "black", fill = "white", binwidth = 0.25)

ggplot(d5, aes(x=sdc_sq)) +
  geom_histogram(color = "black", fill = "white", binwidth = 0.015)

#### 1d. Rename variables for prettier regression tree ####
colnames(d5)[colnames(d5)=="YEAR.x"] <- "year"
colnames(d5)[colnames(d5 )=="FIRE_TYPE2"] <- "class"
colnames(d5)[colnames(d5 )=="AGBUR2"] <- "agency"
colnames(d5)[colnames(d5 )=="MAX_TMMX"] <- "max_high_temp"
colnames(d5)[colnames(d5 )=="MAX_TMMN"] <- "max_low_temp"
colnames(d5)[colnames(d5 )=="MIN_RMAX"] <- "min_high_rh"
colnames(d5)[colnames(d5 )=="MAX_BI"] <- "max_burning_index"
colnames(d5)[colnames(d5 )=="MAX_ERC"] <- "max_erc"
names(d5)

#### 2a. Code for regression trees ####
potential_parms <- 
  c("year", "class","agency","max_high_temp","max_low_temp","min_high_rh","max_burning_index", "max_erc")

####  Model selection algorithm. ####
m_sdc <- glmulti(y="sdc_log",  
                 xr=potential_parms,
                 data=d5,
                 level=1,method="h", plotty = FALSE) 
summary(m_sdc@objects[[1]])
r.squaredGLMM(m_sdc@objects[[1]])

tree <- #Best model according to glmulti. Using this one to create regression tree 
  rpart(m_sdc@objects[[1]]$formula, 
        method = "anova",
        data = d5,
        control = rpart.control(cp = 0.02))


##### 3a. Visualization of analyses - Histogram  ####

h <- 
  ggplot(d5) +
  geom_histogram(aes(sdc_log),binwidth = 0.25, center = -4,
                 fill=c(rep(brewer.pal(9,"RdYlBu")[1],3),  
                        rep(brewer.pal(9,"RdYlBu")[2],2),
                        rep(brewer.pal(9,"RdYlBu")[3],2),
                        rep(brewer.pal(9,"RdYlBu")[5],3),
                        rep(brewer.pal(9,"RdYlBu")[7],1),
                        rep(brewer.pal(9,"RdYlBu")[8],1),
                        rep(brewer.pal(9,"RdYlBu")[9],1)  
                 ),
                 col="black")+
  labs(x = "SDC", y = "Count")+
  coord_cartesian(xlim = c(-7, -3)) +
  theme_bw()+
  theme(plot.title = element_text(family = "Calibri", color = "black", size =12, hjust = 0.5, lineheight = 0.8, face="bold"),
        axis.title.x = element_text(family = "Calibri", color="black", size=12),
        axis.title.y = element_text(family = "Calibri", color="black", size=12),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        plot.margin = unit(c(1, 0.5, 1, 0.5), "cm"))
h

OutPath <- "./R/Climate_Analysis/PDFs/SDC_histogram_test.pdf"
ggsave(OutPath, plot=h, device = cairo_pdf)

#### 3b. Visualization of analyses - Regression Tree ####

tree_viz <- prp(tree, varlen = 0, faclen = 0, type = 3, extra = 101, cex = 0.9, 
    box.col = brewer.pal(9,"RdYlBu")[c(5,5,5,5,2)], #change colors based on histogram
    clip.right.labs = FALSE, #Only applies if using type = 3
    mar = c(2,2,2,2), 
    main = "SDC regression tree")


#### 3c. Merge histogram and regression tree on same plot ####

tree_hist <- print(h, vp=viewport(x = 0.8, y = .20, width = .35, height = .40)) # location and size of histogram figure 
#print(grid.text("a", x=unit(1, "npc"), y= unit (1, "npc"), 
#                vp=viewport(.01, .9, .1, .1) ) )
#dev.off()

OutPath2 <- "./R/Climate_Analysis/PDFs/SDC_histogram_regression_tree.pdf"
ggsave(OutPath2, plot=tree_hist, device = cairo_pdf)
