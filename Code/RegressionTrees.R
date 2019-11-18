library(tidyverse)
library(glmulti) #For glmulti(); version 1.0.7
library(rpart) #For rpart(); version 4.1-11
library(rpart.plot) #for prp(); version 2.1.2
library(RColorBrewer) #for brewer.pal(); version 1.1-2
library(grid) #for viewport(); version 3.4.3


####1. Read and process data####
#d <- read.csv("./Data/master_list_metrics.csv")
d <- read.csv("./Data/master_list_metrics_with_unknowns.csv") #Using unknowns file because same length as climate file so more straightforward merge. Can edit later.
d_climate <- read.csv("./Data/Processed Data/sdc_df_climate.csv")
names(d_climate)[1] <- "FIRE_NAME" #For consistent nomenclature

##Process master data file for merging
d_climate$FIRE_NAME <- as.character(d_climate$FIRE_NAME)
d$FIRE_NAME<-as.character(d$FIRE_NAME)
for(i in 1:nrow(d_climate)){
  d_climate$FIRE_NAME[i] <- gsub( "_1.*$", "", (d_climate$FIRE_NAME[i] )) 
  d_climate$FIRE_NAME[i] <- gsub( "_2.*$", "", (d_climate$FIRE_NAME[i] ))
    #This removes the extra characters after the fire name,
    #from fires in 1900's (first line) and fires in 2000's (second line)
}
for(i in 1:nrow(d)){
  d$FIRE_NAME[i] <- gsub( "_1.*$", "", (d$FIRE_NAME[i] )) 
  d$FIRE_NAME[i] <- gsub( "_2.*$", "", (d$FIRE_NAME[i] ))
  #And do the same for the master data frame, which seems to have some extraneous "_2"'s
  
}

d2 <-merge.data.frame(d,d_climate,by = "FIRE_NAME")
#Notes: 
#1)This is messy because I did it quickly. Adding in the climate data; the merge operation went from 575 lines to 512, need to figure out what happened to the missing files. *TODO
#2)The SDC values in the merged file don't match exactly, I wonder if there was a rounding error somewhere along the way. Would be worth a look. *TODO
#Problem fires:

#az3345810987920150819_CREEK 2015 (not in Jens' file with weather data)
#az3167810952819890704_SWISHELM (not in Jens' file with weather data)
#az3150310908920150617_HOG (mislabeled in Megan's data as "nm"?)

#Fires in Megan's dataset not in Jens':
d$FIRE_NAME[which(!d$FIRE_NAME%in%d_climate$FIRE_NAME)]
View(d[which(!d$FIRE_NAME%in%d_climate$FIRE_NAME),])
#Fires in Jens's dataset not in Megan's:
d_climate$FIRE_NAME[which(!d_climate$FIRE_NAME%in%d$FIRE_NAME)]
View(d_climate[which(!d_climate$FIRE_NAME%in%d$FIRE_NAME),])

####2. Test code for regression trees####
potential_parms <- 
  c("FIRE_TYPE1","FIRE_TYPE2","AGBUR1","max_tmmx","max_tmmn","min_rmax","max_bi")

####2a. Model selection algorithm.####

m_tph <- glmulti(y="SDC",  #m_tph can be changed, an artifact of imported code.
                 xr=potential_parms,
                 data=d,
                 level=1,method="h", plotty = FALSE) 
summary(m_tph@objects[[1]])
r.squaredGLMM(m_tph@objects[[1]])

#TODO: The above code not working on Jens' machine because of rJava issue, need to troubleshoot. Meantime, doing this:
f <- formula(SDC ~ FIRE_TYPE1 + AGBUR1 + max_tmmx + max_tmmn + min_rmax + max_bi)

#tree <- #Best model according to glmulti. Using this one to create regression tree 
#  rpart(m_tph@objects[[1]]$formula, 
#        method = "anova",
#        data = d,
#        control = rpart.control(cp = 0.02))
#levels(tree$frame$var) <- c("<leaf>", "deficit", "precip", "snowpack") #old code

tree <- #Best model according to glmulti. Using this one to create regression tree 
  rpart(f, 
        method = "anova",
        data = d,
        control = rpart.control(cp = 0.02))
#levels(tree$frame$var) <- c("<leaf>", "deficit", "precip", "snowpack")

h <- #Histogram
  ggplot(d)+
  geom_histogram(aes(TPH_Tot), binwidth = 30, center = 15,
                 fill=c(rep(brewer.pal(9,"RdYlBu")[1],2),
                        rep(brewer.pal(9,"RdYlBu")[2],1),
                        rep(brewer.pal(9,"RdYlBu")[3],1),
                        rep(brewer.pal(9,"RdYlBu")[5],1),
                        rep(brewer.pal(9,"RdYlBu")[7],1),
                        rep(brewer.pal(9,"RdYlBu")[8],4),
                        rep(brewer.pal(9,"RdYlBu")[9],6)
                 ),
                 col="black")+
  labs(x = "TPH", y = "Count")+
  coord_cartesian(xlim = c(0, 800)) +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5, size = 9),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 9))

#cairo_pdf("./Figures/El Dorado/Fig2a_v3.pdf", height = 3, width = 6.69)
prp(tree, varlen = 0, faclen = 0, type = 3, ge = " ≥ ", extra = 1, cex = 0.6, 
    #box.col = brewer.pal(9,"RdYlBu")[c(1,3,1,5,1,3,1,7,8)],
    #values match tree$frame[,"yval"]
    clip.right.labs = FALSE, #Only applies if using type = 3
    mar = c(2,2,2,2), 
    main = "SDC regression tree")
#print(h, vp=viewport(.18, .25, .38, .45))
#print(grid.text("a", x=unit(1, "npc"), y= unit (1, "npc"), 
#                vp=viewport(.01, .9, .1, .1) ) )
#dev.off()
  
  