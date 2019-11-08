## This code was developed by Sean Parks and has been modified by Jens Stevens to run on the SW Fires Database assembled by Megan Singleton.
## Date: 7/14/19

## Purpose: This is an attempt to recreate the methodology described in Stevens et al. (2017; Alternative characterization of forest fire regimes: incorperating spatial patterns; Landscape Ecology).
## Stevens et al. propose a fire-based severity metric, called the stand-replacing decay coefficient (SDC) that is concerned with the size and shape of high severity patches.
## However, I instead implement the entire process using the raster package in R. In contrast, the code distributed by Stevens et al. uses a polygon-based approach in R. 
## This code and approach is not necessarily better, but does not depend on converting rasters to polygons within or outside of R.
## I believe this code runs faster than the polygon-based code distributed by Stevens et al.


library(raster)
#library(FNN) #Deprecated?
#library(snow) #Deprecated?
library(tidyverse)
library(gridExtra)
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

## Set resolution of calculation
## 10 meters is slow; 30 meters is reasonably fast
the.res <- 30
options(scipen=999) ## turn off scientific number notation

####JTS specific code###
##Directory naming
wd <- c("../../../GIS/Fires/Southwest Burn Severity Dataset/Forested_Fires_Analysis_1984_2015")
yr <- 1984
#CRS is EPSG 102008 Albers Equal Area (but R doesn't recognize this EPSG code)

## Make a table for the resulting SDC for each fire
sdc.df <- as.data.frame(matrix(nrow=0, ncol=2)); colnames(sdc.df) <- c('fire.name','year')

for (yr in c(1984:2015)){
  wd_yr <- paste(wd,yr,"By_Fire","High_Severity","Forested_High_Severity",sep = "/")
  raster.list <- list.files(wd_yr)
  raster.list <- #subset files to only pull out "tif"s
    raster.list[grep("tif",substrRight(raster.list,3))]

for (j in 1:length(raster.list)) {#Calculate SDC for every fire in a given year
	fire.name <- raster.list[j]

	## Load the high severity raster image
	hs.raster.raw <- raster(paste(wd_yr,fire.name,sep = "/"))
	if(!any(!is.na(getValues(hs.raster.raw)))){#If no high severity (all NA)
	  print(paste(fire.name,"no high severity"))
	} else { #Everything below only happens if there is high severity.
	
	## Create a binary raster representing 'high' and 'other' severity
	## I used an RBR threshold of 298 to define high severity fire.
	## One could easily use dNBR or RdNBR and use a different threshold
	hs.raster <- hs.raster.raw
	hs.raster[is.na(hs.raster.raw)] <- 0  ## other severity
	hs.raster[hs.raster.raw > 0] <- 1  ## high severity 


	#### Start: Get rid of any blob <= 9 pixels  ####
####START HERE#### with e.g. the 2000	Cascabel fire (j=16)
	
	## In theory, this is optional
	## This can and should probably be converted to a function

	## This entire block of code is intended to mimic the 'nibble' function in ArcGIS.
	## I pretty sure someone (with better R skills than I) could clean this up
	## I apologize that I did not make any comments in this block of code

	u <- unique(hs.raster)
	for (i in 1:length(u)) {
		x <- hs.raster==u[i]
		clump.raster <- clump(x)
		rg.df <- na.omit(subset(as.data.frame(freq(clump.raster)), count > 0))
		rg.df$from <- rg.df$value - 0.25
		rg.df$to <- rg.df$value + 0.25
		rg.df$becomes = 1
		rg.df$becomes[which(rg.df$count <= 9)] <- 99   ## change this to value other than 9 if you want a different size threshold
		rg.df$becomes[which(rg.df$count > 9)] <- 1   ## change this to value other than 9 if you want a different size threshold
		rg.df <- as.matrix(rg.df[,c(3:5)]) 
		assign(paste('clump.raster.', u[i], sep=''), reclassify(clump.raster, rg.df))
	}

	the.stack <- stack()
	the.list <- list()
	for (n in 1:length(u)) {
		the.list <- append(the.list, get(noquote(paste0('clump.raster.', u[n]))))
		}
	the.stack <- stack(the.list)

	names(the.stack) <- u
	small.regions <- merge(the.stack)

	cell.num <- Which(small.regions >= 1, cells=T)
	xy <- as.data.frame(xyFromCell(small.regions, cell.num, spatial=F))
	nrow(xy)

	xy$extract <- extract(small.regions, xy)
	xy <- na.omit(xy)
	xy$sev <- extract(hs.raster, xy[,c(1,2)])
	xy$ID <- row.names(xy)

	nibble.df <- subset(xy, extract == 99)
	nibble.df <- nibble.df[,-c(3,4)]
	hs.df <- subset(xy, extract != 99)

	hs.df$ID <- seq(1, nrow(hs.df))
	hs.df.tmp <- hs.df[,-c(1,2,3)]
	hs.df.tmp$id.tmp <- hs.df.tmp$ID

	nibble.nn <- get.knnx(hs.df[,c(1,2)], cbind(nibble.df$x, nibble.df$y), 1, algorithm='brute')
	nibble.df$id.tmp <- nibble.nn$nn.index

	hs.df.tmp <- hs.df.tmp[,-c(2)]

	nibble.df <- merge(nibble.df, hs.df.tmp, by='id.tmp', all.x=T)
	nibble.df <- nibble.df[,-c(1,4)]

	hs.df <- rbind(hs.df[,c(1,2,4)], nibble.df)

	hs.raster <- rasterFromXYZ(hs.df, res=c(30, 30), digits=0, crs=projection(hs.raster))

	#### END: Get rid of any blob <= 9 pixels  ####
	
	if(max(unique(hs.raster))==0){#If no high severity after filter (all 0's)
	  print(paste(fire.name,"no high severity after 9 pixel patch filter"))
	} else { #Everything below only happens if there is high severity after filtering out small patches.
	
	## Makes a raster showing high severity fire only	
	hs.mask <- hs.raster
	hs.mask[hs.mask == 0] <- NA

	## Makes a raster showing 'other severity' fire only
	## This is necessary for buffering the interior of high severity patches
	hs.raster.reverse <- hs.raster
	hs.raster.reverse[hs.raster == 1] <- NA
	hs.raster.reverse[hs.raster == 0] <- 1
	
	## Create the distance raster using a single core
	dist.raster <- distance(hs.raster.reverse)
	
	## Create the distance raster using multiple cores (not working now)
	#beginCluster(4)
	#	dist.raster <- clusterR(hs.raster.reverse, distance)
	#endCluster()
	if(max(unique(dist.raster))<60){#If no high severity >60m from edge, can't fit curve
	  print(paste(fire.name,"not enough high severity to fit SDC curve after 9 pixel patch filter"))
	} else { #Everything below only happens if there is high severity.
	  
	## remove pixels that are not nigh severity	
	dist.raster <- mask(dist.raster, hs.mask)

	## This is the count of pixels and the total area defined as high severity
	hs.count <- freq(dist.raster > 0)[1,2]
	hs.area <- hs.count * 900 * 0.001 #in ha

	## Make a blank 'decay table'. This table, when populated, is used to calcuate SDC.
	decay.table.parks <- as.data.frame(matrix(nrow=0, ncol=2))
	colnames(decay.table.parks) <- c('width', 'prop.hs')

	## This figures out how many unique distance thresholds are necessary
	bins <- round(cellStats(dist.raster, max) / the.res)
	bin.list <- (0:(bins)) * the.res

	## The first entry is always these values
	decay.table.parks[1, 'width'] <- 0
	decay.table.parks[1, 'prop.hs'] <- 1

	## Populate the decay table with proportion of high severity for each distance increment.	
	for (k in 2:length(bin.list)) { 
		tmp.raster <- dist.raster
		tmp.raster[(dist.raster > bin.list[k])] <- 99
		tmp.raster[(dist.raster <= bin.list[k])] <- 0
			
		the.prop <- 1 - (freq(tmp.raster)[1,2] / hs.count)

		decay.table.parks[k, 'width'] <- bin.list[k]
		decay.table.parks[k, 'prop.hs'] <- the.prop
	}


	(SDC <- 
	   round(coef(nls(prop.hs~1/(10^(sdc*width)),data=decay.table.parks,start=list(sdc=0.01))), 4))

	sdc.df[nrow(sdc.df)+1, 'fire.name'] <- fire.name
	sdc.df[nrow(sdc.df), 'sdc'] <- SDC
	sdc.df[nrow(sdc.df), 'year'] <- yr

}#End "Else" portion on initial fire
}#End "Else" portion on post-9 pixel filtered fire
}#End "Else" portion if can't fit curve on post-9 pixel filtered fire
	print(paste(fire.name,"SDC calculated at", Sys.time() ))
}#End SDC "For" loop within year
}#End yearly "For" loop


####2. Post-processing####
for(i in 1:nrow(sdc.df)){
  sdc.df$name.simple[i] <- 
    paste(paste(sapply(strsplit(sdc.df$fire.name[i], "_"), function(x) x[c(2,3)])),collapse = "_")
}
  
write.csv(sdc.df,"./Data/Processed Data/sdc_df.csv")


####3. Plots and analysis####
p1 <-
  ggplot(sdc.df,aes(y = sdc, x = year)) + 
  geom_point() + 
  geom_smooth(method = "lm")+
  labs (title = "all fires")
p1b <-
  ggplot(sdc.df[sdc.df$year<=2013,],aes(y = sdc, x = year)) + 
  geom_point() + 
  geom_smooth(method = "lm")+
  labs (title = "all fires")

for(yr in unique(sdc.df$year)){
  df.tmp <- sdc.df[sdc.df$year==yr,]
  df.tmp25 <- df.tmp[which(df.tmp$sdc<quantile(df.tmp$sdc,0.25)),]
  df.tmp10 <- df.tmp[which(df.tmp$sdc<quantile(df.tmp$sdc,0.10)),]
  if(yr == 1984){
    sdc.yr.10pct <- df.tmp10; 
    sdc.yr.25pct <- df.tmp25; 
  } else{
    sdc.yr.10pct <- rbind(sdc.yr.10pct,df.tmp10)
    sdc.yr.25pct <- rbind(sdc.yr.25pct,df.tmp25)
  }
}  

p2<- 
  ggplot(sdc.yr.25pct,aes(y = sdc, x = year)) + 
  geom_point() + 
  geom_smooth(method = "lm")+
  labs(title = "bottom 25% sdc values per year")

p3<- 
  ggplot(sdc.yr.10pct,aes(y = sdc, x = year)) + 
  geom_point() + 
  geom_smooth(method = "lm")+
  labs(title = "bottom 10% sdc values per year")

p4<- 
  ggplot(sdc.yr.10pct[sdc.yr.10pct$year<=2013,],aes(y = sdc, x = year)) + 
  geom_point() + 
  geom_smooth(method = "lm")+
  labs(title = "bottom 10% sdc values per year\nthrough 2013")

pdf("./Figures/EDA/F1.pdf",height = 8, width = 4)
grid.arrange(p1,p2,p3,p4,ncol=1)
dev.off()
