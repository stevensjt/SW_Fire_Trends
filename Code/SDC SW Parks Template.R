
## Author: Sean Parks; sean_parks@fs.fed.us
## Date: 3/1/2018

## Purpose: This is an attempt to recreate the methodology described in Stevens et al. (2017; Alternative characterization of forest fire regimes: incorperating spatial patterns; Landscape Ecology).
## Stevens et al. propose a fire-based severity metric, called the stand-replacing decay coefficient (SDC) that is concerned with the size and shape of high severity patches.
## However, I instead implement the entire process using the raster package in R. In contrast, the code distributed by Stevens et al. uses a polygon-based approach in R. 
## This code and approach is not necessarily better, but does not depend on converting rasters to polygons within or outside of R.
## I believe this code runs faster than the polygon-based code distributed by Stevens et al.

## One issue to note is that the SDC produced from this code does slightly differ from that of the polygon-based code of Stevens et al.
## I believe this to be because polygons tend to buffer differently (i.e. smoother) than rasters (because rasters are pixels).
## A table showing the SDC values for the polygon vs. raster based approach for four fires in the Northern Rocky mountains is shown below.


library(raster)
library(FNN)
library(snow)

## Set resolution of calculation
## 10 meters is slow; 30 meters is reasonably fast
the.res <- 10

setwd('D:/temp/stevens.et.al/')

options(scipen=999) ## turn off scientific number notation

raster.list <- list.files('sample.fires.sean')


## Make a table for the resulting SDC for each fire
sdc.df <- as.data.frame(matrix(nrow=0, ncol=1)); colnames(sdc.df) <- c('fire.name')

for (j in 1:length(raster.list)) {
	fire.name <- raster.list[j]

	## Load the severity image and mask out areas that are not within the fire perimeter
	mask <- raster(paste0('sample.fires.sean/', fire.name, '/Perimeter/fire_perimeter.tif'))
	rbr <- raster(paste0('sample.fires.sean/', fire.name, '/RBR/RBR.tif'))
	rbr <- crop(rbr, mask)
	mask <- crop(mask, rbr)
	rbr <- mask(rbr, mask)

	## Create a binary raster representing 'high' and 'other' severity
	## I used an RBR threshold of 298 to define high severity fire.
	## One could easily use dNBR or RdNBR and use a different threshold
	hs.raster <- rbr
	hs.raster[rbr < 298] <- 0  ## other severity
	hs.raster[rbr >= 298] <- 1  ## high severity 


	##############################################
	##############################################
	## Start: Get rid of any blob <= 9 pixels  ###
	##############################################
	##############################################

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

	hs.raster <- rasterFromXYZ(hs.df, res=c(30, 30), digits=0, crs=projection(rbr))

	##############################################
	##############################################
	## END: Get rid of any blob <= 9 pixels  ###
	##############################################
	##############################################

	if (the.res != 30) { hs.raster <- projectRaster(hs.raster, crs=projection(hs.raster), res=the.res, method='ngb') }

	## Makes a raster showing high severity fire only	
	hs.mask <- hs.raster
	hs.mask[hs.mask == 0] <- NA

	## Makes a raster showing 'other severity' fire only
	## This is necessary for buffering the interior of high severity patches
	hs.raster.reverse <- hs.raster
	hs.raster.reverse[hs.raster == 1] <- NA
	hs.raster.reverse[hs.raster == 0] <- 1

	## Create the distance raster using multiple cores
	beginCluster(4)
		dist.raster <- clusterR(hs.raster.reverse, distance)
	endCluster()

	## remove pixels that are not nigh severity	
	dist.raster <- mask(dist.raster, hs.mask)

	## This is the count of pixels and the total area defined as high severity
	hs.count <- freq(dist.raster > 0)[1,2]
	hs.area <- hs.count * 0.001

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

	###########################################################################################
	###########################################################################################
	## START of the Stevens et al calculation
	## If one was interested in comparing this code to that of stevens et al, run this block of code
	## Note: the 'decay' function at the bottom of this file needs to be run first
	## Another note: it is kinda slow, partially because converting the raster to poly is slow within R, and partially because the decay function takes a long time to run
	###########################################################################################
	###########################################################################################

	## Slow - I believe Stevens et al conducted this step in ArcMap or elsewhere where this conversion to polygon runs faster
	## hs.poly <- rasterToPolygons(hs.mask, dissolve=TRUE)

	## I used a 30 meter increment here for more of an apples-to-apples comparison
  	## decay.table.stevens <- decay(hs_fire=hs.poly,buf_max=1000, buf_inc=10, name=fire.name) # internal buffering at 10 m increments (from Functions.R) 

	## (SDC.stevens <- round(coef(nls(prop.hs~1/(10^(sdc*width)),data=decay.table.stevens,start=list(sdc=0.01))), 4))


	###########################################################################################
	###########################################################################################
	## END of the Stevens et al calculation
	###########################################################################################
	###########################################################################################

	(SDC.parks <- round(coef(nls(prop.hs~1/(10^(sdc*width)),data=decay.table.parks,start=list(sdc=0.01))), 4))

	sdc.df[j, 'fire.name'] <- fire.name
	sdc.df[j, 'sdc.parks'] <- SDC.parks
##	sdc.df[j, 'sdc.stevens'] <- SDC.stevens

}


## This shows a comparison of the SDC as calculated with polygons (Stevens et al) vs SDC calculated entirely with rasters (sdc.parks)
## using four fires in the northern Rocky Mountains.
## This is for 30-m pixels and 30-m buffer iterations

## sdf.df (30-m pixels and buffer iterations)
##          fire.name sdc.stevens sdc.parks
## 1        AHORN.2007      0.0027    0.0024
## 2       BRIDGE.2007      0.0058    0.0049
## 3 CANYON_CREEK.1988      0.0027    0.0024
## 4       HAMMER.2011      0.0057    0.0047


## Using 10-m pixels and a 10-m buffering increment
          fire.name sdc.parks
1        AHORN.2007    0.0026
2       BRIDGE.2007    0.0055
3 CANYON_CREEK.1988    0.0026
4       HAMMER.2011    0.0052








##  From the code distributed by Stevens et al.
##  Function:  decay
##  Decay profile: Implement internal buffering on high severity patches, create "Long" dataset
decay=function(hs_fire,buf_max,buf_inc,name,cancel=F){
  require(parallel) #for mclapply. version 3.3.2
  require(rgeos) #for gBuffer. version 0.3-21 
  require(raster) #for area. version 2.5-8
  
  #Set up long data frame:
  dist.table.sub <-
    data.frame(name = rep(name,(buf_max/buf_inc)+1),
               width=seq(0,buf_max,by=buf_inc), area_ha=NA 
               )
  #Core operation:
## I put mc.cores=1 because mclapply does not work on Windows (only Unix/Linus and Macs, I think)
## increase this value to speed up this function if you have the correct OS
  buf.list <- #Apply the buffer at every width in X; returns list of length X.
    mclapply(X=-dist.table.sub$width,
             FUN=gBuffer,
             spgeom=hs_fire, 
             byid = FALSE, id = NULL, quadsegs = 5, capStyle = "ROUND", 
             joinStyle = "ROUND", mitreLimit = 1, mc.cores = 1  
             ) 
  
  #Post-processing of long data frame:
  buf.list <- #Remove any NULL values (when buffer eliminated all HS areas)
    Filter(length,buf.list) 
  dist.table.sub$area_ha <- #Calculate area remaining at each buffer distance
    #0.0001 converts m2 to ha
    (as.vector(sapply(buf.list,area))*0.0001)[1:nrow(dist.table.sub)] 
  if(any(is.na(dist.table.sub$area_ha))){
    #If there are NULL values for the area because the buffer got too wide,
    #Set the first null value to 0:
    dist.table.sub$area_ha[which(is.na(dist.table.sub$area_ha))[1]] <- 0 
    if(any(is.na(dist.table.sub$area_ha))){
      #If there were more than one NULL value, remove the rest of them
      dist.table.sub <- #Delete the rest of the table rows with NA in area.
      dist.table.sub[-which(is.na(dist.table.sub$area_ha)),] 
    }
  }
  dist.table.sub$prop.hs <- 
    #Calculate the proportion of the original HS area remaining at each buffer distance
    dist.table.sub$area_ha/dist.table.sub$area_ha[1]
  #East Fire 23 sec w/ 10 m buf_inc; 12 sec w/ 20 m; 8 sec w/ 30 m; all accurate
  #Return the long data frame for the fire in question
  return(dist.table.sub)
}




