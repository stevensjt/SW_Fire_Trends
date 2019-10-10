library(stringi)
library(raster)
library(tidyverse)
library(rgdal) #readOGR, proj4string
library(rgeos)# For gCentroid
library(ncdf4) #for processing functions to handle .nc files

substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

df <- read_csv("./Data/Processed Data/sdc_df.csv")
df <- df[,-1] #remove first column (rownum)
df <- read_csv("./Data/Processed Data/sdc_df_climate.csv")

wd <- c("../../../../GIS/Fires/Southwest Burn Severity Dataset/Forested_Fires_Analysis_1984_2015")

for (r in c(525:575)){ #222 through 2001 works ok; issues in 2002 #nrow(df)
  yr <- as.numeric(gsub("^(?:[^_]+_){2}([^_]+).*", "\\1", df[r,1][[1]]))
  IGNITION_DATE <- substrRight(
    str_extract(df[r,1][[1]], "[^_]+"), 8
  )
  IGNITION_DATE_POSIX <- IGNITION_DATE  
  stri_sub(IGNITION_DATE_POSIX, 5, 4) <- "-"
  stri_sub(IGNITION_DATE_POSIX, 8, 7) <- "-"
  IGNITION_DATE_POSIX <- as.Date(IGNITION_DATE_POSIX)
  CONT_DATE_POSIX <- IGNITION_DATE_POSIX + 14
  CONT_DATE <- as.character(CONT_DATE_POSIX)
  CONT_DATE <- gsub("-", "", CONT_DATE)
  wd_yr <- paste(wd,yr,"By_Fire","High_Severity","Forested_High_Severity",sep = "/")
  fire.name <- df[r,"fire.name"]
  hs.raster.raw <- raster(paste(wd_yr,fire.name,sep = "/"))
  hs.poly.raw <- rasterToPolygons(hs.raster.raw, dissolve = TRUE)
  hs.poly.raw <- disaggregate(hs.poly.raw)
  hs_fire_ll=spTransform(hs.poly.raw,CRS("+proj=longlat +datum=WGS84")) #Reproject to LL to use as NOAA Arg
  
  #Download climate data
  sd=as.character(IGNITION_DATE_POSIX)
  ed=as.character(CONT_DATE_POSIX)
  
  vs=c("tmmx","tmmn","rmax","bi") #Download on "erc" not working properly, so skipping for now.
  vs_full=c("air_temperature","air_temperature","relative_humidity","burning_index_g") 
  #skipping "energy_release_component-g"
  for(index in 1:length(vs)){ #Variable for loop
    v=vs[index]
    v_full=vs_full[index]
    
    #Key function: Create a string with the url for data download
    #This string includes the specific area (lat/long) and time window of interest. 
    #The area is specified by a bounding box with a 0.2 degree buffer around the fire centroid.
    #The resolution of the data is 0.04 degrees per pixel-width.
    v_link=paste0("https://www.reacchpna.org/thredds/ncss/MET/",v,
                  "/",v,"_",as.character(yr),".nc?var=",v_full,
                  "&north=",gCentroid(hs_fire_ll)$y+0.2,
                  "&west=",gCentroid(hs_fire_ll)$x-0.2,
                  "&east=",gCentroid(hs_fire_ll)$x+0.2,
                  "&south=",gCentroid(hs_fire_ll)$y-0.2,"&disableProjSubset=on&horizStride=1",
                  "&time_start=",sd,"T00%3A00%3A00Z",
                  "&time_end=",ed,"T00%3A00%3A00Z&timeStride=1&accept=netcdf")
    
    dest <-  paste0("./Data/climate_nc/",v,".nc" )
    #fire.list[f,"lat"]=gCentroid(hs_fire_ll)$y
    #fire.list[f,"long"]=gCentroid(hs_fire_ll)$x
    tmp=try(download.file(url=v_link,destfile=dest), silent=T)
    if(class(tmp)!="try-error"){ #CHECK download error: 
      #if the download produced an error, will have to do it manually so skip the next bit.
      ncin <- nc_open(dest)
      lat <- ncvar_get(ncin,"lat",verbose=F)
      lat_target <- which.min(abs(lat - gCentroid(hs_fire_ll)$y)) #Find closest latitude to fire centroid
      lon <- ncvar_get(ncin,"lon",verbose=F)
      lon_target <- which.min(abs(lon - gCentroid(hs_fire_ll)$x)) #Find closest latitude to fire centroid
      v_array <- ncvar_get(ncin,v_full)
      if(class(v_array)=="matrix"){ #If there was only one burn day so there's a matrix instead of an array
        #Duplicate the arrayso there's no error produced in maximum calculation
        v_array=replicate(2,v_array,simplify="array")
      }
      if(v=="tmmx"){
        #Get maximum high temperature during the burn window
        df[r,paste0("max_",v)]=max(v_array[lat_target,lon_target,])-273.15 #Convert from K to C
      }
      if(v=="tmmn"){
        #Get maximum low temperature during the burn window
        df[r,paste0("max_",v)]=max(v_array[lat_target,lon_target,])-273.15 #Convert from K to C
      }
      #Get minimum high RH during the burn window
      if(v=="rmax"){
        df[r,paste0("min_",v)]=min(v_array[lat_target,lon_target,]) #Not sure of units
      }
      if(v=="bi"){
        #Get max burn index during burn window
        df[r,paste0("max_",v)]=max(v_array[lat_target,lon_target,]) #Not sure of units
      }
      if(v=="erc"){
        #Get max erc index during burn window
        df[r,paste0("max_",v)]=max(v_array[lat_target,lon_target,]) #Not sure of units
      }
    } #END CHECK download error
  } #END variable for loop
  gc()
  print(fire.name)
}
write_csv(df,"./Data/Processed Data/sdc_df_climate.csv")
