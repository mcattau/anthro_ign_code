
#######################################################################################################################################
### Author: Megan Cattau ##### 
## Earth Lab, Project Forest
## Contact info: megan.cattau@gmail.com or megan.cattau@colorado.edu, 706.338.9436
# Manuscript authors: Megan Cattau, Carol Wessman, Adam Mahood, Jennifer Balch
## Project: Fire trends and anthro ignitions
# Project description: This project aims to elucidate how fire physical characteristics are changing over time in the US, the influence of anthropogenic ignitions on fire physical characteristics, and how these patterns vary spatially

# This code:
# Imports fire data for the entire US and generates annual rasters (50km resolution) 
# All data are clipped to the CONUS and projected to UTM Zone 13 WGS1984
# samples those rasters at 50km

# Data associated with this code, in Data subfolder in this project folder:
# MTBS, MODIS active fire, Short data, States
# raw data Desktop/US_Pyromes/Data/Data

# Required packages:
library(rgdal)
library(sp)
library(raster)
library(sendmailR)
library(tidyr)
library(lubridate)
library(rgeos)
library(ggplot2)

setwd("/Users/meca3122/Dropbox/0_EarthLab/US_Pyromes")
setwd("/Users/megancattau 1/Dropbox/0_EarthLab/US_Pyromes")
getwd()

################################################################
####################  Table of Contents  ############################
################################################################

## 1. Import and process Data
## 2. Create annual rasters for each year, Rasterstack, and Sample @ grid


################################################################
###################  1. Import and Process Data  ####################
################################################################

# Import data - fire data, sampling grid, and states boundaries - all projectd to WGS 1984 UTM Zone 13N
MODIS<-readOGR("Data/MODIS","MODISaf_WGS8413N_c") 					# Spatial Points
MTBS<-readOGR("Data/MTBS","MTBS_WGS84_13N") 									# Spatial Polygons
Short <- read.table("Data/Short_Data/ShortDB.txt", sep=",",
               header=T, fill=T,
               quote = "\"", #quote key
               row.names = NULL, 
               stringsAsFactors = FALSE)          
States<-readOGR("Data/States","CONUS") 															# Spatial Polygons
Fishnet<- raster(ext=extent(States), resolution=50000)										# Raster that's extent of States and 50km resolution
writeRaster(Fishnet,"Data/Fishnet.grd", format="raster", overwrite=TRUE)

projection(Fishnet)<-crs(MODIS)

# MODIS preprocessing
names(MODIS)
MODIS$acq.date<-ymd(MODIS$ACQ_DATE)													# change to date
MODIS$JD<-yday(MODIS$acq.date)																	# add JD col
MODIS$FireYear<-year(MODIS$acq.date)															# add year column


#MTBS preprocessing
names(MTBS)
MTBS$running_JD<-ifelse(MTBS$StartMonth==1, 0,										# add JD col, doesn't account for leapyear
ifelse(MTBS$StartMonth==2, 31,
ifelse(MTBS$StartMonth==3, 59,
ifelse(MTBS$StartMonth==4, 90,
ifelse(MTBS$StartMonth==5, 120,
ifelse(MTBS$StartMonth==6, 151,
ifelse(MTBS$StartMonth==7, 181,
ifelse(MTBS$StartMonth==8, 212,
ifelse(MTBS$StartMonth==9, 243,
ifelse(MTBS$StartMonth==10, 273,
ifelse(MTBS$StartMonth==11, 304,
ifelse(MTBS$StartMonth==12, 334,0
))))))))))))
MTBS$JD<-MTBS$running_JD+MTBS$StartDay
MTBS$FireYear<-MTBS$Year																				# add year column

MTBS$FireID<-1:nrow(MTBS)																				# MTBS into points
MTBS_point1<-gCentroid(MTBS, byid=TRUE, id=MTBS$FireId)
MTBS_point<-SpatialPointsDataFrame(MTBS_point1,MTBS@data)


# Short preprocessing
names(Short)
#remove commas to delineate numeric intervals inserted by Arc - necessary preprocessing?
# Short$fid <- as.numeric(gsub(slim$fid, patt=",", replace=""))

#rounded log of fire size
Short$rlsize <- round(log10(Short$FIRE_SIZE), digits=1)
#size in hectares
Short$ha <- Short$FIRE_SIZE*0.4046856

#classing lightning vs. human ignitions
unique(Short$STAT_CAUSE_DESCR)
ig <- ifelse(Short$STAT_CAUSE_DESCR == "Lightning", "lightning", NA) #lightning
ig <- ifelse(Short$STAT_CAUSE_DESCR == "", "unknown", ig) #unknown class, exclude from analysis
ig <- ifelse(Short$STAT_CAUSE_DESCR == "Missing/Undefined", "unknown", ig) #unknown class, exclude from analysis
ig <- ifelse(Short$STAT_CAUSE_DESCR == "Miscellaneous", "human", ig) #human (Karen Short told Jennifer Balch that miscellaneous was mostly human)
ig <- ifelse(Short$STAT_CAUSE_DESCR == "Debris Burning", "human", ig) #human
ig <- ifelse(Short$STAT_CAUSE_DESCR == "Campfire", "human", ig) #human
ig <- ifelse(Short$STAT_CAUSE_DESCR == "Equipment Use", "human", ig) #human
ig <- ifelse(Short$STAT_CAUSE_DESCR == "Arson", "human", ig) #human
ig <- ifelse(Short$STAT_CAUSE_DESCR == "Children", "human", ig) #human
ig <- ifelse(Short$STAT_CAUSE_DESCR == "Railroad", "human", ig) #human
ig <- ifelse(Short$STAT_CAUSE_DESCR == "Smoking", "human", ig) #human 
ig <- ifelse(Short$STAT_CAUSE_DESCR == "Powerline", "human", ig) #human 
ig <- ifelse(Short$STAT_CAUSE_DESCR == "Structure", "human", ig) #human 
ig <- ifelse(Short$STAT_CAUSE_DESCR == "Fireworks", "human", ig) #human 

Short$ig <- ig
nrow(Short[Short$ig=="unknown",]) # 197906
nrow(Short[Short$ig=="human",]) # 1409412
nrow(Short[Short$ig=="lightning",]) # 273147


#clean up long/lat coordinates (x, y)
Short <- subset(Short, is.na(Short$LATITUDE) == FALSE) #cleaning up a few NAs
Short <- subset(Short, is.na(Short$LONGITUDE) == FALSE)
#remove unknown fires
Short <- subset(Short, Short$ig != "unknown")

Short$JD<-Short$DISCOVERY_DOY																			# add JD col, ddoy = discovery date
Short$FireYear<-Short$FIRE_YEAR																		# add year column

names(Short)
Short_points<-Short
coordinates(Short_points)<-~LONGITUDE + LATITUDE														# make it spatial

# Set the CRS and project it
proj4string(Short_points)<-CRS("+init=epsg:4269")
crs(Short_points)
Short_points_p <- spTransform(Short_points, proj4string(MTBS))
Short<-Short_points_p
Short_human<-Short_points_p[Short_points_p$ig=="human",]

aggregate(Short$ha, by=list(Category=Short$ig), FUN=sum)
# Human: 17,690,441 ha
# Lightning: 34,730,827 ha

aggregate(Short$FIRE_SIZE, by=list(Category=Short$ig), FUN=sum)
#   human 43714035
# lightning 85821753


###################################################################
###  2. Create annual rasters for each year, Rasterstack, and Sample @ grid  ###
###################################################################

######### parse all the fire layers by year into a list of vector objects, calling them "xxxxx_parsed", each vector in the list called xxxx_yyyy
parse_vector <- function(all_data, prefix, year_seq) {
	# The function takes a vector layer (arg: all_data; e.g., MTBS), parses it by year, gives each new object a sensible name with a defined prefix (arg: prefix; e.g., fire), 	
	# returns a list of vector objects
	# args= all_data (the original polygon / points data), prefix (prefix for the name of the parsed vector layers), year_seq (sequence of relevant years)

	# A nested function - parse the data based on Year
	separate_data_by_year<-function(year){
		all_data[all_data@data$FireYear==year,]
	}
	
	# Apply this function over the range of relevant years, resulting in a list of vector objects for each year
	vector_list <- lapply(year_seq, separate_data_by_year)
	# Name each object in the list prefix_year
	names(vector_list) <- paste(prefix, year_seq, sep = "_")
	
	vector_list
}


range(MODIS$FireYear)
MODIS_parsed<-parse_vector(MODIS, "MODIS", 2001:2016)

range(MTBS$FireYear)
MTBS_parsed<-parse_vector(MTBS, "MTBS", 1984:2015)
MTBS_point_parsed<-parse_vector(MTBS_point, "MTBS_point", 1984:2015)

range(Short$FireYear)
Short_parsed<-parse_vector(Short, "Short", 1992:2015)
Short_human_parsed<-parse_vector(Short_human, "Short_human", 1992:2015)



######### convert the a list of vector objects into the relevant rasters
# Generate annual rasters of:
# Number of fires - count of MODIS FRP points, MTBS, Short				
## 1.		Number_fires_MODIS
## 2.		Number_fires_MTBS
## 3.		Number_fires_Short

# Max, mean, std FRP, MODIS																				
## 4. 	Mean_FRP_MODIS
## 5.	 	Max_FRP_MODIS

# Mean, max fire even size MTBS, Short												
## 6. 	Mean_area_MTBS
## 7.	 	Max_area_MTBS
## 8. 	Mean_area_Short
## 9.		Max_area_Short

# Burned area MTBS, Short	
## 10.	Sum_area_MTBS
## 11.	Sum_area_Short

# Season Length (std * 2 JD), MODIS, MTBS, Short									
## 12.	Std_JD_MTBS
## 13.	Std_JD_MODIS
## 14.	Std_JD_Short

# Ignition type (Perc human ignitions) Short
## 15. 	Perc human ignitions Short

annual_rasters <- function(vector_data, template, prefix, year_seq, field, fun, background) {
	# The function takes a list of vector objects, creates annual rasters based on a particular variable, and gives each new raster a sensible name with a defined prefix (arg: prefix; e.g., fire), 		# args= vector_data (a list of vector objects), template (a raster whose extent / res will be used for the putput rasters), prefix (prefix for the name of the raster layers), year_seq (sequence of relevant years), field= field on which to based calculations (e.g., FRP), fun = function to apply on that field (e.g., max), background(what value to give raster if there are no vector elements in there)
		# returns a list of rasters
	
	# A nested function - rasterize a polygon
		raster_each_year<- function(polygon){
		# Set the resolution and extent based on a template raster
		r <- raster(ncol=ncol(template), nrow=nrow(template))
		extent(r)<-extent(template)
		projection(r)<-projection(template)
		
		if (length(polygon)==0)
			new_raster<-setValues(r, 0) 
		else
			new_raster<-rasterize(polygon, r, field=field, fun=fun, background=background)
			# new_raster[is.na(new_raster)]<-0 - remove this because we do want some to be NA
		new_raster
	}

	annual_rasters<-lapply(vector_data, raster_each_year)
	names(annual_rasters) <- paste(prefix, year_seq, sep = "_")
	
	annual_rasters
}

names(MODIS)
names(MTBS)
names(Short)


# args: (vector_data, template, prefix, year_seq, field, fun, background)


# Fire frequency
Number_fires_MODIS<-annual_rasters(MODIS_parsed, Fishnet, "MODIS_Numfires", 2001:2016, "FRP", fun="count", background=0) # for count, it doesn't matter what field you use if points
Number_fires_MTBS<-annual_rasters(MTBS_point_parsed, Fishnet, "MTBS_Numfires", 1984:2015, "FireID", fun="count", background=0) 
Number_fires_Short<-annual_rasters(Short_parsed, Fishnet, "Short_Numfires", 1992:2015, "FireYear", fun="count", background=0) # for count, it doesn't matter what field you use if points

# Fire intensity - give background of NA so that not included in means if no fire in them
Mean_FRP_MODIS<-annual_rasters(MODIS_parsed, Fishnet, "MODIS_meanFRP", 2001:2016, "FRP", fun=mean, background=NA)
Max_FRP_MODIS<-annual_rasters(MODIS_parsed, Fishnet, "MODIS_maxFRP", 2001:2016, "FRP", fun=max,  background=NA)

# Fire event size
Mean_area_MTBS<-annual_rasters(MTBS_point_parsed, Fishnet, "MTBS_meanArea", 1984:2015, "Acres", fun=mean, background=NA) 
Max_area_MTBS<-annual_rasters(MTBS_point_parsed, Fishnet, "MTBS_maxArea", 1984:2015, field="Acres", fun=max,  background=NA) 
Mean_area_Short<-annual_rasters(Short_parsed, Fishnet, "Short_meanArea", 1992:2015, "ha", fun=mean, background=NA) 
Max_area_Short<-annual_rasters(Short_parsed, Fishnet, "Short_maxArea", 1992:2015, "ha", fun=max, background=NA) 

# Burned area
Sum_area_MTBS<-annual_rasters(MTBS_point_parsed, Fishnet, "MTBS_sumArea", 1984:2015, "Acres", fun=sum, background=NA) 
Sum_area_Short<-annual_rasters(Short_parsed, Fishnet, "Short_sumArea", 1992:2015, "ha", fun=sum, background=NA) 

# Fire seasonality
Std_JD_MTBS<-annual_rasters(MTBS_point_parsed, Fishnet, "MTBS_stdJD", 1984:2015, "JD", fun=sd, background=NA)
Std_JD_MODIS<-annual_rasters(MODIS_parsed, Fishnet, "MODIS_stdJD", 2001:2016, "JD", fun=sd, background=NA) 
Std_JD_Short<-annual_rasters(Short_parsed, Fishnet, "Short_stdJD", 1992:2015, "JD", fun=sd, background=NA) 
# Make it Std JD * 2
Std2_JD_MTBS<-calc(stack(Std_JD_MTBS), function(x) x*2, forceapply=TRUE)
Std2_JD_MODIS<-calc(stack(Std_JD_MODIS), function(x) x*2, forceapply=TRUE)
Std2_JD_Short<-calc(stack(Std_JD_Short), function(x) x*2, forceapply=TRUE)

# Perc human ignitions
Number_fires_Short2<-annual_rasters(Short_parsed, Fishnet, "Short_Numfires", 1992:2015, "FireYear", fun="count", background=NA) # background is NA so that if there were no fires in that cell, the perc human ign value won't be calculated
Number_fires_Short_human<-annual_rasters(Short_human_parsed, Fishnet, "Short_Number_fires_human", 1992:2015, "FireYear", fun="count", background=0) #Background = 0 so that if there were fires in that cell but no human fires, the value is still calculated
Perc_fires_Short_human<-stack(Number_fires_Short_human) / stack(Number_fires_Short2)

results_rasterstack<-stack(stack(Number_fires_MODIS), stack(Number_fires_MTBS), stack(Number_fires_Short), stack(Mean_FRP_MODIS), stack(Max_FRP_MODIS), stack(Mean_area_MTBS), stack(Max_area_MTBS), stack(Mean_area_Short), stack(Max_area_Short), stack(Sum_area_MTBS), stack(Sum_area_Short), stack(Std2_JD_MTBS), stack(Std2_JD_MODIS), stack(Std2_JD_Short), stack(Perc_fires_Short_human))

writeRaster(results_rasterstack,"0_Anthro/Data/results_rasterstack.grd", format="raster", overwrite=TRUE)
# results_rasterstack<-stack("0_Anthro/Data/results_rasterstack.grd")							# Import sampled rasters - annual

# stats on each variable across all years rather than annual
Number_fires_MODIS_mean<-calc(results_rasterstack[[1:16]], mean)
Number_fires_MTBS_mean<-calc(results_rasterstack[[17:48]], mean)
Number_fires_Short_mean<-calc(results_rasterstack[[49:72]], mean)

Mean_FRP_MODIS_mean<-calc(results_rasterstack[[73:88]], mean, na.rm=TRUE)
Max_FRP_MODIS_mean<-calc(results_rasterstack[[89:104]], mean, na.rm=TRUE)

Mean_area_MTBS_mean<-calc(results_rasterstack[[105:136]], mean, na.rm=TRUE)
Max_area_MTBS_mean<-calc(results_rasterstack[[137:168]], mean, na.rm=TRUE)
Mean_area_Short_mean<-calc(results_rasterstack[[169:192]], mean, na.rm=TRUE)
Max_area_Short_mean<-calc(results_rasterstack[[193:216]], mean, na.rm=TRUE)

Sum_area_MTBS_mean<-calc(results_rasterstack[[217:248]], mean, na.rm=TRUE)
Sum_area_Short_mean<-calc(results_rasterstack[[249:272]], mean, na.rm=TRUE)

# (2 * SD for season length already calculated above)
Std_JD_MTBS_mean<-calc(results_rasterstack[[273:304]], mean, na.rm=TRUE)
Std_JD_MODIS_mean<-calc(results_rasterstack[[305:320]], mean, na.rm=TRUE)
Std_JD_Short_mean<-calc(results_rasterstack[[321:344]], mean, na.rm=TRUE)

Perc_fires_Short_human_mean<-calc(results_rasterstack[[345:368]], mean, na.rm=TRUE)


results_rasterstack_mean<-stack(Number_fires_MODIS_mean, Number_fires_MTBS_mean, Number_fires_Short_mean, Mean_FRP_MODIS_mean, Max_FRP_MODIS_mean, Mean_area_MTBS_mean, Max_area_MTBS_mean, Mean_area_Short_mean, Max_area_Short_mean,  Sum_area_MTBS_mean, Sum_area_Short_mean, Std_JD_MTBS_mean, Std_JD_MODIS_mean, Std_JD_Short_mean, Perc_fires_Short_human_mean)

writeRaster(results_rasterstack_mean,"0_Anthro/Data/results_rasterstack_mean.grd", format="raster", overwrite=TRUE)


######################### Get data into shape #############################

# if need to reimport:
# results_rasterstack<-stack("0_Anthro/Data/results_rasterstack.grd")							# Import sampled rasters - annual
# results_rasterstack_mean<-stack("0_Anthro/Data/results_rasterstack_mean.grd") # Import sampled rasters - mean
# States<-readOGR("Data/States","CONUS") 															# Import States layer

results_rasterstack_all<-stack(results_rasterstack_mean, results_rasterstack)			# Combine annual and mean sampled pyromes rasters
results_rasterstack_mask<-mask(results_rasterstack_all, States)									# mask combined sampled pyromes rasters

sample_points<-rasterToPoints(results_rasterstack_mask[[1]], spatial=TRUE)			# create sample point locations from one of the rasters
samples<-raster::extract(results_rasterstack_mask, sample_points, sp=TRUE)			# extract values at sample points
names(samples)
samples<-samples[-1]																											# remove repeat layer

# format for the annual layers is datasource_measuredthing_yyyy
# format for the mean layers is layer.x, so rewrite that below
names(samples)<-c("Number_fires_MODIS_mean", "Number_fires_MTBS_mean", "Number_fires_Short_mean", "Mean_FRP_MODIS_mean", "Max_FRP_MODIS_mean", "Mean_area_MTBS_mean", "Max_area_MTBS_mean", "Mean_area_Short_mean", "Max_area_Short_mean", "Sum_area_MTBS_mean", "Sum_area_Short_mean", "Std_JD_MTBS_mean", "Std_JD_MODIS_mean", "Std_JD_Short_mean", "Perc_human_Short_mean", names(results_rasterstack))

samples_p <- SpatialPointsDataFrame(samples, data=samples@data, proj4string=crs(MODIS))		# project 
samples_p$FID<-1:nrow(samples_p)																										# put FID in there


samples_df<-data.frame(samples_p)
samples_df<-samples_df[,-387]																				# remove 'optional' column


# Data/Ecoregion_state/Eco_L1_pclp Ecoregions projected to match States projection and clipped to States extent
Ecoregion<-readOGR("Data/Ecoregion_state", "Eco_L1_pclp")	
compareCRS(Ecoregion, samples_p)											  						#TRUE
proj4string(Ecoregion)<-crs(samples_p)
proj4string(Ecoregion)<-crs(samples_p)
overlay <- fortify(Ecoregion, region="NA_L1NAME")
write.csv(overlay, "0_Anthro/Data/overlay.csv")
overlay<-read.csv("0_Anthro/Data/overlay.csv")

eco_data<-sp::over(samples_p, Ecoregion[,"NA_L1NAME"])
samples_df$ecoregion<-eco_data$NA_L1NAME
samples_p$ecoregion<-eco_data


### Write and retrieve samples_df dataframe
write.csv(samples_df, "0_Anthro/Data/samples_df.csv")
# samples_df<-read.csv("0_Anthro/Data/samples_df.csv")
# names(samples_df)
# samples_df<-samples_df[,-1]	









