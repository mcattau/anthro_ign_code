
# Anthropogenic and lightning-started fires are becoming larger and more frequent over a longer season length in the U.S.
# Global Ecology and Biogeography
# Code author: Megan Cattau
# Contact info: megan.cattau@gmail.com or megancattau@boisestate.edu
# Manuscript authors: Megan Cattau, Carol Wessman, Adam Mahood, Jennifer Balch
# Project description: This project aims to elucidate how fire physical characteristics are changing over time in the US, the influence of anthropogenic ignitions on fire physical characteristics, and how these patterns vary over space and time.

# This code:
# Imports fire data for the entire US 
# generates annual rasters (50km resolution) 
# samples those rasters at points 
# All data are clipped to the CONUS and projected to UTM Zone 13 WGS1984


# Required packages:
library(rgdal)
library(sp)
library(sf)
library(dplyr)
library(raster)
library(sendmailR)
library(tidyr)
library(lubridate)
library(rgeos)
library(ggplot2)
library(assertthat)


################################################################
####################  Table of Contents  ############################
################################################################

## 1. Import and process Data
## 2. Create annual rasters for each year, Rasterstack, and Sample @ grid


################################################################
###################  1. Import and Process Data  ####################
################################################################

# Data will include US States boundaries and three fire data sources: MTBS fire perimeters, MODIS active fire data, and the FPA-FOD dataset.
# create a Data folder in your working directory 


setwd("/Users/megancattau 1/Dropbox/0_EarthLab/US_Pyromes")
setwd("Data/")

# Projection for layers
#EPSG:32613
data_crs<- " +proj=utm +zone=13 +datum=WGS84 +units=m +no_defs +ellps=WGS84 "


### US States
States_download <- file.path('States', 'cb_2016_us_state_20m.shp')
if (!file.exists(States_download)) {
	from <- "https://www2.census.gov/geo/tiger/GENZ2016/shp/cb_2016_us_state_20m.zip"
	to <- paste0('States', ".zip")
	download.file(from, to)
	unzip(to, exdir = 'States')
	unlink(to)
	assert_that(file.exists(States_download))
}   
               
States <- st_read(dsn = 'States', layer = "cb_2016_us_state_20m", quiet = TRUE) %>%
filter(!(NAME %in% c("Alaska", "Hawaii", "Puerto Rico"))) %>%
dplyr::select(STATEFP, STUSPS) %>%
setNames(tolower(names(.))) %>% 
st_transform(., data_crs)
          

# Create a raster that's extent of States and 50km resolution and write into Data folder
Fishnet<- raster(ext=extent(States), resolution=50000)		
projection(Fishnet)<-crs(data_crs)	
writeRaster(Fishnet,"Fishnet.grd", format="raster", overwrite=TRUE)


### MTBS fire perimeters 
# will be downloaded directly to Data folder
MTBS_download <- file.path('MTBS', 'mtbs_perims_DD.shp')
if (!file.exists(MTBS_download)) {
	from <- "https://edcintl.cr.usgs.gov/downloads/sciweb1/shared/MTBS_Fire/data/composite_data/burned_area_extent_shapefile/mtbs_perimeter_data.zip"
	to <- paste0('MTBS', ".zip")
	download.file(from, to)
	unzip(to, exdir = 'MTBS')
	unlink(to)
	assert_that(file.exists(MTBS_download))
}

MTBS <- st_read(dsn = 'MTBS', layer = "mtbs_perims_DD", quiet = TRUE) %>%
st_transform(., data_crs)

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
MTBS$FireID<-1:nrow(MTBS)																				


### MODIS active fire data
# Submit a request via the link below for Country->United States from Jan 2003 - Dec 2016, MODIS C6, and download into Data folder and unzip
# https://firms.modaps.eosdis.nasa.gov/download/

MODIS <- st_read(dsn = 'DL_FIRE_M6_88165', layer = "fire_archive_M6_88165", quiet = TRUE) %>%
st_transform(., data_crs)

# MODIS preprocessing
names(MODIS)
MODIS$acq.date<-ymd(MODIS$ACQ_DATE)													# change to date
MODIS$JD<-yday(MODIS$acq.date)																	# add JD col
MODIS$FireYear<-year(MODIS$acq.date)															# add year column
MODIS<-MODIS[MODIS$FireYear>=2003,]														# Remove all before 2003


 ### FPA-FOD
 # will be downloaded directly to Data folder
 # Can be directly downloaded from here https://www.fs.usda.gov/rds/archive/catalog/RDS-2013-0009.4 as a Geodatabase: https://www.fs.usda.gov/rds/archive/products/RDS-2013-0009.4/RDS-2013-0009.4_GDB.zip
 # The below has been processed to be easier to work with:
 https://www.dropbox.com/s/xiyzlme6rt3ierx/shortfireseco073015.txt?dl=0
   
   
Short <- read.table("Short_Data/ShortDB.txt", sep=",",
               header=T, fill=T,
               quote = "\"", #quote key
               row.names = NULL, 
               stringsAsFactors = FALSE)      
  
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
Short_points_p <- spTransform(Short_points, (data_crs))
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
		all_data[all_data$FireYear==year,]
	}
	
	# Apply this function over the range of relevant years, resulting in a list of vector objects for each year
	vector_list <- lapply(year_seq, separate_data_by_year)
	# Name each object in the list prefix_year
	names(vector_list) <- paste(prefix, year_seq, sep = "_")
	
	vector_list
}

range(MODIS$FireYear)
MODIS_parsed<-parse_vector(MODIS, "MODIS", 2003:2016)

range(MTBS$FireYear)
MTBS_parsed<-parse_vector(MTBS, "MTBS", 1984:2015)
# MTBS_point_parsed<-parse_vector(MTBS_point, "MTBS_point", 1984:2015)

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
Number_fires_MODIS<-annual_rasters(MODIS_parsed, Fishnet, "MODIS_Numfires", 2003:2016, "FRP", fun="count", background=0) # for count, it doesn't matter what field you use if points
Number_fires_MTBS<-annual_rasters(MTBS_parsed, Fishnet, "MTBS_Numfires", 1984:2015, "FireID", fun="count", background=0) 
Number_fires_Short<-annual_rasters(Short_parsed, Fishnet, "Short_Numfires", 1992:2015, "FireYear", fun="count", background=0) # for count, it doesn't matter what field you use if points

# Fire intensity - give background of NA so that not included in means if no fire in them
Mean_FRP_MODIS<-annual_rasters(MODIS_parsed, Fishnet, "MODIS_meanFRP", 2003:2016, "FRP", fun=mean, background=NA)
Max_FRP_MODIS<-annual_rasters(MODIS_parsed, Fishnet, "MODIS_maxFRP", 2003:2016, "FRP", fun=max,  background=NA)

# Fire event size
Mean_area_MTBS<-annual_rasters(MTBS_parsed, Fishnet, "MTBS_meanArea", 1984:2015, "Acres", fun=mean, background=NA) 
Max_area_MTBS<-annual_rasters(MTBS_parsed, Fishnet, "MTBS_maxArea", 1984:2015, field="Acres", fun=max,  background=NA) 
Mean_area_Short<-annual_rasters(Short_parsed, Fishnet, "Short_meanArea", 1992:2015, "ha", fun=mean, background=NA) 
Max_area_Short<-annual_rasters(Short_parsed, Fishnet, "Short_maxArea", 1992:2015, "ha", fun=max, background=NA) 

# Burned area
Sum_area_MTBS<-annual_rasters(MTBS_parsed, Fishnet, "MTBS_sumArea", 1984:2015, "Acres", fun=sum, background=NA) 
Sum_area_Short<-annual_rasters(Short_parsed, Fishnet, "Short_sumArea", 1992:2015, "ha", fun=sum, background=NA) 

# Fire seasonality
Std_JD_MTBS<-annual_rasters(MTBS_parsed, Fishnet, "MTBS_stdJD", 1984:2015, "JD", fun=sd, background=NA)
Std_JD_MODIS<-annual_rasters(MODIS_parsed, Fishnet, "MODIS_stdJD", 2003:2016, "JD", fun=sd, background=NA) 
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

writeRaster(results_rasterstack,"results_rasterstack.grd", format="raster", overwrite=TRUE)
# results_rasterstack<-stack("results_rasterstack.grd")							# Import sampled rasters - annual

# stats on each variable across all years rather than annual
Number_fires_MODIS_mean<-calc(results_rasterstack[[1:14]], mean)
Number_fires_MTBS_mean<-calc(results_rasterstack[[15:46]], mean)
Number_fires_Short_mean<-calc(results_rasterstack[[47:70]], mean)

Mean_FRP_MODIS_mean<-calc(results_rasterstack[[71:84]], mean, na.rm=TRUE)
Max_FRP_MODIS_mean<-calc(results_rasterstack[[85:98]], mean, na.rm=TRUE)

Mean_area_MTBS_mean<-calc(results_rasterstack[[99:130]], mean, na.rm=TRUE)
Max_area_MTBS_mean<-calc(results_rasterstack[[131:162]], mean, na.rm=TRUE)
Mean_area_Short_mean<-calc(results_rasterstack[[163:186]], mean, na.rm=TRUE)
Max_area_Short_mean<-calc(results_rasterstack[[187:210]], mean, na.rm=TRUE)

Sum_area_MTBS_mean<-calc(results_rasterstack[[211:246]], mean, na.rm=TRUE)
Sum_area_Short_mean<-calc(results_rasterstack[[243:266]], mean, na.rm=TRUE)

# (2 * SD for season length already calculated above)
Std_JD_MTBS_mean<-calc(results_rasterstack[[267:298]], mean, na.rm=TRUE)
Std_JD_MODIS_mean<-calc(results_rasterstack[[299:312]], mean, na.rm=TRUE)
Std_JD_Short_mean<-calc(results_rasterstack[[313:336]], mean, na.rm=TRUE)

Perc_fires_Short_human_mean<-calc(results_rasterstack[[337:360]], mean, na.rm=TRUE)


results_rasterstack_mean<-stack(Number_fires_MODIS_mean, Number_fires_MTBS_mean, Number_fires_Short_mean, Mean_FRP_MODIS_mean, Max_FRP_MODIS_mean, Mean_area_MTBS_mean, Max_area_MTBS_mean, Mean_area_Short_mean, Max_area_Short_mean,  Sum_area_MTBS_mean, Sum_area_Short_mean, Std_JD_MTBS_mean, Std_JD_MODIS_mean, Std_JD_Short_mean, Perc_fires_Short_human_mean)


writeRaster(results_rasterstack_mean,"results_rasterstack_mean.grd", format="raster", overwrite=TRUE)


######################### Get data into shape #############################

# if need to reimport:
# results_rasterstack<-stack("results_rasterstack.grd")							# Import sampled rasters - annual
# results_rasterstack_mean<-stack("results_rasterstack_mean.grd") # Import sampled rasters - mean
# States<-readOGR("States","CONUS") 															# Import States layer

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
samples_df<-samples_df[,-379]																				# remove 'optional' column


### EPA Level I Ecoregions 
# Available here: https://www.epa.gov/eco-research/ecoregions-north-america
# will be downloaded directly to Data folder
Ecoregion_download <- file.path('Ecoregion', 'NA_CEC_Eco_Level1.shp')
if (!file.exists(Ecoregion_download)) {
	from <- "ftp://newftp.epa.gov/EPADataCommons/ORD/Ecoregions/cec_na/na_cec_eco_l1.zip"
	to <- paste0('Ecoregion', ".zip")
	download.file(from, to)
	unzip(to, exdir = 'Ecoregion')
	unlink(to)
	assert_that(file.exists(Ecoregion_download))
}

Ecoregion <- st_read(dsn = 'Ecoregion', layer = "NA_CEC_Eco_Level1", quiet = TRUE) %>%
st_transform(., data_crs)

Ecoregion2<-as(Ecoregion, 'Spatial')

eco_data<-sp::over(samples_p, Ecoregion2[,"NA_L1NAME"])
samples_df$ecoregion<-eco_data$NA_L1NAME
samples_p$ecoregion<-eco_data


### Write and retrieve samples_df dataframe
write.csv(samples_df, "samples_df.csv")
# samples_df<-read.csv("0_Anthro/Data/samples_df.csv")
# names(samples_df)
# samples_df<-samples_df[,-1]	









