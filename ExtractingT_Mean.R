##############################################################################################################
# This script allows to extract the raw mean temperature and the mean anomaly (in °C) at sites for given dates.
# Those sites have to be in metropolitan France (and Corsica)
#
# The anomaly is defined as the difference between the mean temperature of the site at the given date and 
# the mean of the mean temperature at this site for the given day and month from 1980 to 2010.
#
# It returns the data frame of each day and site (given in arguments) with a new column for the T° of 
# the day (tg_day), and as many columns as the number of days before the survey specified in arguments 
# (e.g. T° of the day before: tg_day_before1, T° two days before: tg_day_before2, etc) and the same 
# things for anomalies (tg_day_ano, tg_day_ano_before1, tg_day_ano_before2, etc.).
#
##############################################################################################################
#
# Meteorological data (from E-OBS) as to be download at this link : 
# https://surfobs.climate.copernicus.eu/dataaccess/access_eobs.php#datafiles
# Choose either  : Ensemble Mean TG from the 0.25 deg. regular grid 
# or the 0.1 deg. regular grid (more precise but heavier)
#
##############################################################################################################
#
# Arguments : 
#
# path_to_meteo_nc : path to access to the meteo .nc file
# 
# tableDaysSites : a table with (minimum):
# -> a column named "Latitude" with the site latitude in WGS82
# -> a column named "Longitude" with the site longitude in WGS82
# -> a column with the sites unique IDs
# -> a column named "Date" with the date of the survey for each sites
#       Dates should be formatted this way : "2020-05-26" if characters
#       They might also be c("POSIXct","POSIXt") objects
#
# nbrPreviousDays : number of previous days for which T° should be returned
# (for example "3", to obtain the temperature for the survey date but for also the three previous days)

##############################################################################################################

extract_t_mean <- function(path_to_meteo_nc,tableDaysSites,nbrPreviousDays){

  ####Preliminary####
  
  #Install and open required packages
    load <- function(pkg){
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg))
      install.packages(new.pkg, dependencies = TRUE)
    sapply(pkg, require, character.only = TRUE)
  } 
  packages <- c("ncdf4","sf","stringr","lubridate","ggplot2","tidyr","data.table")
  load(packages)
  
  #creation of a site layer
  sites_SHP <- st_as_sf(tableDaysSites, coords = c("Longitude","Latitude"), crs = 4326)
  sites_SHP <- st_transform(sites_SHP,crs=2154)
  
  if (all(class(tableDaysSites[,which(colnames(tableDaysSites)=="Date")])!=c("POSIXct","POSIXt"))){
    tableDaysSites[,which(colnames(tableDaysSites)=="Date")] <- parse_date_time(
      as.character(tableDaysSites[,which(colnames(tableDaysSites)=="Date")]),"%Y-%m-%d", tz="Europe/Paris")
  }
  
  ####I. Meteo data formatting####
  
  #Open meteo file
  tg.ncdf <- nc_open(path_to_meteo_nc)
  
  #This file has three dimensions : longitude, latitude, date
  #Extraction of the values taken by each of these dimensions
  lon <- tg.ncdf$dim$longitude$vals
  lat <- tg.ncdf$dim$latitude$vals
  date <- tg.ncdf$dim$time$vals
  date <- as_date(date,origin="1950-01-01") #conversion in dates object
  
  #lim date : from 1980-01-01 to today
  start_date = which(date=="1980-01-01")
  end_date = length(date)
  len_date <- end_date-start_date+1
  
  #lim longitude (for France)
  start_lon <- which(lon+6==min(abs(lon+6)))
  stop_lon <- which(lon-10==min(abs(lon-10)))
  len_lon <- stop_lon-start_lon+1
  
  #lim longitude (for France)
  start_lat <- which(lat-41==min(abs(lat-41)))
  stop_lat <- which(lat-52==min(abs(lat-52)))
  len_lat <- stop_lat-start_lat+1
  
  #data extraction from the .nc file into an array
  tg_array <- ncvar_get(tg.ncdf,"tg",start=c(start_lon,start_lat,start_date),count=c(len_lon,len_lat,len_date))
  dimnames(tg_array) <- list(c(lon[start_lon:stop_lon]),c(lat[start_lat:stop_lat]),c(as.character(date[start_date:end_date])))
  
  print("1. Data formatting OK")
  
  ####II. Link between sites and meteo points####
  
  ####____A. Detect longitude * latitude for which data are available####
    
  #Take one date not to recent (data has to exist) but not old either
  #Here arbitrarily : 2019-05-15
  #And extract temperature for each longitude * latitude at this date (array with 2 dimensions)
    
  date_map = which(date=="2019-05-15")
  tg_link_array <- ncvar_get(tg.ncdf,"tg",start=c(start_lon,start_lat,date_map),count=c(len_lon,len_lat,1)) #extraction
  dimnames(tg_link_array) <- list(c(lon[start_lon:stop_lon]),c(lat[start_lat:stop_lat])) #longitude and latitude coord
  tg_link_arrayDT <- as.data.frame(tg_link_array) #array to data.frame
  tg_link_arrayDT$lon <- rownames(tg_link_arrayDT) #longitude in a column
  tg_link_arrayDT <- gather(tg_link_arrayDT,key="lat",value="tg",c(1:((dim(tg_link_arrayDT)[2])-1))) #reshape the data frame in a data frame with a column longitude, a column latitude and a column with the mean temperature
  tg_link_arrayDT <- tg_link_arrayDT[which(!(is.na(tg_link_arrayDT$tg))),] #remove longitude * latitude with NA
  tg_link_arrayDT$coord <- paste0(tg_link_arrayDT$lon,"_",tg_link_arrayDT$lat) #create a column with the concatenation of longitude and latitude
  
  ####____B. Create a shapefile with the points of the meteo data####
  
  coord_meteo_SHP <- st_as_sf(tg_link_arrayDT, coords = c("lon","lat"), crs = 4326)
  coord_meteo_SHP <- st_transform(coord_meteo_SHP,crs=2154)
    
  ####____C. create the link between meteo data and sites####
    
  #Creation of a new column in the shapefile with the sites that indicate, for each site, the closest point
  #of the meteo data
    
  sites_SHP$corresp_EObs <- coord_meteo_SHP$coord[st_nearest_feature(sites_SHP,coord_meteo_SHP)]
   
  print("2. Link between sites and meteo points OK")
  
  ####III. Temperature and anomaly extraction for each night####
  
  tableDaysSites$corresp_EObs <- as.data.frame(sites_SHP$corresp_EObs)
  
  tableDaysSites$tg_day <- rep(NA,dim(tableDaysSites)[1])
  
  if (nbrPreviousDays !=0){
    for (i in 1:nbrPreviousDays){
      tableDaysSites[,dim(tableDaysSites)[2]+1] <- rep(NA,dim(tableDaysSites)[1])
      colnames(tableDaysSites)[dim(tableDaysSites)[2]] <- paste0("tg_day_before",i)
    }
  }
  
  tableDaysSites$tg_day_ano <- rep(NA,dim(tableDaysSites)[1])
  
  if (nbrPreviousDays !=0){
    for (i in 1:nbrPreviousDays){
      tableDaysSites[,dim(tableDaysSites)[2]+1] <- rep(NA,dim(tableDaysSites)[1])
      colnames(tableDaysSites)[dim(tableDaysSites)[2]] <- paste0("tg_day_ano_before",i)
    }
  }
  
  for (index in c(1:dim(tableDaysSites)[1])){
    
    date_day <- tableDaysSites[index,which(colnames(tableDaysSites)=="Date")]
    
    if (as.character(date_day) %in% (dimnames(tg_array)[[3]])==T){
      
      coord_meteo <- sites_SHP$corresp_EObs[index]
      coord_lon <- str_split(coord_meteo,"_",simplify = T)[1]
      coord_lat <- str_split(coord_meteo,"_",simplify = T)[2]
      
      #raw temperature
      tableDaysSites$tg_day[index] <- tg_array[coord_lon,coord_lat,as.character(date_day)]
      
      #anomaly
      date_1980_2010 <- paste0(c(1980:2010),str_sub(date_day,start=5))
      anomaly_1980_2010 <- mean(tg_array[coord_lon,coord_lat,date_1980_2010])
      tableDaysSites$tg_day_ano[index] <- tableDaysSites$tg_day[index]-anomaly_1980_2010
      
      if (nbrPreviousDays!=0){
        for (i in (1:nbrPreviousDays)){
          tableDaysSites[index,which(colnames(tableDaysSites)==paste0("tg_day_before",i))] <- tg_array[coord_lon,coord_lat,as.character(date_day-days(i))]
          tableDaysSites[index,which(colnames(tableDaysSites)==paste0("tg_day_ano_before",i))] <- tableDaysSites[index,which(colnames(tableDaysSites)==paste0("tg_day_before",i))]-anomaly_1980_2010
        }
      }
      
    }
    
    cat(paste0(round(index/(dim(tableDaysSites)[1])*100,digits=1),"%\r"))
  }
  
  print("3. Temperature and anomaly extraction for each night OK")
  
  return(tableDaysSites)
}