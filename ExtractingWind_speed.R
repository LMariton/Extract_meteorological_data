
##############################################################################################################
# This script allows to extract the wind speed (in m.s-1) at sites for given nights at 18, 0 or 6 o'clock.
# Those sites have to be in metropolitan France (and Corsica)
#
# It returns the data frame of each night and site (given in arguments) with 3 columns for the wind speed
# for the day studied (wind_night_18 for the wind speed at 18 o'clock, i.e beginning of night, wind_night_00
# for the wind speed at midnight, wind_night_06 for the wind speed at 6 o'clock, i.e. end of the night)
# and as many columns as the number of days before the survey specified in arguments x 3 
# (e.g. wind speed of the day before at o'clock: wind_day_18_before1, wind speed two days before at 18 
# o'clock: wind_day_18_before2, etc).
#
##############################################################################################################
#
# /!\ The data used have a low spatial resolution (2°) but a high temporal one (every 6 hours)
# Package RNCEP - data from the NOAA
#
# An internet connection is required.
#
##############################################################################################################
#
# Arguments : 
# 
# tableNightsSites : a table with (minimum):
# -> a column named "Latitude" with the site latitude in WGS82
# -> a column named "Longitude" with the site longitude in WGS82
# -> a column with the sites unique IDs
# -> a column named "Date" with the date of the survey for each sites.
#       The date for a night should be the date when it began.
#       Dates should be formatted this way : "2020-05-26" if characters.
#       They might also be c("POSIXct","POSIXt") objects.
#
# Col_ID : column name of the table containing the site unique IDs (ex : "CODE_SITE")
#
# nbrPreviousDays : number of previous nights for which T° should be returned
# (for example "3", to obtain the temperature for the survey night but for also the three previous ones)
#
# year_beg : year when the surveys began
#
# year_end : year when the surveys ended
#
##############################################################################################################

extract_WindSpeed <- function(tableNightsSites,Col_ID,nbrPreviousDays,year_beg,year_end){

  #Install and open required packages
  load <- function(pkg){
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg))
      install.packages(new.pkg, dependencies = TRUE)
    sapply(pkg, require, character.only = TRUE)
  } 
  packages <- c("RNCEP","sf","lubridate","data.table","stringr","dplyr")
  load(packages)
  
  #Data formatting
  if (all(class(tableNightsSites[,which(colnames(tableNightsSites)=="Date")])!=c("POSIXct","POSIXt"))){
    tableNightsSites[,which(colnames(tableNightsSites)=="Date")] <- parse_date_time(
      as.character(tableNightsSites[,which(colnames(tableNightsSites)=="Date")]),"%Y-%m-%d", tz="Europe/Paris")
  }
  
  #####I. Extraction####
  Years=c(year_beg,year_end)
  LatWindow=c(41,52)
  LongWindow=c(-6,9)
  
  point_wind_U=NCEP.gather("uwnd.sig995","surface"
                           ,months.minmax=c(1,12)
                           ,years.minmax=Years
                           ,lat.southnorth=LatWindow, lon.westeast=LongWindow)
  
  
  point_wind_V=NCEP.gather("vwnd.sig995","surface"
                           ,months.minmax=c(1,12)
                           ,years.minmax=Years
                           ,lat.southnorth=LatWindow, lon.westeast=LongWindow)
  
  point_wind_UDT <- as.data.frame(as.data.table(point_wind_U))
  point_wind_VDT <- as.data.frame(as.data.table(point_wind_V))
  all(point_wind_UDT$V1==point_wind_VDT$V1)
  all(point_wind_UDT$V2==point_wind_VDT$V2)
  all(point_wind_UDT$V3==point_wind_VDT$V3)
  point_wind <- point_wind_UDT
  colnames(point_wind) <- c("latitude","longitude","year_date_time","wind_ms_U")
  point_wind$wind_ms_V <- point_wind_VDT$value
  point_wind$wind_ms <- sqrt((point_wind$wind_ms_U^2 + point_wind$wind_ms_V^2))
  
  point_wind$coord <- paste0(point_wind$latitude,"_",point_wind$longitude)
  point_wind$hour <- as.numeric(str_sub(point_wind$year_date_time,start=-2))
  point_wind$date <- parse_date_time(str_sub(point_wind$year_date_time,end=-4),"%Y_%m_%d", tz="Europe/Paris")
  
  ####II. Correspondance between survey points and RNCEP points####
  
  site_Nights <- distinct(tableNightsSites[,which(colnames(tableNightsSites) %in% c(Col_ID,"Longitude","Latitude"))])
  
  #Creation of a survey site layer
  sites_SHP <- st_as_sf(site_Nights, coords = c("Longitude","Latitude"), crs = 4326)
  sites_SHP <- st_transform(sites_SHP,crs=2154)
  
  #Creation of RNCEP points layer
  site_RNCEP <- distinct(point_wind[,which(colnames(point_wind) %in% c("coord","longitude","latitude"))])
  RNCEP_SHP <- st_as_sf(site_RNCEP, coords = c("longitude","latitude"), crs = 4326)
  RNCEP_SHP <- st_transform(RNCEP_SHP,crs=2154)

  #Correspondance
  site_Nights$corresp_RNCEP <- site_RNCEP$coord[st_nearest_feature(sites_SHP,RNCEP_SHP)]
  tableNightsSites$corresp_RNCEP <- site_Nights$corresp_RNCEP[match(tableNightsSites[,which(colnames(tableNightsSites)==Col_ID)],site_Nights[,which(colnames(site_Nights)==Col_ID)])]
  
  ####III. Wind speed extraction####
  
  tableNightsSites$wind_night_18 <- rep(NA,dim(tableNightsSites)[1])
  tableNightsSites$wind_night_00 <- rep(NA,dim(tableNightsSites)[1])
  tableNightsSites$wind_night_06 <- rep(NA,dim(tableNightsSites)[1])
  
    for (j in (1:nbrPreviousDays)){
      for (k in c("18","00","06")){
        tableNightsSites[,dim(tableNightsSites)[2]+1] <- rep(NA,dim(tableNightsSites)[1])
        colnames(tableNightsSites)[dim(tableNightsSites)[2]] <- paste0("wind_night_",k,"_before",j)
      }
    }

    for (index in c(1:dim(tableNightsSites)[1])){
      
      #only data for this site
      point_wind_subset <- point_wind[which(point_wind$coord==tableNightsSites$corresp_RNCEP[index]),]
      
      if(length(which((point_wind_subset$date==tableNightsSites$Date[index])&(point_wind_subset$hour==18)))>0){
        tableNightsSites$wind_night_18[index] <- point_wind_subset$wind_ms[
          which((point_wind_subset$date==tableNightsSites$Date[index])&(point_wind_subset$hour==18))]
        }
      
      if(length(which((point_wind_subset$date==tableNightsSites$Date[index]+days(1))&(point_wind_subset$hour==00)))>0){
        tableNightsSites$wind_night_00[index] <- point_wind_subset$wind_ms[
          which((point_wind_subset$date==tableNightsSites$Date[index]+days(1))&(point_wind_subset$hour==00))]
        }
      if(length(which((point_wind_subset$date==tableNightsSites$Date[index]+days(1))&(point_wind_subset$hour==06)))>0){
        tableNightsSites$wind_night_06[index] <- point_wind_subset$wind_ms[
          which((point_wind_subset$date==tableNightsSites$Date[index]+days(1))&(point_wind_subset$hour==06))]
      }
      
      
      if(nbrPreviousDays!=0){
        
        for (nbrNight in (1:nbrPreviousDays)){
          
          if (length(which((point_wind_subset$date==tableNightsSites$Date[index]-days(nbrNight))&(point_wind_subset$hour==18)))>0){
            tableNightsSites[
              index,which(colnames(tableNightsSites)==paste0(
                "wind_night_18_before",nbrNight))] <- point_wind_subset$wind_ms[
                  which((point_wind_subset$date==tableNightsSites$Date[index]-days(nbrNight))&(point_wind_subset$hour==18))]
          }
          
          if (length(which((point_wind_subset$date==tableNightsSites$Date[index]-days(nbrNight-1))&(point_wind_subset$hour==00)))>0){
          tableNightsSites[
            index,which(colnames(tableNightsSites)==paste0(
              "wind_night_00_before",nbrNight))] <- point_wind_subset$wind_ms[
                which((point_wind_subset$date==tableNightsSites$Date[index]-days(nbrNight-1))&(point_wind_subset$hour==00))]
          }
          
          if (length(which((point_wind_subset$date==tableNightsSites$Date[index]-days(nbrNight-1))&(point_wind_subset$hour==06)))>0){
            tableNightsSites[
              index,which(colnames(tableNightsSites)==paste0(
                "wind_night_06_before",nbrNight))] <- point_wind_subset$wind_ms[
                  which((point_wind_subset$date==tableNightsSites$Date[index]-days(nbrNight-1))&(point_wind_subset$hour==06))]
          }
        }
      }
      
      cat(paste0(round(index/dim(tableNightsSites)[1]*100,digits=1),"%\r"))
      
    }
  
  tableNightsSites <- tableNightsSites[,which(colnames(tableNightsSites)!="corresp_RNCEP")]
  return(tableNightsSites)
}
