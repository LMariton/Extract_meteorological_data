# Extract_meteorological_data
Scripts to extract meteorological data

### Extract_eobs.R
This script extracts at sites in Europe for given dates:
- raw mean temperature and mean anomaly (in °C)
- or precipitations sum (in mm)
- or mean wind speed (in m.s-1)

It uses the E-Obs meteorological data: https://surfobs.climate.copernicus.eu/dataaccess/access_eobs.php#datafiles

### ExtractingWind_speed.R
This script extracts the wind speed (in m.s-1) at sites for given nights at 18, 0 or 6 o'clock (minor changes can allow to transform this script for diurnal surveys).
Those sites have to be in metropolitan France (and Corsica) (if needed, minor changes will allow to use this script for other sites throughout the word).
/!\ The data used have a low spatial resolution (2°) but a high temporal one (every 6 hours).
They are extracted thanks to the package RNCEP (an internet connection is required), data come from from the NOAA data base.

### ExtractingCloud_cover.R
This script extracts the cloud cover (in %) at sites for given nights at 18, 0 or 6 o'clock (minor changes can allow to transform this script for diurnal surveys).
Those sites have to be in metropolitan France (and Corsica) (if needed, minor changes will allow to use this script for other sites throughout the word).
/!\ The data used have a low spatial resolution (2°) but a high temporal one (every 6 hours).
They are extracted thanks to the package RNCEP (an internet connection is required), data come from from the NOAA data base.


### Warning:
**These scripts are « functions »**, to use them you have to download them and keep the R. files somewhere on your computer.  

In theory, you don’t even need to open them with R, however I strongly advise you to do so to read the first lines presenting what the function does and how to use it (in particular what the arguments of the function are and how your dataset should be structured). 

In your script, simply add the following command:  
`source(“pathToTheRFile”)`  
Then you can use the function!  

Here is an example of what your script should look like: 

`path <- "D:/Data/Data_meteo/tg_ens_mean_0.1deg_reg_v23.1e.nc"`  
`surveysData <- read.csv("C:/Users/Name/Desktop/Surveys.csv")`  
`nbr <- 2`

`source("C:/Users/Name/Document/Extract_eobs.R")`  
`surveysDataNew <- extract_eobs(path,surveysData,nbr)`

