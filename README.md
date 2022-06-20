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
