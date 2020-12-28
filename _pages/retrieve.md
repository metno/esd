---
layout: default 
title: Retrieve
parent: Utilities
nav_order: 2
permanent: \retrieve\ 
---

# Retrieving data: I/O
There are different types of data that can be handled by the ‘esd’ tool: station data, gridded data, and (storm) trajectories. Station data contain meteorological observations recorded at weather (or hydrological) stations, while gridded data can comprise various analyses (e.g. E-OBS gridded version of the European and Climate Assessment data set), reanalyses (e.g. NCEP, ERA, ...) or global/regional climate model results (e.g. CMIP3/5 experiment). Trajectories are
mainly used for analysis of storm tracks (e.g. IMILAST13). There are two main methods for retrieving data in the ‘esd’ tool: ‘station’ and ‘retrieve’.
It is also possible to read data using the R-built in functions and convert it to esd data format. This requires more effort from the user using the ‘esd’ pre-defined functions `as.station`, `as.field`), and `as.trajectory`. The package comes with a set of sample data mainly included for demonstration purposes,
testing, or troubleshooting. For gridded data, these are filtered data (mainly for reducing the size of the data objects) and should not be used as predictors in the actual analysis. For instance, air temperature reanalysis data are stored as a set of 20 empirical orthogonal functions (EOFs) with global coverage, which are then transformed to a field object upon a retrieval. However, the sample station data can be used without restrictions and corresponds to observed values at a specific location (e.g. `data(Oslo)` or `nacd=station(src='nacd')`).

11http://figshare.com/articles/esd\_for\_Mac\_amp\_Linux/1160493
12http://figshare.com/articles/esd\_for\_windows/1160494
13http://www.proclim.ch/imilast/index.html

## Function `retrieve()`
The ‘retrieve’ is designed to read data from NetCDF files following the standard Climate and Forecast (‘CF’) conventions and ordered on a longitude-latitude grid. The latter could be regular or irregular grid (e.g. the output of the Weather Research and Forecasting model (WRF) or any outputs on a rotated grid from RCMs). The function ‘retrieve’ also performs quick checks of the data itself to verify that the meta data and data are consistent. For instance, in case the frequency of the data is missing, `retrieve()` will try to detect the frequency automatically form the data itself. The function `retrieve` returns two types of objects depending on the type of the spatial grid. For instance, data stored on a regular grid is returned as `field` objects including attributes containing meta data, and data stored on irregular grid - such as rotated longitude-latitude grids - are returned as a `station` objects. 

The function `retrieve()` has also been adapted to read global climate data from the CMIP3/5 experiments, most of the global reanalysis such as those provided by the European Centre for Medium-Range Weather Forecasts (ECMWF) known as ERA-40 and ERA-INTERIM, the National Center for Environmental Prediction (NOAA) known as NCEP/NCAR reanalysis, the Japan Meteorological Agency (Japanese Reanalysis JRA-25,55), and the NASA GSFC Global Modeling and Assimilation Office (GMAO) known as MERRA reanalysis. A full overview of all available reanalysis can be found at http://reanalysis.org/atmosphere/overview-current-reanalyses.
As for now, `retrieve()`, does not display a list of missing attributes that are mandatory for further post-processing of the data. The user must add the missing attributes manually. The strength of `retrieve()` is that it can read and return formatted objects with common attributes for post-processing, visualising and making outputs comparable (e.g. re-gridding the field objects into the same grid resolution). Basically, all reanalysis, general circulation models (GCMs), and regional climate models (RCMs) can be read using the ‘esd’ tool and further combined into one object, or analysed, albeit with some limitations due to the size of
the returned object.

## Function `station()`
The package ‘esd’ includes the function ‘station’ for obtaining historical climate data by querying various web portals, for instance, the MET Norway archive (KDVH) provided by the Norwegian Meteorological Institute14. Data from MET climate web service ‘eKlima’ needs to be adapted manually, the Global Historical Climate Network (GHCN, (Peterson and Vose,1997)) provided by the American National Climatic Data Center15, and the European Climate
Assessment and Dataset (ECA&D, (Klein Tank et al., 2002)) made available by the Royal Netherlands Meteorological Institute (KNMI) (http://eca.knmi.nl/). Some of the data that is 14http://eklima.met.no; however, this function only works within the firewall included in the package has been pre-formatted within the ‘clim.pact’ package and adapted to meet the new ‘esd’ data format requirements.

15http://www1.ncdc.noaa.gov/pub/

SOURCE(S) :  ECAD/GHCNM/METNOD/METNOM/NACD/NORDKLIM
T2M/
14632
2014 / 1773
ESD package − map.station() − MET Norway 2014 (www.met.no)
Figure 1: Map of available weather stations recording temperature that are included in the meta-data
of the ‘esd’ package.

## Quick search - Function `select.station()`
The sample data includes also a meta-data object (`stationmeta`) that can be loaded directly and contains meta data such as name of the location (`loc`) and its standard ientification number (`stid`), geographical coordinates such as longitude (`lon`), latitude (`lat`), and altitude (`alt`), country name (`country`), parameter name (`param`) of recorded weather variables (e.g. temperature and precipitation), data source (`src`) or provider for thousands of stations all over the world (Figures 1 and 2).
These meta-data have been imported from existing meta data from the different data source providers. It has to be noted also that most of the available stations are managed by the World Meteorological Organisation (WMO) and made available by the American National Climate Data Centre (NCDC) for scientific research only. The meta data has been merged from the different sources mentioned earlier. Also, other additional data sources can be easily included into the package.

Figure 2: Map of available weather stations recording precipitation that are included in the meta-data of the ‘esd’ package.

There are two ways of obtaining station data using `station()` method. The first option, which we recommend, is to select a subset of stations from the meta data using the function `select.station()` from a set of criteria. Among these criteria, the minimum length of the recorded data (‘nmin’) can be specified to get, for instance, a subset of stations recording for at least a minimum number of (e.g. 100) years of data (e.g. `select.station(nmin=100)`). Thus, a subsample of the meta data is returned and can be checked for duplications of stations from the various sources. Afterwards, users can start downloading the data for the selected stations. It has to be mentioned that the downloading process for the GHCN and ECA&D is different as the data is stored differently. For the GHCN data sets, each station is stored separately and a direct retrieval can be made. We highly recommend that users perform a prior selection of stations for the GHCN data sets before starting the downloading process. For the ECA&D data sets on the other hand, all stations are downloaded and stored locally as zip files at first call, after which data is extracted for the selection of stations. Other data sets, such as the NACD (Frich et al., 1996), NARP (Førland ), and Nordklim (Tuomenvirta et al., 2001) are stored entirely in ‘esd’ (only holds monthly values of a limited set).

The ‘station’ method can retrieve station data from a range of different sources, many over the web (GHCN and ECA&D). Example 2.1 demonstrates how to select, retrieve and plot temperature observations from a single station. The ‘esd’ tool holds meta-data for various data sources, which makes it quick and easy to
search for weather data recorded at stations according to the parameter, the geographical location (region, country, and coordinates (single or a range of longitude and latitude values) and altitude, and time interval.
