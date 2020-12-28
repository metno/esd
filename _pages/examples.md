---
layout: default
title: Examples
nav_order: 9
permalink: /examples/
---

# Examples
{: .no_toc }

<details open markdown="block">
  <summary>
    Table of contents
  </summary>
  {: .text-delta }
1. TOC
{:toc}
</details>


## Example 2.1 on how to select weather stations and retrieve recorded values
```R
# Select a station across India recording daily maximum temperature
# from the global historical climate network-daily
ss <- select.station(cntr=’India’,param=’tmax’,src=’ghcnd’)
ss.new <- subset(ss,subset=!duplicated(ss$location))
map(ss.new,cex=1.2,col="red",bg="pink",add.text=TRUE)
ss <- select.station(stid=’IN022021900’,cntr=’india’,param=’tmax’)
y <- station(ss)
#>[1] "Retrieving data ..."
#>[1] "1 TMAX IN022021900 NEW DELHI/S INDIA GHCND"
# Display the name of the location
loc(y)
# Highlights the location on the map
points(lon(y),lat(y),pch=19,col="blue",cex=2)
# aggregate daily values to annual values
ya <- annual(y,FUN=’mean’,nmin=100)
# Subset for the period 1970 to 2012
ya <- subset(ya,it=c(1970,2012))
# plot the time series including the error bar
plot(ya,ylim=c(29.5,32.5))
# Add the linear trend as
lines(trend(ya),col="red",lwd=2)
```

Figure 3: Map showing available stations across India from the GHCN-D dataset and b) a plot of
the annual maximum temperature recorded at New Delhi weather station (blue point) including
linear trend line.

## Example 2.2 on the different esd objects and clasees
```R
# Example of monthly station data:
> data(Oslo)
> class(Oslo)
[1] "station" "month"
"zoo"
# Example of daily station data:
> data(ferder)
> class(ferder)
[1] "station" "day"
"zoo"
# Example of annual station data:
> class(annual(ferder))
[1] "station" "annual"
"zoo"
# Example of a field object
> t2m <- t2m.NCEP(lon=c(-30,30),lat=c(40,70))
> class(t2m)
[1] "field" "month" "zoo"
> class(EOF(t2m))
[1] "eof"
"field" "month" "zoo"
Example 2.3. # Load the data for Ferder weather station
data(ferder)
# Display the structure of the data as
> str(ferder)
’zoo’ series from 1900-01-01 to 2013-12-04
Data: atomic [1:41611] 2.1 -1.8 -0.9 -3 -7.2 -6.5 -2.6 -2.4 -1.6 -0.3 ...
- attr(*, "location")= chr "Ferder lighthouse"
- attr(*, "station_id")= num 27500
- attr(*, "wmo_id")= num 1482
- attr(*, "longitude")= num 10.5
- attr(*, "latitude")= num 59
- attr(*, "altitude")= num 6
- attr(*, "variable")= chr "t2m"
- attr(*, "unit")= chr "deg C"
- attr(*, "country")= chr "Norway"
- attr(*, "source")= chr "MET Norway eklima"
- attr(*, "long name")= chr "daily mean temperature"
- attr(*, "URL")= chr "http://eKlima.met.no"
- attr(*, "calendar")= chr "gregorian"
- attr(*, "quality")= logi NA
- attr(*, "aspect")= chr "original"
- attr(*, "type")= chr "observation"
- attr(*, "reference")= chr "MET Norway climate archive"
- attr(*, "info")= logi NA
- attr(*, "method")= chr "mean estimated from thermometer measurements"
- attr(*, "history")=List of 3
..$ call
: chr "1 0"
..$ timestamp
: chr "Fri Dec 13 15:51:49 2013"
..$ sessioninfo:List of 3
.. ..$ R.version
: chr "R version 3.0.2 (2013-09-25)"
.. ..$ esd.version: chr "esd_0.2-1"
.. ..$ platform
: chr "x86_64-pc-linux-gnu (64-bit)"
Index:
Date[1:41611], format: "1900-01-01" "1900-01-02" "1900-01-03" "1900-01-04" ...
```

## Example 2.4 on how to load some data sets and display the summary of statistics
```R
## Load data for Ferder
> data(ferder)
## Display the summary of statistics
> summary(ferder)
```
## Example 2.5 on how to load a weather station object and display the location on a map
```R
# Load bjornholt data set
data(bjornholt)
# plot the time series
plot(bjornholt)
# map the location
map(bjornholt)
```

Figure 4: Example of a station (a) and a map plot (b) for a single station ‘Bjornholt’.

## Example 2.6 on how to read a field object and plot the results
```R
#Get NCEP 2m air temperature for the selected spatial window defined by lon and lat
t2m <- t2m.NCEP(lon=c(-30,30),lat=c(40,70))
# Computes the EOFs
X <- EOF(t2m)
# Plot the result
plot(X)
```

Figure 5: An example of a plot results for an EOF object.

## Example 2.7 on how to retrieve wind data sets and visualise the results
```R
# Load 10m zonal and meridional wind components
u10 <- retrieve(’data/ERAINT/eraint_elnino.nc’,param=’u10’)
v10 <- retrieve(’data/ERAINT/eraint_elnino.nc’,param=’v10’)
# Map the data
map(u10,colorbar=FALSE)
# Display the vectors
vec(u10,v10,new=FALSE,a=2,length=0.05)
```

## Example 2.8. test
```R
# Get 2m temperature data for Oslo
# (works within MET Norway firewall only)
x <- station(stid=18700,param=’t2m’,src=’metnod’)
x <- subset(x,it=c(1837,2014))
## Or gets data from ECA&D as
## x <- station(loc="Oslo Blindern",param=’t2m’,src=’ecad’)
## x <- subset(x,is=duplicated(loc(x))) # to remove duplicated stations
# Cumulative average
cumugram(x)
# Seasonal wheel
wheel(x)
# seasonal variations of year-to-year variance
climvar(x)
# daily seasonal cycle for all years
diagram(x)
```

Figure 7: Examples of ‘cumugram’, ‘wheel’, ‘climvar’, and ‘diagram’ plots. The ‘cumugram’ shows how the mean value of some variable has evolved from the start of the year and compares this curve to previous years. The graphics produced by ‘wheel’, on the other hand, emphasises how the seasonal variations affect the variable, e.g. whether some extremes tend to be associated with a specific season. Panel c shows results produced by ‘climvar’ shows the year-to-year statistics for a variable, e.g. the standard deviation of the temperature on February 1st. The ‘diagram’ method can be used in different context, and for a ‘station’ object, it produces graphics that compare the day-to-day values with those of previous years.

## Example 2.9 on how to process the data and visualise the trend of the time series
```R
# Get 2m temperature data for Ferder and calculate annual mean
data(ferder)
x <- annual(ferder)
# Plot the time series and add a trend line
plot(x,ylim=c(4,11))
lines(trend(x),col="red",lwd=2)
# Visualise trends for various periods 40 years or longer (minlen)
# and mark trends that are significant at the 1% level (pmax)
vis.trends(x,minlen=40,pmax=0.01)
t2m trend (deg C/decade)
```
Figure 8: Visualising the trend in the temperature time series from Ferder for the full period (left panel) and simultaneously for different periods of various length (right panel). The trend analysis indicates a significant increase in temperature at Ferder since the 1900s, but also shows that the strength and significance of the trend is sensitive to the period considered because the increase is not monotone.
