---
title: "Examples"
layout: home
toc: true
classes: wide
permalink: /examples/
---

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

## Example 2.3. # Load the data for Ferder weather station
```R
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
..$ timestamp : chr "Fri Dec 13 15:51:49 2013"
..$ sessioninfo:List of 3
.. ..$ R.version : chr "R version 3.0.2 (2013-09-25)"
.. ..$ esd.version: chr "esd_0.2-1"
.. ..$ platform
: chr "x86_64-pc-linux-gnu (64-bit)"
Index: Date[1:41611], format: "1900-01-01" "1900-01-02" "1900-01-03" "1900-01-04" ...
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

## Example 2.6 on how to read a field object and plot the results
```R
#Get NCEP 2m air temperature for the selected spatial window defined by lon and lat
t2m <- t2m.NCEP(lon=c(-30,30),lat=c(40,70))
# Computes the EOFs
X <- EOF(t2m)
# Plot the result
plot(X)
```

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

## Example 3.1. # Load data for Bjornholt 
```R
data(bjornholt)
y <- bjornholt
# Annual wet-mean:
mu <- annual(y,FUN=’wetmean’)
# Plot the result for the period 1883 to 2014
plot(mu,ylim=c(4,14))
# Annual wet-freq:
fw <- annual(y,FUN=’wetfreq’)
# Plot the result for the period 1883 to 2014
plot(subset(fw,it=c(1883,2014)),ylim=c(0.2,0.5))
# Annual number of events with y > 10mm/day:
nw <- annual(y,FUN=’count’,threshold=10)
# Plot the result for the period 1883 to 2014
plot(subset(nw,it=c(1883,2014)),ylim=c(10,60))
# Compute wet/dry spells
sp <- spell(y,threshold=1)
# Plot the result
plot(sp)
```

## Example 3.2.
```R
# Load ECAD dataset over Norway recording a 100 years of daily data
ecad <- station(param="precip",src="ecad",cntr="Norway",nmin=100)
# Map the vulnerability index:
map(ecad,FUN=’precip.vul’,cex=1.2,xlim=c(-10,40), col="black",colbar=list(col=heat.colors(20)))
# Map approximate 10-year return values:
map(ecad,FUN=’precip.rv’,cex=1.2,xlim=c(-10,40), col="black",colbar=list(col=heat.colors(20)))
```

## Example 3.3.
```R
# Load NCEP 2m air temperature
t2m <- t2m.NCEP(lon=c(-30,30),lat=c(40,70))
# map the original field
map(t2m)
# Regrid based on the new lon and lat values
y <- regrid(t2m,is=list(lon=seq(-5,15,by=0.5),lat=seq(55,65,by=0.5)))
# Map on the new grid
map(y)
```

## Example 3.4.
```R
# Get NCEP data
> ncep <- t2m.NCEP(lon=c(-15,45),lat=c(35,70))
# Display the latitude values
> lat(ncep)
# Compute the grid resolution
> diff(lat(ncep))
# Display the longitude values
> lon(ncep)
# Compute the resolution in the lon axis
> diff(lon(ncep))
# Get GCM data
> gcm <- t2m.NorESM.M(lon=c(-15,45),lat=c(35,70))
> lat(gcm)
> diff(lat(gcm))
> lon(gcm)
> diff(lon(gcm))
# Do the re-gridding # this line does not work
> gcm.regrid <- regrid(gcm, is=ncep)
# Make sure that both longitudes are set to data line (i.e. greenwich=FALSE)
> ncep <- g2dl(ncep,greenwich=FALSE)
> gcm <- g2dl(gcm,greenwich=FALSE)
# Do the re-gridding
> gcm.regrid <- regrid(gcm, is=ncep)
> lat(gcm.regrid)
> lon(gcm.regrid)
```

## Example 3.5.
```R
data(Oslo)
# Extract an interval:
y <- subset(Oslo,it=as.Date(c("1883-01-01","2013-12-05")))
# Extract only the winter data (use aggregate for winter statistics)
djf <- subset(y,it=’djf’)
# Extract data for May:
may <- subset(y,it=’May’)
```

## Example 3.6.
```R
# Retrieve stations across Scandinavian regions from the
# ECA$\&$D dataset with a minimum of 50 years of data
y <- station(src="ecad",cntr="norway",nmin=50)
# Show the selected stations including all available stations
map(y,cex=1.4,col="red",bg="pink",showall=TRUE)
# Subset stations with altitude higher than 100m
y1 <- subset(y,is=list(alt=100))
# Show the stations with elevation greater than 100m above sea level:
map(y1,cex=1.4,col="darkred",bg="red",showall=TRUE)
# Show the stations with elevation below 100m above sea level:
y2 <- subset(y,is=list(alt=-100))
map(y2,cex=1.4,col="darkred",bg="orange",showall=TRUE)
```

## Example 3.7.
```R
# Load data for "Bjornholt" station
data(bjornholt)
# Check the class of the object
class(bjornholt)
## [1] "station" "day"
"zoo"
# Aggregate on annual maximum values
bjornholt <- annual(bjornholt,FUN=’max’)
# Check the class of the aggregated object
class(bjornholt)
## [1] "station" "annual"
"zoo"
# Plot the results
plot(bjornholt)
```

## Example 3.8.
```R
# Load 2m air temperature from NCEP reanalysis on a resoltuion of 2.5 deg.
t2m <- t2m.NCEP()
# Do the spatial aggregation
x <- aggregate(t2m,by=list(lon=seq(0,360,by=10),lat=seq(0,360,by=10)),FUN=’mean’)
# Map the results
map(x)
```

## Example 3.9.
```R
# Load 2m air temperature from the NorESM.M global climate model
t2m <- t2m.NorESM.M()
# Compute the areal mean over the whole domain
T2m <- aggregate.area(t2m,FUN=’mean’)
# Plot the annual aggregated values
plot(annual(T2m),ylim=c(11.5,17.5))
```

```R
ds <- DS(Oslo,ceof)
# Subset the calibration
z.calibration <- predict(ds,newdata=EOF(t2m))
# Extract projected data sets
z.projection <- project.ds(ds)
```

## Example 4.1.
```R
# Load 2m air temperature from NCEP reanalysis
t2m <- t2m.NCEP()
# Compute the EOFs
eof <- EOF(t2m)
# Plot the eof
plot(eof)
```

## Example 4.2.
```R
# Retrieve NACD temperature weather stations
nacd <- station(src=’nacd’,param=’t2m’)
# Compute the annual mean values
NACD <- annual(nacd)
# Compute the number of valid data points for each station
ok<- apply(NACD,2,FUN=’nv’)
# Retain only the stations with minimum of 100 valid data points
NACD <- subset(NACD,is=(1:length(ok))[ok > 100])
# Do the PCA
pca <- PCA(NACD)
# Visualise the results
vis(pca)
```

## Example 4.3.
```R
# Get NACD stations
nacd <- station(src=’nacd’,param=’t2m’)
# Aggregate to annual values
NACD <- annual(nacd)
# Check for missing values
ok<- apply(NACD,2,FUN=’nv’)
# Subset NACD stations with more than 100 data points
NACD <- subset(NACD,is=(1:length(ok))[ok > 100])
# Compute the PCAs
pca <- PCA(NACD)
# Retrieve slp field and aggregate to annual values
slp <- annual(slp.DNMI())
# compute EOF of slp
eof <- EOF(slp)
# compute CCA on both pca and eof
cca <- CCA(pca,eof)
plot(cca)
slp(hPa)
```

## Example 4.4. # Sample temperature data for Oslo
```R
data(Oslo)
# Get ERA40 2m air temperature
t2m <- t2m.ERA40(lon=c(0,20),lat=c(50,70))
# Get NorESM.M 2m air temperature
T2m <- t2m.NorESM.M(lon=c(lon=c(0,20),lat=c(50,70)))
# Combine the two fields
comb <- combine(t2m,T2m)
# Compute the common eof t2m
ceof <- EOF(comb)
# Do the downscaling of y based on X
ds <- DS(Oslo,ceof)
# Subset the calibration
z.calibration <- predict(ds,newdata=EOF(t2m))
# Extract projected data sets
z.projection <- project.ds(ds)
```

## Example 4.5.
```R
# Sample trajectory object
data(imilast.M03)
# Select storm trajectories in region 40-60N
x <- subset(imilast.M03,is=list(lat=c(40,60)))
# Select winter storms
djf <- subset(x,it=’djf’)
# Calculate the seasonal storm count
n <- count.trajectory(djf,by=’year’)
# Plot the annual storm count
plot(x,new=TRUE,ylim=c(0,60),col=’black’)
# Add storm count for winter (djf)
lines(n,col=’blue’,lty=2)
legend(’topright’,c(’all year’,’djf’),lty=c(1,2),col=(’black’,’blue’))
```

## Example 4.6. # Sample trajectory object
```R
data(imilast.M03)
#Map storm trajectories for the winter season (djf)
map(imilast.M03,it=’djf’,projection=’sphere’)
# Map number density of storms (per year and 1000km2)
map.density.trajectory(imilast.M03,it=’djf’, xlim=c(-90,60),ylim=c(30,90),dx=4,dy=2)
```

## Example 4.7.
```R
# Sample trajectory object
data(imilast.M03)
# Perform PCA of storm tracks
pca <- PCA.trajectory(imilast.M03,anomaly=TRUE)
# Plot PCA
plot.pca.trajectory(pca)
# Show PCA on map
map.pca.trajectory(pca)
```

## Example 5.1.
```R
# Load eofs of NCEP 2m air temperature
data(eof.t2m.NCEP)
# Load eofs of NCEP sea level pressure
data(eof.slp.NCEP)
# Load temperature data for Oslo
data(Oslo)
# Do the downscaling
z <- DS(Oslo,list(t2m=eof.t2m.NCEP,slp=eof.slp.NCEP),mon=1)
# Plot the results
plot(z)
```

## Example 5.2.
```R
# Load ERAINT 2m air temperature from file
T2m <- t2m.ERAINT(lon=c(-10,30),lat=c(50,70))
# Compute the EOFs on annual values
eof <- EOF(annual(T2m))
# Retrieve Nordklim data sets
y <- station(src=’nordklim’)
# Get rid of stations with little valid data
Y <- allgood(annual(y))
# Do the PCA
pca <- PCA(Y)
# DO the downscaling
z <- DS(pca,eof)
# Plot the result
plot(z)
# Get the observations corresponding to the first station in the list
obs <- subset(Y,is=1)
# Get the predictions
pre <- subset(as.station(z),is=1)
# Match dates between obs and pre
obs <- matchdate(obs,pre)
# Update the attributes
attr(obs,’location’) <- ’Observed’
attr(pre,’location’) <- ’Predicted’
# Combine the two
obspre <- combine.stations(obs,pre)
# Plot the result in a single plot
plot(obspre,plot.type="single")
```

## Example 5.3.
```R
# Get the predictand: Maximum temperature for some Norwegian sttions
download.file(url="http://files.figshare.com/2073466/Norway.Tx.rda",destfile="NorwayTx.rda")
load("NorwayTx.rda")
# Process the predictand (annual mean)
ma <- annual(Tx)
pca <- PCA(ma)
# Get the ERAINT T(2m):
t2m <- t2m.ERAINT(lon=c(-30,30),lat=c(50,75))
# Process the predictor: Compute EOFs of annual means
eof <- EOF(annual(X))
# Test: CCA:
cca <- CCA(pca,eof)
# Do the downscaling
z <- DS(pca,eof,rmtrend=FALSE)
# Plot the results for one station:
plot(subset(ma,is=1))
lines(subset(as.station(pca),is=is),lty=2,col="darkred")
lines(subset(as.station(z),is=is),lwd=2)
```

## Example 5.4.
```R
# Load storm trajectories
data(imilast.M03)
# Load NCEP sea level pressure data for the northern hemisphere
slp <- slp.NCEP(lat=c(0,90))
# Compute eof
eof.slp <- EOF(slp)
# Do downscaling for the storm count
(default trajectory downscaling)
ds <- DS(imilast.M03,eof.slp)
# Do downscaling from slp
ds.slp <- DS(imilast.M03,eof.slp,param=’slp’,FUN=’min’,unit=’hPa’)
# Plot downscaled results
plot(ds)
plot(ds.slp)
```

## Example 6.1.
```R
# Retrieve NACD data sets
nacd <- station(src=’nacd’)
# Diagnose for data availability
diagnose(nacd)
Example 6.2.
# Sample 2m air temperature ERAINT field from EOF’s data
t2m <- t2m.ERAINT(lon=c(-30,40),lat=c(50,70))
# Retrieve GCM air temperature from file
# Sample 2m air temperature NorESM field from EOF’s data
gcm <- t2m.NorESM.M()
# Combine annual fields
X <- combine(annual(t2m),annual(gcm))
# Diagnose the result
diagnose(X)
Example 6.3.
# Compute EOF, X is taken from previous example
eof <- EOF(X)
# Diagnose eof
a <- diagnose(eof)
# Plot the results
plot(a)
```

## Example 6.4. # Load temperature data for Oslo
```R
data(Oslo)
# Get the 2m air temperature from NCEP reanalysis
t2m <- t2m.NCEP(lon=c(-5,20),lat=c(50,65))
# Compute the EOF on annual temperature values
eof <- EOF(annual(t2m))
# Downscale the annual values
z <- DS(annual(Oslo),eof)
# Do the diagnostic
diagnose(z,plot=TRUE)
Example 6.5. # Download the file from the figshare link as
download.file(url="http://files.figshare.com/2081407/dse.oslo.rda",destfile="dseOslo.rda")
# Load the file into R session
load(’dseOslo.rda’)
# Do the diagnostic
diagnose(dse.oslo)
```
