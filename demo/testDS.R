## Test script for DS on one station data 

## ------------------- Example 1 --------------------------
library(esd)
lon <- c(-12,37)
lat <- c(52,72)

## If not stored locally, fetch the predictor data from FTP
if (!file.exists('air.mon.mean.nc')) {
  url <- 'ftp://ftp.cdc.noaa.gov/Datasets/ncep.reanalysis.derived/surface/air.mon.mean.nc'
  download.file(url,'air.mon.mean.nc')
}

## Retrieve the NCEP/NCAR reanalysis
t2m <- retrieve('air.mon.mean.nc',lon=lon,lat=lat,type='ncdf4')

## Get the predictand
data(Oslo)
## Prepare the data for downscling of winter temperature
y <- subset(as.4seasons(Oslo),it='djf')
X <- EOF(subset(as.4seasons(t2m),it='djf'))
z <- DS(y,X)
plot(z)

## ------------------- Example 2 --------------------------

## Downscaling using common EOF
ylim <- c(-6,6)
t2m <- t2m.NCEP(lon=lon,lat=lat)
T2m <- t2m.NorESM.M(lon=lon,lat=lat)
X <- combine(t2m,T2m)

eof <- EOF(X,it=7) # common EOF = EOF of matrix containing combined data (anomalies)
ds <- DS(Oslo,eof)
plot(ds)

## ------------------- Example 3 --------------------------

## Downscale seasonal mean and standard deviation
data(ferder)
t2m <- t2m.NCEP(lon=c(-30,50),lat=c(40,70))
slp <- slp.NCEP(lon=c(-30,50),lat=c(40,70))
T2m <- as.4seasons(t2m)
SLP <- as.4seasons(slp)
X <- EOF(T2m,it='jan')
Z <- EOF(SLP,it='jan')
y <- anomaly(ferder)
ym <- as.4seasons(y,FUN="mean")
ys <- as.4seasons(y,FUN="sd")
dsm <- DS(ym,X)
dss <- DS(ys,Z)

