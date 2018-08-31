library(esd)

## Get the data for the 2-meter temperature from NCEP/NCAR reanalysis 1:
if (!file.exists('air.mon.mean.nc')) download.file('ftp://ftp.cdc.noaa.gov/Datasets/ncep.reanalysis.derived/surface/air.mon.mean.nc','air.mon.mean.nc')
t2m <- retrieve('air.mon.mean.nc')
sigma <- 5.67e-8
## Estimate the longwave radiation from the surface using the Stefan-Boltzman's law
slr <- sigma*(t2m+273.15)^4
attr(slr,'variable') <- 'slr'
attr(slr,'unit') <- 'W/m^2'
map(slr,projection='sphere',colbar=colbar)


## Get the outgoing longwave radiation from NOAA Earth Systems Research Laboratory
if (!file.exists('olr.mon.mean.nc')) download.file('http://www.esrl.noaa.gov/psd/repository/entry/get/olr.mon.mean.nc?entryid=synth%3Ae570c8f9-ec09-4e89-93b4-babd5651e7a9%3AL2ludGVycF9PTFIvb2xyLm1vbi5tZWFuLm5j','olr.mon.mean.nc')
olr <- retrieve('olr.mon.mean.nc')
colbar <- list(pal='heat',breaks=seq(100,500,by=10))
map(olr,projection='sphere',colbar=colbar)
