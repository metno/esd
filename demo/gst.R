## R-script to plot global mean temperature
library(esd)

## Download the latest reanalysis: near-surface temperature
download.file('ftp://ftp.cdc.noaa.gov/Datasets/ncep.reanalysis.derived/surface/air.mon.mean.nc','air.mon.mean.nc')
t2m <- retrieve('air.mon.mean.nc')

## Estimate annual mean of global mean
gst <- zoo(annual(aggregate.area(t2m,FUN='mean')))
gst <- gst - mean(window(gst,start=1961,end=1990))
data("global.t2m.cmip5")
data("global.t2m.cmip3")
par(bty='n')
plot(gst,lwd=3,type='b',main='Global mean temperature',
     sub='NCEP/NCAR reanalysis 1: baselin=1961-1990',
     ylab=expression(T*phantom(0)*(degree*C)),
     xlab='year',ylim=c(-0.5,1))
for (i in 1:dim(global.t2m.cmip5)[2]) 
  lines(global.t2m.cmip5[,i] - mean(window(global.t2m.cmip5[,i],start=1961,end=1990)),
        col=rgb(0.5,0.2,0.1,0.1))
for (i in 1:dim(global.t2m.cmip3)[2]) 
  lines(global.t2m.cmip3[,i] - mean(window(global.t2m.cmip3[,i],start=1961,end=1990)),
        col=rgb(0.1,0.2,0.5,0.1))
lines(gst,lwd=3,type='b')
grid()
legend(1950,1,c('reanalysis','CMIP5 RCP4.5','CMIP3 SRESA1b'),lty=1,
       lwd=c(3,1,1),col=c('black',rgb(0.5,0.2,0.1),rgb(0.1,0.2,0.5)),
       bty='n')
lines(c(1961,1990),rep(0,2),lwd=3,col=rgb(0.5,0.5,0.1,0.25))
