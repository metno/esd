## Huricanes and tropical cyclones
## http://www.aoml.noaa.gov/hrd/hurdat/hurdat2-1851-2014-022315.html
## Review and update
## Benestad, R.E. 'On Tropical Cyclone Frequency and the Warm Pool Area'
## Nat. Hazards Earth Syst. Sci., 9, 635-645, 2009
## http://www.nat-hazards-earth-syst-sci.net/9/635/2009/nhess-9-635-2009.html
## N_{Atl} \propto A^{5.06 \pm 0.25}
rm(list=ls())
library(esd)

update <- TRUE
process <- TRUE
est.area <- TRUE
temperature.weighted <- FALSE
exceedance <- TRUE
T.crit <- 26.5
#last.year <-2010
natl.lon <- c(-80,10); natl.lat <- c(0,40)
car.lon <- c(-100,-80); car.lat <- c(15,30)
fname <- 'sst.mnmean.v4.nc'
fname2 <- 'sst.mon.mean.nc'

## Function that returns standardised data for comparison:
stand <- function(x,m,s) (x - m)/s
## Only estimate maximum warm area for seasons with more than 6 months
maxtcs <- function(x) {if (length(x)==6) y <- max(x) else y <- NA; return(y)}

#source("strip.R")

print("The Tropical North Atlantic:")
## Get the latest sea surface temperatures (SST) from NOAA:
if (!file.exists(fname)) {
  download.file('ftp://ftp.cdc.noaa.gov/Datasets/noaa.ersst/sst.mnmean.v4.nc',fname)
}
if (!file.exists(fname2)) {
  download.file('https://psl.noaa.gov/repository/entry/get/sst.mon.mean.nc?entryid=synth%3Ae570c8f9-ec09-4e89-93b4-babd5651e7a9%3AL0NPQkUvc3N0Lm1vbi5tZWFuLm5j',fname2)
}

## NINO3.4 index is defined as the average of SST anomalies over the region 5°N - 5°S and 170° - 120°W. 
## NCC classifies the NINO3. 4 temperature anomaly as "warm" if it exceeds 0.8°C, which is about one 
## standard deviation above average.
#nino3.4 <- aggregate(subset(NINO3.4(),it=c('Sep','Oct','Nov')),year,'mean')
nino3.4 <- retrieve(fname2,lon=c(-170,-120),lat=c(-5,5))
nino3.4 <- aggregate.area(annual(anomaly(nino3.4,ref=1981:2010)),FUN='mean')
index(nino3.4) <- year(nino3.4)

sst <- retrieve(fname,lon=natl.lon,lat=natl.lat)
sst2 <- retrieve(fname2,lon=natl.lon,lat=natl.lat)

## Area of boxes used to estimate area of warm surface
boxarea <- rep(cos(pi*lat(sst)/180)* pi*min(diff(lon(sst)))/180*
               pi*min(diff(lat(sst)))/180*(6.378e03)^2,length(lon(sst)))
dim(boxarea) <- c(length(lat(sst)),length(lon(sst)));
boxarea2 <- rep(cos(pi*lat(sst2)/180)* pi*min(diff(lon(sst2)))/180*
                 pi*min(diff(lat(sst2)))/180*(6.378e03)^2,length(lon(sst2)))
dim(boxarea2) <- c(length(lat(sst2)),length(lon(sst2)));
# image(lon(sst),lat(sst),t(boxarea))
# contour(lon(sst),lat(sst),t(boxarea),add=TRUE)

## Estimate the area of the warm ocean: first by masking land areas
sst <- mask(sst,land=TRUE)
sst <- aggregate(subset(sst,it=month.abb[6:11]),year,'maxtcs')
sstw <- coredata(sst)
sstw[sstw <T.crit] <- NA
sstw[sstw >=T.crit] <- 1
sstw <- t(t(sstw)*c(t(boxarea)))
coredata(sst) <- sstw
warmarea <- zoo(apply(sst,1,sum,na.rm=TRUE),order.by=index(sst))
index(warmarea) <- year(warmarea)
## Repeat for the sst.mon.mean.nc data:
sst2 <- mask(sst2,land=TRUE)
sst2 <- aggregate(subset(sst2,it=month.abb[6:11]),year,'maxtcs')
sstw2 <- coredata(sst2)
sstw2[sstw2 <T.crit] <- NA
sstw2[sstw2 >=T.crit] <- 1
sstw2 <- t(t(sstw2)*c(t(boxarea2)))
coredata(sst2) <- sstw2
warmarea2 <- zoo(apply(sst2,1,sum,na.rm=TRUE),order.by=index(sst2))
## Fudge - vmissing values rather than zero
warmarea[warmarea == 0] <- NA
warmarea2[warmarea2 == 0] <- NA
index(warmarea2) <- year(warmarea2)

# Data on tropical cyclones:
## Read the data from the Internet
url <- 'http://www.aoml.noaa.gov/hrd/hurdat/hurdat2-1851-2014-022315.html'
hurdat2 <- readLines(url)
writeLines(hurdat2,con='hurdat2-1851-2014-022315.html')

## Extract the storms
itc <- grep('^AL',hurdat2,perl=TRUE)
ihu <- grep(', HU,',hurdat2,perl=TRUE)
its <- grep(', TS,',hurdat2,perl=TRUE)
## Ignore events which did not become hurricanes or tropical storms, last event is neither
lhu <- rep(FALSE,length(itc)); lts <- lhu
for (ii in 1:(length(itc)-1)) {
  lhu[ii] <- (sum( (itc[ii] < ihu) & (itc[ii+1] > ihu)) >0)
  lts[ii] <- (sum( (itc[ii] < its) & (itc[ii+1] > its)) >0)
}
itc <- itc[lhu | lts]
storms <- hurdat2[itc]
yr <- as.integer(substr(storms,5,8))
stormstats <- table(yr)
ntc <- zoo(as.numeric(stormstats),order.by=as.numeric(rownames(stormstats)))

## Use the relationship from Benestad (2009):
y <- warmarea^5.06
y2 <- warmarea2^5.06
scl <- mean(window(ntc,start=1961,end=1990))/mean(window(y,start=1961,end=1990))
scl2 <- mean(window(ntc,start=1961,end=1990))/mean(window(y2,start=1961,end=1990))
## Use a smoothing filter
y <- stats::filter(y*scl,rep(1,7)/7)
y2 <- y2 * scl2
yu <- stats::filter(warmarea^5.31*scl,rep(1,7)/7)
yl <- stats::filter(warmarea^4-81*scl,rep(1,7)/7)

caldat <- data.frame(ntc=window(ntc,start=1900,end=1960),
                      nino3.4=window(nino3.4,start=1900,end=1960),
                      warmarea=window(y,start=1900,end=1960))
i1 <- is.element(year(nino3.4),year(y2))
i2 <- is.element(year(y2),year(nino3.4))
predat <- data.frame(nino3.4=nino3.4[i1],warmarea=y2[i2])
model <- glm(ntc ~ warmarea + nino3.4 + I(nino3.4^2) +
             I(nino3.4^3) + I(nino3.4^4),
             data=caldat,family='poisson')
print(summary(model))
cntc <- zoo(exp(predict(model)),order.by=1900:1960)
xntc <- zoo(exp(predict(model,newdata=predat)),order.by=index(nino3.4))

par(bty='n',xpd=TRUE,las=3)

## Update the record of number of named cyclones manually: 2015-2023
## From Wikipedia
ntc2 <- zoo(c(11,15,17,15,17,30,21,14,20),order.by=2015:2023)
ntc <- c(ntc,ntc2)

plot(ntc,lty=2,xlim=range(index(y)),
     ylab=expression(paste(N[TC],' per season')),xlab="",
     pch=15,cex=1.2,type='b',
     main="Number of tropical cyclones in the North Atlantic",
     sub='Tropical Storms and Hurricanes')
grid()
rect(1900,0,1960,25,col=rgb(0.7,0.7,0.7,0.2),border=NA)
text(1905,24,'calibration',pos=4,col=rgb(0.3,0.3,0.7),font=2)
text(1880,25.5,'Nat. Hazards Earth Syst. Sci., 9, 635-645, 2009',
     side=4,cex=0.7,col='green')
lines(cntc,lwd=3,col=rgb(0.5,0.5,0.5,0.5))
lines(xntc,lwd=3,col=rgb(0.5,0.5,1,0.5))
lines(y,lwd=7,col=rgb(0,1,0,0.5))
lines(y2,lwd=2,col=rgb(0,1,0,0.5))
polygon(c(index(yu),rev(index(yl))),
       c(coredata(yu),rev(coredata(yl))),
       col=rgb(0,1,0,0.7))
mtext(url,side=4,cex=0.7)
legend(1960,3.5,c('HURDAT',expression(n%prop% A^5.06),
          expression(ln(hat(n))==beta[0]+beta[1]*A^5.06+
              beta[2]*phantom(0)*f(NINO3.4))),lty=c(2,1,1),
       lwd=c(1,7,5),pch=c(15,26,26),
       col=c('black','green',rgb(0.5,0.5,1,0.7)),bty='n',cex=0.7)       
legend(1850,24,c('observed','Estimated','Estimated'),lty=c(2,1,1),
       lwd=c(1,7,5),pch=c(15,26,26),
       col=c('black','green',rgb(0.5,0.5,1,0.7)),bty='n',cex=0.7)


