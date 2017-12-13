## ---- eval = FALSE-------------------------------------------------------
## ## Extract just the R-code
## library(knitr)
## purl('~/gmst-outlook.Rmd', output='~/git/esd/demo/gmst-outlook.R')

## ------------------------------------------------------------------------
library(esd)
download.file('ftp://ftp.cdc.noaa.gov/Datasets/ncep.reanalysis.derived/surface/air.mon.mean.nc',
              destfile = 'air.mon.mean.nc')

## ------------------------------------------------------------------------
t2m <- retrieve('~/Downloads/air.mon.mean.nc')
gmst <- aggregate.area(t2m,FUN='mean')
gmsta <- anomaly(gmst,ref=1961:1990)

## ------------------------------------------------------------------------
t <- index(gmst)
mon <- month(t[length(t)])
Y <- annual(gmsta)
x <- aggregate(zoo(subset(gmsta,it=month.abb[1:mon])),year,FUN='mean')
z <- lm(coredata(Y) ~ coredata(x))
y <- predict(z,newdata=x)
n <- length(y)

## ------------------------------------------------------------------------
ylim<- range(c(coredata(Y),coredata(y)),na.rm=TRUE) + 0.25*c(-1,1)
plot(zoo(x),zoo(Y),ylim=ylim,
     xlab=paste('Jan-',month.abb[mon],' mean',sep=''),
     ylab='annual mean', sub='base line: 1961-1990; source NCEP/NCAR reanalysis 1',
     main='Global mean temperature',pch=19,cex=1.5)
lines(rep(coredata(x)[n],2),ylim,lwd=7,col='grey90')
arrows(coredata(x)[n],coredata(y)[n],min(x,na.rm=TRUE),coredata(y)[n],lty=2)
points(coredata(x)[n],coredata(y)[n],cex=1.5,lwd=2)
text(coredata(x)[n],min(y),max(year(x)),sub=1)
grid()

## ------------------------------------------------------------------------
map(subset(anomaly(t2m,ref=1961:1990),it=length(index(t2m))),colbar=list(breaks=seq(-20,20,by=1)))

