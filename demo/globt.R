## Rasmus Benestad
## Script to plot the global mean temperature from GCMs, reanalysis, and
## observations on top of picture of the earth

## THIS SCRIPT ASSUMES ACCESS TO THE INTERNET

## This analysis uses a free analytical package written for R. The script
## will check whether it is installed; if not, it will install it
## for you. Its web address is http://github.com/metno/esd

## Check if you need to get the esd-package:
install.esd <- ("esd" %in% rownames(installed.packages()) == FALSE)

if (install.esd) {
  print('Need to install the esd package')
  ## Need online access.
  ## Use the devtools-package for simple facilitation of installing.
  if ("devtools" %in% rownames(installed.packages()) == FALSE)
    install.packages('devtools')
  library(devtools)
  ## Install esd directly from github
  install_github('metno/esd')
  print('The latest version of esd has been installed')
}

 if ("jpeg" %in% rownames(installed.packages()) == FALSE)
    install.packages('jpeg')

library(esd)
library(jpeg)

## Download image and data from the internet

url <- 'http://vignette2.wikia.nocookie.net/septimusheap/images/6/68/638831main_globe_east_2048.jpg/revision/latest?cb=20150827184925'
img <- "638831main_globe_east_2048.jpg"


if (!file.exists(img))
    download.file(url,img)

ftp <- 'ftp://ftp.cdc.noaa.gov/Datasets/ncep.reanalysis.derived/surface/air.mon.mean.nc'

if (!file.exists('air.mon.mean.nc'))
    download.file(ftp,'air.mon.mean.nc')

balls <- function(x) {
  for (i in 1:20) points(x,cex=seq(2,0.1,length=20)[i],
  col=rgb(i/20,i/20,i/20))
}

img <- readJPEG(img)

data(global.t2m.cmip5)

hadcrut <- annual(HadCRUT4())
ncep <- anomaly(aggregate.area(retrieve('air.mon.mean.nc'),FUN='mean'))
## anomalies:
ncep <- ncep - mean(subset(ncep,it=c(1961,1990)))

## Produce the graphics:
dev.new()
plot(hadcrut)
par(mar=rep(0,4))
plot(c(0,1),c(0,1),type='n')
rasterImage(img, -0.05, -0.05, 1.05, 1.05)

par(new=TRUE,col.axis='white',col.lab='white',xaxt='n',yaxt='n')
plot(hadcrut,lwd=5,col='black',ylim=c(-0.8,1.2),xlim=c(1950,2020),
     ylab=expression(T[2*m]*~(degree*C)))
for (i in 1:dim(global.t2m.cmip5)[2]) {
  offs <- mean(subset(hadcrut,it=c(1961,1990))) -
          mean(window(global.t2m.cmip5[,i],start=1961,end=1990))
  lines(global.t2m.cmip5[,i]+offs,lwd=7,col=rgb(1,0.7,0.7,0.1))
}
lines(hadcrut,lwd=5,type='b',pch=19)
lines(annual(ncep),lty=2,lwd=3)
balls(annual(ncep))
balls(hadcrut)

## Repeat the last anomaly
ncep.2015 <- subset(ncep,it=max(year(ncep)))
nm <- length(ncep.2015)
fcst <- zoo(rep(coredata(ncep.2015)[nm],12-nm),
                order.by=as.Date(paste(year(ncep.2015)[1],
		                       (nm+1):12,'01',sep='-')))
ncep.2015 <- c(zoo(ncep.2015),fcst)
points(annual(ncep.2015),col=rgb(0.9,0.9,0.9,0.3),cex=3,pch=19)
points(annual(ncep.2015),lwd=2,cex=3)
points(annual(ncep.2015),pch='?',cex=1.25,col='grey50',font=2)

par(xaxt='s',yaxt='s')
axis(1,col='white')
axis(2,col='white')

dev.copy2pdf(file='globalwarming.pdf')
