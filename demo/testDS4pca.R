##  Test DS for PCA data.

library(esd)
it <- 'djf'

## Test data - Norwegian temperature series.
## Only download if the data is not stored locally
if (!file.exists('Norway.Tx.rda')) {
  url <- 'http://files.figshare.com/2073466/Norway.Tx.rda'
  download.file(url,'Norway.Tx.rda')
}

## If not stored locally, fetch the predictor data from FTP
if (!file.exists('air.mon.mean.nc')) {
  url <- 'ftp://ftp.cdc.noaa.gov/Datasets/ncep.reanalysis.derived/surface/air.mon.mean.nc'
  download.file(url,'air.mon.mean.nc')
}

load('Norway.Tx.rda')
## Problem with the unit - quick fix.
attr(Tx,'unit') <- rep('deg C',length(unit(Tx)))

## Retrieve the NCEP/NCAR reanalysis
t2m <- retrieve('air.mon.mean.nc',lon=c(-20,20),lat=c(55,70))

## Get the winter season
print('Extract winter temperatures')
Y <- as.4seasons(Tx)  # Y = predictand - Station data
Y <- subset(Y,it=it)
X <- as.4seasons(t2m) # X = predictor - Reanalysis
X <- EOF(subset(X,it=it))

## Do the PCA
pca <- PCA(Y,n=4,verbose=TRUE)

#3 Extract one PC for probing
pc1 <- zoo(pca[,1])

## Do the downscaling with reduced ensemble:
zpca <- DS(pca,X,eofs=1:2,verbose=FALSE)

print('------ Diagnostics -----------')

## Compare original PC1 with the downscaled version of PC1:
zpc1 <- zoo(zpca[,1])
plot(pc1,lwd=4,xlim=range(index(zpc1)),
     ylim=range(c(coredata(zpc1),coredata(pc1)),na.rm=TRUE),
     main=paste('leading',toupper(it),'PC'))

lines(zpc1,col='red',lwd=3,lty=2)

## Plot the PCA-based DS-results
dev.new()
plot(zpca)
