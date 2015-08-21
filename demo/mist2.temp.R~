##  Test DSensemble for PCA data.

library(esd)
it <- 'djf'

## Test data - Norwegian temperature series.
## Only download if the data is not stored locally
if (!file.exists('Norway.Tx.rda')) {
  url <- 'http://files.figshare.com/2073466/Norway.Tx.rda'
  download.file(url,'Norway.Tx.rda')
}

load('Norway.Tx.rda')

## Problem with the unit - quick fix.
attr(Tx,'unit') <- rep('deg C',length(unit(Tx)))

## Get the winter season
print('Extract winter temperatures')
X <- as.4seasons(Tx)
X <- subset(X,it=it)

## Extract one station series for probing
y.0 <- zoo(X[,1])

## Do the PCA
pca <- PCA(X,n=4,verbose=TRUE)


## Do the downscaling with reduced ensemble:
zpca <- DSensemble.pca(pca,biascorrect=TRUE,select=1:10,verbose=TRUE)
zobs <- DSensemble.t2m(subset(Tx,is=1),biascorrect=TRUE,select=1:10,verbose=TRUE)

## Extract one PC for probing
plot(zoo(pca[,1]),lwd=3,col='grey')
lines(zoo(zpca[[2]][,1]),lwd=2)

zpc1 <- zoo(zpca[[12]])
plot(zpc1,lwd=3,xlim=range(index(zpc1)),
     ylim=range(c(coredata(zpc1),coredata(zobs)),na.rm=TRUE))

z <- zoo(zobs)
lines(z,col='blue',lty=2)

lines(zoo(subset(Tx,is=1)))
