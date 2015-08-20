##  Run DSensemble for PCA data for monthly temperature in Norway.

library(esd)

## Test data - Norwegian temperature series.
## Only download if the data is not stored locally
if (!file.exists('Norway.Tx.rda')) {
  url <- 'http://files.figshare.com/2073466/Norway.Tx.rda'
  download.file(url,'Norway.Tx.rda')
}

dswork <- function(Tx,it='djf') {

## Get the winter season
  print(paste('Extract monthly temperatures for ',it))
  X <- as.monthly(Tx)
  X <- subset(X,it=it)

  ## Do the PCA
  pca <- PCA(X)

  ## Do the downscaling with reduced ensemble:
  zpca <- DSensemble.pca(pca,biascorrect=TRUE)

  save(file=paste('mist2.temp.',it,'.rda',sep=''),zpca)
}


load('Norway.Tx.rda')
## Problem with the unit - quick fix.
attr(Tx,'unit') <- rep('deg C',length(unit(Tx)))

for (it in month.abb) dswork(Tx,it)
