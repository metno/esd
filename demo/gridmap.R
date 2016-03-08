gridmap <- function(Y,FUN='mean',colbar=NULL,project='lonlat',xlim=NULL,ylim=NULL,verbose=FALSE) {
  if (verbose) print(paste('gridmap',FUN))
  require(LatticeKrig)

  if (is.null(xlim)) xlim <- range(lon(Y))
  if (is.null(ylim)) ylim <- range(lat(Y))
  if (!is.null(dim(y)))
      y <- apply(Y,2,FUN,na.rm=TRUE)
  else
     y <- Y  ## single specific date

  ## Get data on the topography on the 5-minute resolution
  if (verbose) print('Use etopo5 elevation data')
  data(etopo5)
  etopo5 <- subset(etopo5,
                   is=list(lon=range(lon(Y))+c(-1,1),
                           lat=range(lat(Y))+c(-1,1)))
  ## Mask the sea: elevations below 1m below sea level is masked.
  etopo5[etopo5<=-1] <- NA

  ## Set the grid to be the same as that of etopo5:
  if (verbose) print('Use same structure as etopo5')
  grid <- structure(list(x=lon(etopo5),y=lat(etopo5)),class='gridList')

  ## Flag dubplicated stations:
  if (verbose) print('Check for duplicates')
  ok <- !(duplicated(lon(Y)) & duplicated(lat(Y)))

  ## Kriging
  if (verbose) print('Apply kriging')
  browser()
  obj <- LatticeKrig( x=cbind(lon(Y)[ok],lat(Y)[ok]),
                      y=y[ok],Z=alt(Y)[ok])

  ##  obj <- LatticeKrig( x=cbind(lon[ok],lat[ok]), y=z[2,ok],Z=alt[ok])
  if (verbose) print('Predict surface')
  w <- predictSurface(obj, grid.list = grid,Z=etopo5)
  w$z[is.na(etopo5)] <- NA

## Get rid of packages that have functions of same name:
  detach("package:LatticeKrig")
  detach("package:fields")
  detach("package:spam")
  detach("package:grid")
  detach("package:maps")

  ## Convert the results from LatticeKrig to esd:
  W <- w$z
  attr(W,'variable') <- varid(Y)[1]
  attr(W,'unit') <- unit(Y)[1]
  attr(W,'longitude') <- w$x
  attr(W,'latitude') <- w$y
  class(W) <- class(etopo5)

  ## Make a projection that zooms in on the Barents region

  map(W,xlim=xlim,ylim=ylim,
      colbar=colbar,project=project)
  invisible(W)
}
