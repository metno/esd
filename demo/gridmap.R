gridmap <- function(Y,breaks=NULL,pal=NULL,verbose=FALSE) {
  require(LatticeKrig)

  y <- apply(annual(Y,FUN='sum'),2,'mean',na.rm=TRUE)

  ## Get data on the topography on the 5-minute resolution
  data(etopo5)
  etopo5 <- subset(etopo5,
                   is=list(lon=range(lon(Y))+c(-1,1),
                           lat=range(lat(Y))+c(-1,1)))
  ## Mask the sea: elevations below 1m below sea level is masked.
  etopo5[etopo5<=-1] <- NA

  ## Set the grid to be the same as that of etopo5:
  grid <- structure(list(x=lon(etopo5),y=lat(etopo5)),class='gridList')

  ## Flag dubplicated stations:
  ok <- !(duplicated(lon(Y)) & duplicated(lat(Y)))

  ## Spread in the  90-percente interval changing
  obj <- LatticeKrig( x=cbind(lon(Y)[ok],lat(Y)[ok]),
                      y=y[ok],Z=alt(Y)[ok])

  ##  obj <- LatticeKrig( x=cbind(lon[ok],lat[ok]), y=z[2,ok],Z=alt[ok])
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

  rev <- switch(varid(Y)[1],'t2m'=FALSE,'precip'=TRUE)
  Wx <- max(abs(W),na.rm=TRUE)
  if (is.null(breaks)) breaks <- round(seq(-Wx,Wx,length=31),2) 
  if (is.null(pal)) pal <- varid(Y)[1]
  map(W,xlim=range(lon(W)),ylim=range(lat(W)),
      colbar=list(pal=pal))
  invisible(W)
}
