#' aggregate.grid
#' 
#' The aggregation functions are based on the S3 method for \code{zoo} objects,
#' but takes care of extra house keeping, such as attributes with meta data.
#' Use given longitude-latitude coordinates to find the grid-boxes which are within these bins and
#' aggregate them as with a data.frame object. This is a bit similar to regrid, but does not 
#' use bilinear interpolation. but the mean of grid boxes within a larger grid-box. 
#' Rasmus E. Benestad
#'
#' \code{aggregate.grid} is used for aggregating spatial statistics of a finer grid-mesh onto a coarser grid.
#'
#' 
#' @seealso aggregate.area regrid
#' @aliases aggregate.grid
#' 
#' @param x A \code{\link{field}} object with finer grid mesh.
#' @param is A list or \code{\link{field}} with the coarser grid with aggregated data from finer mesh.
#' @param FUN A function, e.g., 'sum' or 'mean'
#' @param na.rm a boolean; if TRUE ignore NA, see \code{\link{mean}}
#' @param verbose a boolean; if TRUE print information about progress
#' @param \dots additional arguments
#'
#' @return The call returns a field object
#'
#' @author R.E. Benestad
#' @keywords utilities
#' @examples
#' 
#' test.aggregate.grid()
#' 
#' @export aggregate.grid
aggregate.grid <- function(x,...,is,FUN='mean',na.rm=TRUE,verbose=FALSE) {
  ## The coordinates of the new aggregated results
  if (verbose) print('aggregate.grid')
  if (is.field(is)) {
    lons <- lon(is)
    lats <- lat(is)
  } else if (is.list(is)) {
    lons <- is[[1]]
    lats <- is[[2]]
  }
  ## The oroginal coordinates
  Lons <- lon(x)
  Lats <- lat(x)
  
  
  ## Organise the indices
  dx <- diff(lons)[1]
  ix <- (lons - min(lons))/dx
  dy <- diff(lats)[1]
  iy <- (lats - min(lats))/dy
  dX <- diff(Lons)[1]
  iX <- (Lons - min(lons))/dx
  dY <- diff(Lats)[1]
  iY <- (Lats - min(lats))/dy
  
  ## length of lons and lats: the dimensions of the resulting aggregated field:
  ny <- length(lats); nY <- length(Lats)
  nx <- length(lons); nX <- length(Lons)
  if (verbose) print(c(nx,ny,nX,nY,min(ix),max(ix),min(iy),max(iy),min(iX),max(iX),min(iY),max(iY)))
  
  xy <- paste(rep(round(ix),ny),sort(rep(round(iy),nx)))
  XY <- paste(rep(round(iX),nY),sort(rep(round(iY),nX)))
  if (verbose) {print(xy); print(table(XY))}
  if (sum(!is.element(XY,xy))!=0) { 
    XY[!is.element(XY,xy)] <- NA
    #print(XY[!is.element(XY,xy)])
  }
  nt <- length(index(x))
  z <- matrix(rep(NA,nx*ny*nt),nt,nx*ny)
  for (it in 1:nt) {
    zzz <- data.frame(x=c(coredata(x)[it,]))
    ZZZ <- aggregate(zzz,by=list(XY),FUN=FUN, na.rm=na.rm, ...)
    z[it,match(ZZZ$Group.1,xy)] <- ZZZ$x
    if (verbose) cat('.')
  } 
  z <- zoo(x=z,order.by=index(x))
  z <- as.field(z,lon=lons,lat=lats,param=varid(x),unit=esd::unit(x))
  attr(z,'history') <- history.stamp()
  class(z) <- class(x)
  invisible(z)
}

#' @export test.aggregate.grid
test.aggregate.grid <- function(x=NULL,verbose=FALSE) {
  print("Test aggregate.grid")
  if (is.null(x)) {
    print('Use default: t2m.DNMI')
    x <- t2m.DNMI()
  }
  map(x,main='Original')
  y <- regrid(x,is=list(lon=seq(min(lon(x)),max(lon(x)),by=1),lat=seq(min(lat(x)),max(lat(x)),by=1)))
  map(y,main='High-res version')
  print('aggregate gridboxes...')
  z <- aggregate.grid(y,is=x,verbose=verbose)
  print('Aggregate.grid finished')
  map(z,main='Grid-aggegated results')
  print('***')
}