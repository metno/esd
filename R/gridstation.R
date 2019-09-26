#' Function for gridding station data Y.
#'
#' \code{gridstation} transforms a \code{station} object into a \code{field} object
#' by interpolation, using the package 'LatticeKrig'.
#'
#' @param Y a \code{station} object
#' @param i index
#' @param verbose a boolean; if TRUE print information about progress 
#' @param xlim longitude range
#' @param ylim latitude range
#'
#' @return a \code{field} object
#'
#' @export
gridstation <- function(Y,i=1,verbose=FALSE,xlim=NULL,ylim=NULL) {
  if (verbose) print(paste('gridstation'))
  if (!requireNamespace("LatticeKrig", quietly = TRUE)) {
    stop("Package 'LatticeKrig' needed to use 'gridstation'. Please install it.")
  } else {
    if (is.null(xlim)) xlim <- range(lon(Y))
    if (is.null(ylim)) ylim <- range(lat(Y))
    
    ## Get data on the topography on the 5-minute resolution
    data(etopo5, envir = environment())
    etopo5 <- subset(etopo5, is=list(lon=range(lon(Y))+c(-1,1),
                             lat=range(lat(Y))+c(-1,1)))
    ## Mask the sea: elevations below 1m below sea level is masked.
    etopo5[etopo5<=-1] <- NA
    
    ## Set the grid to be the same as that of etopo5:
    grid <- structure(list(x=lon(etopo5),y=lat(etopo5)),class='gridList')
    
    ## Flag dubplicated stations:
    ok <- !(duplicated(lon(Y)) & duplicated(lat(Y)))
    
    if (verbose) print(paste('Use',sum(ok),'locations from',length(lon(Y))))
    
    obj <- LatticeKrig::LatticeKrig( x=cbind(lon(Y)[ok],lat(Y)[ok]),
                        y=Y[ok,i],Z=alt(Y)[ok])
    
    w <- fields::predictSurface(obj, grid.list = grid,Z=etopo5)
    w$z[is.na(etopo5)] <- NA
    
    ## Convert the results from LatticeKrig to esd:
    W <- w$z
    attr(W,'variable') <- varid(Y)[1]
    attr(W,'unit') <- unit(Y)[1]
    attr(W,'longitude') <- w$x
    attr(W,'latitude') <- w$y
    
    ## Make a projection that zooms in on the Barents region
    
    invisible(W)
  }
}