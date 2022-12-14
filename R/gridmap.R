#' Creates a griddded map
#' 
#' A function that uses \code{LatticeKrieg} and elevation data to grid station
#' based data and present a map.
#' @aliases gridmap gridmap.default gridmap.station gridmap.pca gridstations
#'  
#' @param Y A station object or a PCA object. 
#' @param FUN A function or name of a function, e.g, "mean" or "trend"
#' @param colbar A list specifying the color bar, e.g., list(col="precip",
#' breaks=seq(1,10), rev=FALSE)
#' @param project projection: "lonlat" or "sphere"
#' @param xlim range of x-axis
#' @param ylim range of y-axis
#' @param zlim range of color axis
#' @param verbose if TRUE print information about progress
#' @param plot if TRUE display results as plots
#' @param new if TRUE plot in new window
#'
#' @examples
#' 
#' data("precip.NORDKLIM")
#' precip.gp <- gridmap(precip.NORDKLIM, plot=TRUE)
#' map(precip.gp)
#' @export
gridmap <- function(Y,FUN='mean',colbar=list(pal='t2m'),project='lonlat',xlim=NULL,ylim=NULL,
                    zlim=NULL,verbose=FALSE,plot=FALSE,new=TRUE) { UseMethod("gridmap") }

#' @exportS3Method
#' @export 
gridmap.default <- function(Y,FUN='mean',colbar=list(pal='t2m'),project='lonlat',xlim=NULL,ylim=NULL,
                    zlim=NULL,verbose=FALSE,plot=FALSE,new=TRUE) {

  if (verbose) print(paste('gridmap',FUN))
  if (is.null(Y)) {warning('Empty station object'); return(NULL)}
  if (!requireNamespace("LatticeKrig", quietly = TRUE)) {
    stop("Package 'LatticeKrig' needed to use 'gridmap'. Please install it.")
  } else {
    
    if (is.null(xlim)) xlim <- range(lon(Y))
    if (is.null(ylim)) ylim <- range(lat(Y))
    if (!is.null(dim(Y))) {
      y <- apply(Y,2,FUN,na.rm=TRUE)
    } else {
      y <- Y  ## single specific date
    }
  
    ## Get data on the topography on the 5-minute resolution
    if (verbose) print('Use etopo5 elevation data')
    data(etopo5, envir = environment())
    etopo5 <- subset(etopo5,is=list(lon=range(lon(Y))+c(-1,1),
                                    lat=range(lat(Y))+c(-1,1)))
    
    ## Mask the sea: elevations below 1m below sea level is masked.
    etopo5[etopo5<=-1] <- NA
    if (!is.null(zlim)) {etopo5[(etopo5<min(zlim)) | ((etopo5>max(zlim)))] <- NA}

    ## Set the grid to be the same as that of etopo5:
    if (verbose) print('Use same structure as etopo5')
    grid <- structure(list(x=lon(etopo5),y=lat(etopo5)),class='gridList')

    ## Flag duplicated stations:
    if (verbose) print('Check for duplicates')
    ok <- !(duplicated(lon(Y)) & duplicated(lat(Y)))

    ## Kriging
    if (verbose) print(paste('Apply kriging to',sum(ok),'locations'))
    obj <- LatticeKrig::LatticeKrig(x=cbind(lon(Y)[ok],lat(Y)[ok]),
                                    y=y[ok],Z=alt(Y)[ok])

    ##  obj <- LatticeKrig::LatticeKrig( x=cbind(lon[ok],lat[ok]), y=z[2,ok],Z=alt[ok])
    if (verbose) print('Predict surface')
    w <- fields::predictSurface(obj, grid.list = grid, Z=etopo5)
    w$z[is.na(etopo5)] <- NA

    ## Convert the results from LatticeKrig to esd:
    W <- w$z
    attr(W,'variable') <- varid(Y)[1]
    attr(W,'unit') <- esd::unit(Y)[1]
    attr(W,'longitude') <- w$x
    attr(W,'latitude') <- w$y
    class(W) <- class(etopo5)
    ## Make the graphics
    if(plot) {
      if (verbose) print("make the map")
      map(W,xlim=xlim,ylim=ylim,zlim=zlim,colbar=colbar,
          project=project,new=new,verbose=verbose)
    }
    invisible(W)
  }
}

#' @exportS3Method
#' @export 
gridmap.station <- function(Y,FUN='mean',colbar=list(pal='t2m'),project='lonlat',xlim=NULL,ylim=NULL,
                            zlim=NULL,verbose=FALSE,plot=FALSE,new=TRUE) {
  if (verbose) print('gridmap.station')
  ## KMP 2021-05-19: Calling gridmap.station from within itself creates an infinite loop
  #x <- gridmap.station(Y=Y,FUN=FUN,colbar=colbar,project=project,xlim=xlim,ylim=ylim,zlim=zlim,verbose=verbose,plot=plot,new=new)
  x <- gridmap.default(Y=Y,FUN=FUN,colbar=colbar,project=project,xlim=xlim,ylim=ylim,zlim=zlim,verbose=verbose,plot=plot,new=new)
  return(x)
}
## REB 2021-05-07: added a method to grid PCAs through kriging.
#' @exportS3Method
#' @export 
gridmap.pca <- function(Y,FUN='mean',colbar=list(pal='t2m'),project='lonlat',xlim=NULL,ylim=NULL,
                            zlim=NULL,verbose=FALSE,plot=FALSE,new=TRUE) {
  ## Convert a PCA to EOF
  if (verbose) print('gridmap.pca')
  d <- dim(Y)
  for (id in 1:d[2]) { 
    if (verbose) print(paste('PCA pattern',id))
    y <- attr(Y,'pattern')[,id]
    attr(y,'longitude') <- lon(Y)
    attr(y,'latitude') <- lat(Y)
    attr(y,'altitude') <- alt(Y)
    attr(y,'unit') <- 'weight'
    attr(y,'variable') <- paste0(varid(y),'.pca')
    x <- gridmap.default(Y=y,FUN=FUN,colbar=colbar,project=project,xlim=xlim,ylim=ylim,zlim=zlim,verbose=verbose,plot=plot,new=new)
    D <- dim(x)
    if (id==1) X <- c(x) else X <- cbind(X,c(x))
  }
  ## Also grid the mean values
  ym <-  attr(Y,'mean')
  attr(ym,'longitude') <- lon(Y)
  attr(ym,'latitude') <- lat(Y)
  attr(ym,'altitude') <- alt(Y)
  attr(ym,'unit') <- 'weight'
  attr(ym,'variable') <- paste0(varid(y),'.pca')
  xm <- gridmap.default(Y=ym,FUN=FUN,colbar=colbar,project=project,xlim=xlim,ylim=ylim,zlim=zlim,verbose=verbose,plot=plot,new=new)
  
  dim(X) <- c(D,d[2])
  attr(Y,'longitude') <- lon(x)
  attr(Y,'latitude') <- lat(x)
  attr(Y,'pattern') <- X
  attr(Y,'dimensions') <- c(D,d[1])
  attr(Y,'mean') <- xm
  class(Y)[1] <- 'eof'
  return(Y)
}



