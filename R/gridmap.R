#' Creates a griddded map %% ~~function to do ... ~~
#' 
#' A function that uses \code{LatticeKrieg} and elevation data to grid station
#' based data and present a map.
#' 
#' 
#' @param Y A station object
#' @param FUN A function or name of a function, e.g, "mean" or "trend"
#' @param colbar A list specifying the color bar, e.g., list(col="precip",
#' breaks=seq(1,10), rev=FALSE)
#' @param project projection: "lonlat" or "sphere"
#' @param xlim range of x-axis
#' @param ylim range of y-axis
#' @param zlim range of color axis
#' @param verbose if TRUE print information about progress
#' @param plot if TRUE display results as plots
#' @examples
#' 
#' data("precip.NORDKLIM")
#' precip.gp <- gridmap(precip.NORDKLIM, plot=TRUE)
#' 
#' @export gridmap
gridmap <- function(Y,FUN='mean',colbar=list(pal='t2m'),project='lonlat',xlim=NULL,ylim=NULL,zlim=NULL,verbose=FALSE,plot=FALSE,new=TRUE) {

  if (verbose) print(paste('gridmap',FUN))
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
    w <- fields::predictSurface(obj, grid.list = grid,Z=etopo5)
    w$z[is.na(etopo5)] <- NA

    ## Convert the results from LatticeKrig to esd:
    W <- w$z
    attr(W,'variable') <- varid(Y)[1]
    attr(W,'unit') <- unit(Y)[1]
    attr(W,'longitude') <- w$x
    attr(W,'latitude') <- w$y
    class(W) <- class(etopo5)
  
    ## Make the graphics
    if(plot) {
      if (verbose) print("make the map")
      map(W,xlim=xlim,ylim=ylim,zlim=zlim,colbar=colbar,project=project,new=new)
    }
    invisible(W)
  }
}
