## field-station comparison
## Extract time series from field objects with the same coordinates as the stations
## Estimate correlation, RMSE, max and min for the extracted series. Store as 
## attributes in the extracted station object. Use regression to calibrate to adjust for the
## fact that a grid-box volume is not the same as a point measurement.
## Rasmus E. Benestad
#' aggregate
#' 
#' The aggregation functions are based on the S3 method for \code{zoo} objects,
#' but takes care of extra house keeping, such as attributes with meta data.
#'
#' \code{howsimilar} is used for assessing reanalyses against station data.
#'
#' 
#' @seealso regrid
#' @aliases howsimilar
#' 
#' @param x A \code{\link{field}} object representing the reanalysis
#' @param y A \code{\link{station}} object
#' @param plot a boolean; if TRUE ignore NA, see \code{\link{mean}}
#' @param regress a boolean defaulting to FALSE
#' @param verbose a boolean; if TRUE print information about progress
#'
#' @return The call returns a station object
#'
#' @author R.E. Benestad
#' @keywords utilities
#' @examples
#' test.howsimilar()
#'
#' @export 

howsimilar <- function(x,y,plot='cor',FUN=NULL,regress=FALSE,verbose=FALSE) {
  if (verbose) {print('howsimilar'); print(class(x)); print(class(y))}
  if (is.null(FUN) & !is.precip(x)) FUN <- 'mean' else
    if (is.null(FUN) & is.precip(x)) FUN <- 'sum'
    if (class(x)[2] != class(y)[2]) {
      if (verbose) print('Ensure same timescales')
      if (class(x)[2]=='annual') y <- annual(y,FUN=FUN)
      if (class(x)[2]=='season') y <- as.4seasons(y,FUN=FUN)
      if (class(x)[2]=='month') y <- as.monthly(y,FUN=FUN)
      if (class(x)[2]=='day') { 
        if (class(y)[2]=='annual') x <- annual(x,FUN=FUN)
        if (class(y)[2]=='season') x <- as.4seasons(x,FUN=FUN)
        if (class(y)[2]=='month') x <- as.monthly(x,FUN=FUN)
      }
    }
    y <- subset(y,is=list(lon=range(lon(x)),lat=range(lat(x))))
    y2 <- regrid(x,is=y)
    if (verbose) {str(y);str(y2); print(class(x)); print(class(y))}
    y <- zoo(y); y2 <- zoo(y2)
    xy <- merge(y,y2)
    nt <- length(index(xy))
    ns <- dim(y)[2]
    r <- rep(NA,ns); md <- r; rmse <- r
    for (is in 1:ns) { 
      ok <- is.finite(coredata(xy[,i])) & is.finite(coredata(xy[,i+ns]))
      if (sum(ok)>30) { 
        r[is] <- cor( coredata(xy)[ok,i], coredata(xy)[ok,i+ns] )
        md[is] <- max( abs( coredata(xy[ok,i]) - coredata(xy[ok,i+ns]) ) )
        rmse[is] <- sum( (coredata(xy[ok,i]) - coredata(xy[ok,i+ns]) )^2, na.rm=TRUE)/nt
      } else browser()
    }
    
    attr(y,'cor') <- r
    attr(y,'rmse') <- rmse
    attr(y,'max.diff') <- md
    browser()
    if (is.character(plot)) map(y,FUN=plot)
    invisible(y)
}

#' @export
test.howsimilar <- function(x=NULL,y=NULL,verbose=FALSE) {
  if (verbose) print('test.howsimilar')
  if (is.null(x)) x <- annual(t2m.DNMI())
  if (is.null(y)) y <- annual(station(param='t2m',src='nacd'))
  nv <- apply(y,2,'nv')
  y <- subset(y,is= (nv >= 60))
  print(names(y))
  if (verbose) print('Got the test data.')
  z <- howsimilar(x,y,verbose=verbose)
}