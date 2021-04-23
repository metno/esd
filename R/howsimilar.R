#' howsimilar
#' 
#' field-station comparison
#' Extract time series from field objects with the same coordinates as the stations
#' Estimate correlation, RMSE, max and min for the extracted series. Store as 
#' attributes in the extracted station object. Use regression to calibrate to adjust for the
#' fact that a grid-box volume is not the same as a point measurement.
#' Rasmus E. Benestad
#'
#' \code{howsimilar} is used for assessing reanalyses against station data.
#'
#' 
#' @seealso regrid
#' @aliases howsimilar
#' 
#' @param x A \code{field} object representing the reanalysis
#' @param y A \code{station} object
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
#' @export howsimilar

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
    y0 <- y
    y2 <- regrid(x,is=y)
    if (verbose) {str(y);str(y2); print(class(x)); print(class(y))}
    y <- zoo(y); y2 <- zoo(y2)
    xy <- merge(y,y2)
    ns <- dim(y)[2]
    r <- rep(NA,ns); md <- r; rmse <- r
    for (is in 1:ns) { 
      ok <- is.finite(coredata(xy[,is])) & is.finite(coredata(xy[,is+ns]))
      if (sum(ok)>30) { 
        if (regress) {
          cal <- data.frame(x=coredata(xy)[ok,is],y= coredata(xy)[ok,is+ns])
          z1 <- predict(lm(y ~x, data=cal))
          z2 <- cal$x
        } else {
          z1 <- coredata(xy)[ok,is+ns]
          z2 <- coredata(xy)[ok,is]
        }
        nt <- sum(ok)
        r[is] <- cor(z1,z2)
        md[is] <- max( abs( z1 - z2 ) )
        rmse[is] <- sum( (z1 - z2 )^2, na.rm=TRUE)/ nt
      } 
    }
    
    y <- y0; rm(y0)
    attr(y,'cor') <- r
    attr(y,'rmse') <- rmse
    attr(y,'max.diff') <- md
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
  print(loc(y))
  if (verbose) print('Got the test data.')
  z <- howsimilar(x,y,verbose=verbose)
}