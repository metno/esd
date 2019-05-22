#' aggregate
#' 
#' The aggregation functions are based on the S3 method for \code{zoo} objects,
#' but takes care of extra house keeping, such as attributes with meta data.
#'
#' \code{aggregate.area} is used for aggregating spatial statistics, such as
#' the global mean or the global area of some phenomenon.
#' 
#' @seealso aggregate aggregate.size
#' 
#' @param x A \code{\link{station}} object
#' @param it A list or data.frame providing time index, see \code{\link{subset}}
#' @param is A list or data.frame providing space index, see \code{\link{subset}}
#' @param FUN A function, e.g., 'sum' or 'mean'
#' @param na.rm a boolean; if TRUE ignore NA, see \code{\link{mean}}
#' @param smallx a boolean defaulting to FALSE
#' @param verbose a boolean; if TRUE print information about progress
#' @param a radius of earth (unit: km)
#' @param threshold threshold to be used if FUN is 'area','exceedance', or 'lessthan'
#'
#' @return The call returns a station object
#'
#' @author R.E. Benestad
#' @keywords utilities
#' @examples
#' 
#' ## S3 method for class 'station'
#' data(Svalbard)
#' x <- aggregate(Svalbard, month, FUN='mean', na.rm=TRUE)
#' plot(x)
#'
#' ## S3 method for class 'field'
#' slp <- slp.DNMI()
#' y <- aggregate(slp, year, FUN='mean', na.rm=FALSE)
#'
#' ## Aggregate area
#' w <- aggregate.area(y)
#' plot(w)
#'
#' @export
aggregate.area <- function(x,is=NULL,it=NULL,FUN='sum',
                           na.rm=TRUE,smallx=FALSE,verbose=FALSE,
                           a=6378, threshold=NULL) {
  # Estimate the area-aggregated values, e.g. the global mean (default)
  if (verbose) print(paste("aggregate.area",FUN))
  if (verbose) {
    if (FUN=='sum') print(rowSums(coredata(x),na.rm=TRUE)) else
                    print(rowMeans(coredata(x),na.rm=TRUE))
  }
  if (inherits(x,'eof')) {
    if (verbose) print('aggregate.area for EOF')
    y <- as.pattern(x)
    ya <- aggregate.area(y,is=is,FUN=FUN,na.rm=na.rm,smallx=smallx,verbose=verbose,
                         a=a,threshold=threshold)
    if (verbose) {print(length(ya)); print(length(attr(x,'eigenvalues'))); print(t(dim(coredata(x))))}
    z <- apply(diag(ya*attr(x,'eigenvalues')) %*% t(coredata(x)),2,FUN='sum')
    if (is.zoo(x)) z <- zoo(x=z,order.by=index(x))
    attr(z,'history') <- history.stamp(x)
    return(z)
  }
  x <- subset(x,is=is,it=it,verbose=verbose)
  if ( (verbose) & (!is.null(is) | !is.null(it)) ) {
    if (FUN=='sum') print(rowSums(coredata(x),na.rm=TRUE)) else
                    print(rowMeans(coredata(x),na.rm=TRUE))
  }
  if (inherits(FUN,'function')) FUN <- deparse(substitute(FUN)) # REB140314
  if (!is.null(attr(x,'dimensions'))) d <- attr(x,'dimensions') else d <- c(dim(x),1)
  if (verbose) print(paste('dimensions',paste(d,collapse='-')))
  if (inherits(x,'pattern')) {
    if (verbose) print('need to make the pattern look like field')
    dim(x) <- c(d[1]*d[2],d[3])
    x <- t(x)
  }

  srtlat <- order(rep(lat(x),d[1]))
  dY <- a*diff(pi*lat(x)/180)[1]
  dtheta <- diff(pi*lon(x)/180)[1]
  ## The first assumes a global field and the second is for a limited longitude range
  #if (diff(range(lon(x)))> 350) {
  #  aweights <- rep(dY * 2*pi/d[1] * a*cos(pi*lat(x)/180),d[1])[srtlat]
  #} else {
  aweights <- rep(dY * dtheta * a*cos(pi*lat(x)/180),d[1])[srtlat]
  #}
  if (verbose) print(sum(aweights))
  if (FUN=='mean') {
    aweights <- aweights/sum(aweights,na.rm=TRUE)
    FUN <- 'sum'
  }
  if (verbose) print(paste('Sum of aweights should be area or 1:',round(sum(aweights))))
  
  ## REB: For sum, we also need to consider the area:
  if (FUN %in% c('sum','area','exceedance','exceedence','lessthan')) {
    if (FUN=='area') {
      ## Estimate the area of the grid boxes
      coredata(x) -> cx
      if (is.null(threshold)) {
        cx[is.finite(cx)] <- 1; cx[!is.finite(cx)] <- 0
      } else {
        cx[cx<threshold] <- 0; cx[cx >= threshold] <- 1
      }
      coredata(x) <- cx; rm('cx'); gc(reset=TRUE)
      FUN <- 'sum'
    } else if ( (FUN %in% c('exceedance','exceedence')) & !is.null(threshold) ) {
      # Estimate the sum of grid boxes with higher value than threshold
      coredata(x) -> cx
      cx[cx < threshold] <- NA
      coredata(x) <- cx; rm('cx'); gc(reset=TRUE)
      FUN <- 'sum'
    } else if ( (FUN == 'lessthan') & !is.null(threshold) ) {
      # Estimate the sum of grid boxes with lower value than threshold
      coredata(x) -> cx
      cx[cx >= threshold] <- NA
      coredata(x) <- cx; rm('cx'); gc(reset=TRUE)
      FUN <- 'sum'
    } 
    attr(x,'unit') <- paste(attr(x,'unit'),' * km^2')
  }
  
  if (smallx) {
    if (sum(!is.finite(x))==0) {
      X <- coredata(x)%*%diag(aweights)
    } else {
      if (verbose) print('Need to account for missing data in the area weighting')
      Aweights <- rep(aweights,length(index(x)))
      dim(Aweights) <- dim(x)
      print('This is incomplete - needs checking!')
      Aweights[!is.finite(coredata(x))] <- NA
      Aweights <- Aweights/apply(Aweights,1,FUN='sum',na.rm=TRUE)
      X <- coredata(X)*Aweights
    }
    y <- zoo(apply(X,1,FUN,na.rm=na.rm),order.by=index(x))
  } else {
    X <- coredata(x) 
    if (d[3]==1) dim(X) <- c(1,length(X)) ## If only one map, then set the dimensions right to get a matrix.
    if (verbose) {print(dim(X)); print(length(aweights))}
    for (i in 1:d[3]) {
      ## Temporary weights to account for variable gaps of missing data
      aweights2 <- aweights
      aweights2[!is.finite(coredata(x[i,]))] <- NA
      aweights2 <- aweights2/sum(aweights2,na.rm=TRUE)
      X[i,] <- X[i,]*aweights2
    }
    y <- zoo(apply(X,1,FUN,na.rm=na.rm),order.by=index(x))
  }
  if (verbose) print(y)
  
  Y <- as.station(y,loc=paste('area',FUN,'of',src(x)),
                  param=attr(x,'variable'),
                  unit=attr(x,'unit'),
                  lon=range(lon(x)),lat=range(lat(x)),alt=NA,cntr=NA,
                  longname=paste(FUN,attr(x,'longname')),stid=NA,quality=NA,
                  src=attr(x,'source'),url=attr(x,'url'),
                  reference=attr(x,'reference'),info=attr(x,'info'),
                  method=paste(FUN,attr(x,'method')),type='area aggregate',
                  aspect=attr(x,'aspect'))
  if (verbose) attr(Y,'aweights') <- aweights
  attr(Y,'history') <- history.stamp(x)
  return(Y)
}

