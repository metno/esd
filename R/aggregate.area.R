#' Aggregate Functions for Climate Data
#' 
#' These aggregation functions are based on the S3 method for \code{zoo} objects
#' and handle additional tasks such as maintaining attributes with metadata.
#'
#' The function \code{aggregate.area} is specifically designed for aggregating
#' spatial statistics, such as calculating the global mean or the total area of a phenomenon.
#' 
#' The function \code{aggregateArea} is an alias for \code{aggregate.area}.
#' 
#' @seealso \code{\link{aggregate}}, \code{\link{aggregate.size}}
#' @aliases aggregateArea
#' 
#' @param x A \code{\link{station}} object or similar climate data object.
#' @param it A list or data.frame specifying the time index, see \code{\link{subset}}.
#' @param is A list or data.frame specifying the spatial index, see \code{\link{subset}}.
#' @param FUN A function to apply, such as 'sum' or 'mean'.
#' @param na.rm Logical; if TRUE, missing values (\code{NA}) are ignored, see \code{\link{mean}}.
#' @param smallx Logical; if TRUE, reduces the size of the data, default is FALSE.
#' @param verbose Logical; if TRUE, prints detailed progress messages for debugging.
#' @param a Radius of the Earth (in km), used for spatial weighting.
#' @param threshold Numeric; a threshold value used if \code{FUN} is 'area', 'exceedance', or 'lessthan'.
#' @param \dots Additional arguments passed to the function.
#'
#' @return Returns a \code{station} object or a similar climate data object with aggregated results.
#'
#' @author 
#' R.E. Benestad
#'
#' @keywords utilities
#' @examples
#' 
#' ## Aggregation for a station object
#' data(Svalbard)
#' x <- aggregate(Svalbard, month, FUN = 'mean', na.rm = TRUE)
#' plot(x)
#'
#' ## Aggregation for a field object
#' slp <- slp.DNMI()
#' y <- aggregate(slp, year, FUN = 'mean', na.rm = FALSE)
#'
#' ## Spatial aggregation
#' w <- aggregate.area(y)
#' plot(w)
#'
#' @export aggregate.area

aggregate.area <- function(x,...,is=NULL,it=NULL,FUN='sum',
                           na.rm=TRUE,smallx=FALSE,verbose=FALSE,
                           a=6378, threshold=NULL) {
  if (verbose) message("Entering aggregate.area")
  # Call the aggregateArea function with the provided arguments
  y <- aggregateArea(x,...,is=is,it=it,FUN=FUN,
                           na.rm=na.rm,smallx=smallx,verbose=verbose,
                           a=a, threshold=threshold)

  # Ensure indices align if lengths match
  if (length(index(y)) == length(index(x)) && all(index(x) %in% index(y))) {
    index(y) <- index(x)
  } else if (verbose) {
    message("Warning: Index lengths do not match.")
  }
  
  return(y)
}
 
aggregateArea <- function(x,...,is=NULL,it=NULL,FUN='sum',
                           na.rm=TRUE,smallx=FALSE,verbose=FALSE,
                           a=6378, threshold=NULL) {
  # Estimate the area-aggregated values, e.g. the global mean (default)
  if (verbose) {
    message("Entering aggregateArea with function: ", FUN)
  }
  ## REB 20201-02-15: fix the class
  cls0 <- class(x)
  if (verbose) print("Class of x is: ",cls0)
  
  if (verbose) {
    if(is.character(FUN)) {
      if (FUN=='sum') print(rowSums(coredata(x),na.rm=TRUE)) else
                    print(rowMeans(coredata(x),na.rm=TRUE))
    } else print(rowMeans(coredata(x),na.rm=TRUE))
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
    if(is.character(FUN)) {
      if (FUN=='sum') print(rowSums(coredata(x),na.rm=TRUE)) else
                    print(rowMeans(coredata(x),na.rm=TRUE))
    } else print(rowMeans(coredata(x),na.rm=TRUE))
  }
  if (!is.null(attr(x, 'dimensions'))) {
    d <- attr(x, 'dimensions')
  } else {
    d <- c(dim(x), 1)
  }
  if (verbose) message("Data dimensions: ", paste(d, collapse = "-"))
  if (inherits(x,'pattern')) {
    if (verbose) print('need to make the pattern look like field')
    dim(x) <- c(d[1]*d[2],d[3])
    x <- t(x)
  }
  srtlat <- order(rep(lat(x),d[1]))
  dY <- a*diff(pi*lat(x)/180)[1]
  dtheta <- diff(pi*lon(x)/180)[1]
  aweights <- rep(dY * dtheta * a*cos(pi*lat(x)/180),d[1])[srtlat]
  #}
  if (verbose) {
    message("Sum of weights: ", round(sum(aweights, na.rm = TRUE), 4))
  }
  if(is.character(FUN)) if (FUN=='mean') {
    aweights <- aweights/sum(aweights,na.rm=TRUE)
    FUN <- 'sum'
  }
  if (verbose) print(paste('Sum of aweights should be area or 1:',round(sum(aweights))))
  
  ## REB: For sum, we also need to consider the area:
  if(is.character(FUN)) {
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
      aweights2[!is.finite(coredata(X[i,]))] <- NA
      aweights2 <- aweights2/sum(aweights2,na.rm=TRUE)
      X[i,] <- X[i,]*aweights2
    }
    y <- zoo(apply(X,1,FUN,na.rm=na.rm),order.by=index(X))
  }
  nok <- apply(X,1,nv)
  y[nok==0] <- NA
  if (verbose) print(y)
  if(is.character(FUN)) {
    loc <- paste('area',FUN,'of',src(x)) 
    longname <- paste(FUN, attr(x,'longname'))
    method <- paste(FUN,attr(x,'method'))
  } else {
    loc <- paste('area aggregate of',src(x))
    longname <- paste("aggregated", attr(x,'longname'))
    method <- paste("area aggregated",attr(x,'method'))
  }
  Y <- as.station(y,loc=loc,#paste('area',FUN,'of',src(x)),
                  param=attr(x,'variable'),
                  unit=attr(x,'unit'),
                  lon=range(lon(x)),lat=range(lat(x)),alt=NA,cntr=NA,
                  longname=longname,
                  stid=NA,quality=NA,
                  src=attr(x,'source'),url=attr(x,'url'),
                  reference=attr(x,'reference'),info=attr(x,'info'),
                  method=method,
                  type='area aggregate',
                  aspect=attr(x,'aspect'))
  
  attr(Y,'aweights') <- aweights
  ## REB 20201-02-15: fix the class
  class(Y) <- c('station',cls0[-1])
  attr(Y,'history') <- history.stamp(x)
  return(Y)
}

