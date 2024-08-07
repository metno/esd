#' \code{nearest} selects the time series in \code{x} that are closest to the locations specified in \code{is}
#'
#' @param x a \code{station} or \code{field} object (if x is a field it will first be transformed to a station object)
#' @param is a spatial index, either a list or data frame containing coordinates (e.g., list(lon=c(...), lat=c(...))) or
#'        a \code{station} or \code{field} object whos longitudes and latitudes will be used
#'
#' @aliases nearest nearest.station nearest.field
#'
#' @export
nearest <- function(x,is) UseMethod("nearest")

#' @exportS3Method
#' @export nearest.station
nearest.station <- function(x,is) {
  
  mind <- function(x) {
    # Minimum distance returns the index pointin to the closest location
    n <- length(x)/2
    x1 <- x[1];   xn <- x[2:n]
    y1 <- x[n+1]; yn <- x[(n+2):(2*n)]
    #print(c(x1,y1,n,length(x),length(yn)))
    
    d <- distAB(x1,y1,xn,yn)
    #print(summary(d))
    i <- (1:length(d))[d == min(d,na.rm=TRUE)][1]
    #print(d[i])
    i
  }

  # The old coordinates
  xo <- lon(x); yo <- lat(x)
  # The wanted locations/new coordinates:
  if (inherits(is,c('station','field'))) {
    xn <- lon(is); yn <- lat(is)
  } else if ( (is.data.frame(is)) | (is.list(is)) ) {xn <- is[[1]]; yn <- is[[2]]}
  # Make matrices containing the coordinates:
  Xo <- matrix(rep(xo,length(xn)),length(xo),length(xn))
  Yo <- matrix(rep(yo,length(yn)),length(yo),length(yn))
  Xo <- rbind(xn,Xo); Yo <- rbind(yn,Yo)
  X <- rbind(Xo,Yo)
  i <- apply(X,2,FUN=mind)
  
  y <- subset(x,is=i)
  return(y)
}

#' @exportS3Method
#' @export nearest.field
nearest.field <- function(x,is) {
  # For field, turn the object into a station class:
  xo <- rep(lon(x),length(lat(x)))
  yo <- sort(rep(lat(x),length(lon(x))))
  y <- x
  ##lon(y) <- xo; lat(y) <- yo
  class(y) <- c('station',class(x)[-1])
  z <- nearest.station(y,is)
  return(z)
}
