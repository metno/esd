#' Normalise data
#'
#' \code{as.stand} returns normalised values:
#' If the input contains precipitation data the data are normalised by the mean value.
#' If the input contains temperature data, the data are stanardised by subtracting the mean and dividing by the standard deviation.
#' \code{as.original} transforms normalised data to its original values.
#'
#' @aliases as.stand as.stand.station as.original as.original.station
#'
#' @param x a station object
#' @param verbose a boolean; if TRUE print information about progress
#' @param na.rm a boolean; if TRUE remove NA values
#'
#' @export
as.stand <- function(x,...) UseMethod("as.stand")

#' @export
as.stand.station <- function(x,...,verbose=FALSE,na.rm=TRUE) {
  if(verbose) print("as.stand.station")
  if (is.precip(x)) {
    mu <- apply(x,2,mean,na.rm=na.rm)
    X <- 100*x/mu
    attr(X,'clim') <- mu
    attr(X,'aspect') <- 'proportional'
    attr(X,'unit') <- '%'
    attr(X,'oldunit') <- attr(x,'unit')
  } else if (is.T(x)) {
    mu <- apply(x,2,mean,na.rm=na.rm)
    sigma <- apply(x,2,sd,na.rm=na.rm)
    X <- (x - mu)/sigma
    attr(X,'mean') <- mu
    attr(X,'sigma') <- sigma
    attr(X,'aspect') <- 'standardised'    
  }
  attr(X,'history') <- history.stamp(x)
  return(X)
}

#' @export
as.original <- function(x) UseMethod("as.original")

#' @export
as.original.station <- function(x) {
  if (attr(x,'aspect')=='proportional') {
    X <- attr(x,'clim')*x/100
    attr(X,'clim') <- NULL
    attr(X,'unit') <- attr(x,'oldunit')
    attr(X,'oldunit') <- NULL
    attr(X,'aspect') <- 'original'
  } else if (attr(x,'aspect')=='standardised') {
    X <- x * attr(x,'sigma') + attr(x,'mean')
    attr(X,'mean') <- NULL
    attr(X,'sigma') <- NULL
    attr(X,'aspect') <- 'original'     
  } else X <- x
  attr(X,'history') <- history.stamp(x)
  return(X)
}
