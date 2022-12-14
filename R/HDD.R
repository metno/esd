#' Heating, growing and cooling degree days
#' 
#' Functions to estimate \code{HDD}, \code{GDD}, and \code{CDD}. 
#'
#' @aliases HDD GDD CDD
#' @seealso hotsummerdays coldwinterdays coldspells heatwavespells nwetdays plot count spell
#'
#' @importFrom stats glm qqline
#'
#' @param x station or field object
#' @param threshold threshold value
#' @param na.rm TRUE - remove NAs.

#'
#' @return Station or field objects
#'
#' @keywords utilities
#' @examples
#' 
#' # Growing degree days:
#' data(ferder)
#' plot(as.seasons(ferder,FUN='GDD'), new=FALSE)
#'
#' @export
HDD <- function(x,threshold=18,na.rm=TRUE) {
  cold <- x < threshold
  hdd <- sum(threshold - x[cold],na.rm=na.rm)
  return(hdd)
}


#' @export
CDD <- function(x,threshold=22,na.rm=TRUE) {
  warm <- x > threshold
  cdd <- sum(x[warm] - threshold,na.rm=na.rm)
  return(cdd)
}

#' @export
GDD <- function(x,threshold=10,na.rm=TRUE) {
  gdd <- CDD(x,threshold=threshold)
  attr(gdd,'url') <- 'http://en.wikipedia.org/wiki/Growing_degree-day'
  return(gdd)
}
