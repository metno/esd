#' Extention of \code{\link{approx}} for longitude and latitude
#' 
#' Linearly interpolate longitudes and latitudes to n equally spaced points spanning 
#' the interval (min(lon), max(lon)) and (min(lat), max(lat)).
#' 
#' @param lon longitudes
#' @param lat latitudes
#' @param n length of output
#' @param a the radius of the earth (unit: m)
#' @param verbose a boolean; if TRUE print information about progress
#'
#' @export 
approx.lonlat <- function(lon,lat,n=20,a=6.378e06,verbose=FALSE) {
  if (verbose) print("approx.lonlat")
  x <- a * cos( lat*pi/180 ) * cos( lon*pi/180 )
  y <- a * cos( lat*pi/180 ) * sin( lon*pi/180 )
  z <- a * sin( lat*pi/180 )
  xa <- approx(x,n=n)$y
  ya <- approx(y,n=n)$y
  za <- approx(z,n=n)$y
  lona <- atan2( ya, xa )*180/pi
  lata <- asin( za/sqrt( xa^2 + ya^2 + za^2 ))*180/pi
  invisible(cbind(lona,lata))
}