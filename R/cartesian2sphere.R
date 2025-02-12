#' Convert Cartesian coordinates to spherical coordinates 
#' 
#' @param lon A numeric vector of longitudes
#' @param lat A numeric vector of latitudes
#' @param lonR Center longitude of viewing angle
#' @param latR Center latitude of viewing angle
#' @param verbose a boolean; if TRUE print information about progress
#' 
#' @returns A list(X=X, Y=Y, Z=Y) where X is the radial distance, Y is the polar angle, and Z is the azimuthal angle. 
#' When displaying data on a spherical grid, the longitude/latitude coordinate system is replaced with X/Z, i.e., 
#' plot(X, Z) instead of plot(lon, lat).
#' 
#' @export
cartesian2sphere <- function(lon, lat, lonR=NULL, latR=NULL, verbose = FALSE) {
  if(verbose) print("cartesian2sphere")
  if(is.null(lonR)) lonR <- mean(lon, na.rm=TRUE)
  if(is.null(latR)) latR <- mean(lat, na.rm=TRUE)
  
  Theta <- degrees2radians(lon)
  Phi <- degrees2radians(lat)
  
  # Transform -> (X,Y,Z):
  X <- sin(Theta)*cos(Phi)
  Y <- cos(Theta)*cos(Phi)
  Z <- sin(Phi)
  
  # Grid coordinates:
  d <- dim(X)
  
  # Rotate data grid:  
  A <- rotM(x=0, y=0, z=lonR) %*% rbind(c(X), c(Y), c(Z))
  A <- rotM(x=latR, y=0, z=0) %*% A
  X <- A[1,]; Y <- A[2,]; Z <- A[3,]
  if(!is.null(d)) {dim(X) <- d; dim(Y) <- d; dim(Z) <- d}
  return(list(X=X, Y=Y, Z=Z))
}