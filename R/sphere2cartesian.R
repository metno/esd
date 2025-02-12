#' Convert rotated spherical coordinates back to Cartesian coordinates
#'
#' @param X A numeric vector of X coordinates, where X is the radial distance
#' @param Y A numeric vector of Y coordinates, where Y is the polar angle
#' @param Z A numeric vector of Z coordinates, where Z is the azimuthal angle
#' @param lonR Center longitude of viewing angle
#' @param latR Center latitude of viewing angle
#' @param verbose a boolean; if TRUE print information about progress
#'
#' @returns A list(Lon=longitudes, Lat=latitudes)
#' 
#' @export
sphere2cartesian <- function(X, Y, Z, lonR, latR, verbose = FALSE) {
  if (verbose) print("sphere2cartesian")
  
  # Convert rotation angles to radians
  lonR_rad <- degrees2radians(lonR)
  latR_rad <- degrees2radians(latR)
  
  # Inverse rotation to undo the viewing angle rotation
  inv_rotation_matrix_lat <- solve(rotM(x = latR, y = 0, z = 0))
  inv_rotation_matrix_lon <- solve(rotM(x = 0, y = 0, z = lonR))
  
  # Apply inverse rotation in reverse order
  A <- inv_rotation_matrix_lat %*% rbind(X, Y, Z)
  A <- inv_rotation_matrix_lon %*% A
  X_inv <- A[1, ]
  Y_inv <- A[2, ]
  Z_inv <- A[3, ]
  
  # Convert (X, Y, Z) back to spherical coordinates (longitude, latitude)
  Theta <- atan2(Y_inv, X_inv)
  Phi <- asin(Z_inv)
  
  Lon <- radians2degrees(Theta)
  Lat <- radians2degrees(Phi)
  
  return(list(Lon = Lon, Lat = Lat))
}