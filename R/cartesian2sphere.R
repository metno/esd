#' Convert Cartesian coordinates to spherical coordinates 
#' 
#' @param lon A numeric vector of longitudes
#' @param lat A numeric vector of latitudes
#' @param lonR Center longitude of viewing angle
#' @param latR Center latitude of viewing angle
#' @param stereographic if TRUE return a stereographic projection, else skip the last projection step and return rotated Cartesian coordinates
#' @param verbose a boolean; if TRUE print information about progress
#' 
#' @returns A list(X=X, Y=Y, Z=Y) where X is the radial distance, Y is the polar angle, and Z is the azimuthal angle. 
#' When displaying data on a spherical grid, the longitude/latitude coordinate system is replaced with X/Z, i.e., 
#' plot(X, Y) instead of plot(lon, lat).
#' 
#' @export
cartesian2sphere <- function(lon, lat, lonR = 120, latR = 90, stereographic = TRUE,
                             verbose = FALSE) {
  
  if (verbose) print("cartesian2sphere: Projecting to a rotated sphere.")
  
  # --- 1. Convert to Cartesian (X, Y, Z) on a unit sphere ---
  # Standard spherical to Cartesian conversion.
  # Theta (lon) is azimuth, Phi (lat) is elevation.
  Theta <- degrees2radians(lon)
  Phi   <- degrees2radians(lat)
  
  # Note: Standard GIS (lon, lat) to (X,Y,Z) often is:
  # X = cos(Phi) * cos(Theta)
  # Y = cos(Phi) * sin(Theta)
  # Z = sin(Phi)
  # Your original function had X=sin(Theta)*cos(Phi) and Y=cos(Theta)*cos(Phi). 
  # We will stick to the standard convention for clarity, but you may need to 
  # adjust back if your 'rotM' depends on your original definition.
  # **We'll use the standard convention for a unit sphere:**
  X <- cos(Phi) * cos(Theta)
  Y <- cos(Phi) * sin(Theta)
  Z <- sin(Phi)
  
  # Store original dimensions
  d <- dim(X)
  
  # Flatten coordinates for matrix operations
  coords <- rbind(c(X), c(Y), c(Z))
  
  # --- 2. Define and Apply Rotation ---
  # The goal is to rotate the coordinate system such that the chosen
  # center of projection (lonR, latR) is moved to the North Pole (0, 0, 1)
  
  # Rotation 1: Rotate about Z-axis by -(lonR) to bring lonR to the X-Z plane (lon=0)
  # lonR is the longitude that will run along the prime meridian in the new system.
  # Rotation matrices often use radians, so ensure rotM handles this.
  lonR_rad <- degrees2radians(-lonR) # Rotate opposite to lonR
  M1 <- rotZ(lonR_rad)
  
  # Rotation 2: Rotate about Y-axis by -(90 - latR) to bring latR to the Z-axis (pole)
  # The angle is the co-latitude. Rotate opposite to co-latitude.
  # Rotate along the X-Z plane (around Y-axis).
  latR_co_rad <- degrees2radians(-(90 - latR)) 
  M2 <- rotY(latR_co_rad)
  
  # Combined rotation: M2 %*% M1
  A <- M2 %*% M1 %*% coords
  
  # New (rotated) coordinates:
  X_prime <- A[1, ]
  Y_prime <- A[2, ]
  Z_prime <- A[3, ]
  
  # --- 3. Apply Stereographic Projection (if requested) ---
  if (stereographic) {
    # The new North Pole is now at (0, 0, 1) in the X_prime, Y_prime, Z_prime system.
    # We project from the South Pole (0, 0, -1) onto the plane Z_prime=0.
    
    # Formula for projection from South Pole (0, 0, -1) onto Z'=0:
    # x = X' / (1 + Z')
    # y = Y' / (1 + Z')
    # The 'visible' mask is for points on the near side of the projection point.
    # A point is projected if it's NOT the projection point itself (Z' != -1).
    divisor <- 1 + Z_prime
    # Mask out points at or very near the projection pole (South Pole in this setup)
    mask_visible <- abs(divisor) > 1e-10 
    
    # Initialize projected coordinates
    x <- X_prime * NA
    y <- Y_prime * NA
    
    # Apply projection only to visible points
    x[mask_visible] <- X_prime[mask_visible] / divisor[mask_visible]
    y[mask_visible] <- Y_prime[mask_visible] / divisor[mask_visible]
    
    # Restore original dimensions
    if (!is.null(d)) {
      dim(x) <- d
      dim(y) <- d
      dim(mask_visible) <- d
    }
    
    output <- list(X = x, Y = y, visible = mask_visible)
    
  } else {
    # --- 4. Return Rotated Cartesian Coordinates ---
    # Restore original dimensions
    if (!is.null(d)) {
      dim(X_prime) <- d
      dim(Y_prime) <- d
      dim(Z_prime) <- d
    }
    # For a general rotation, a simple near-side/far-side mask is less meaningful.
    # The Z_prime coordinate is now the distance from the projection plane.
    output <- list(X = X_prime, Y = Y_prime, Z = Z_prime)
  }
  
  return(output)
}
# cartesian2sphere <- function(lon, lat, lonR=NULL, latR=NULL, verbose = FALSE,
#                              mask=TRUE, stereographic=TRUE) {
#   if(verbose) print("cartesian2sphere")
#   if(is.null(lonR)) lonR <- mean(lon, na.rm=TRUE)
#   if(is.null(latR)) latR <- mean(lat, na.rm=TRUE)
#   
#   Theta <- degrees2radians(lon)
#   Phi <- degrees2radians(lat)
#   
#   # Transform -> (X,Y,Z):
#   X <- sin(Theta)*cos(Phi)
#   Y <- cos(Theta)*cos(Phi)
#   Z <- sin(Phi)
#   
#   # Grid coordinates:
#   d <- dim(X)
#   
#   # Rotate data grid:
#   A <- rotM(x=0, y=0, z=lonR) %*% rbind(c(X), c(Y), c(Z))
#   A <- rotM(x=latR, y=0, z=0) %*% A
#   X <- A[1,]; Y <- A[2,]; Z <- A[3,]
#   
#   if(stereographic) {
#     # Apply stereographic projection for northern hemisphere
#     visible <- Y > 1e-16  # Logical mask for points on the near side
#     if(latR>0) {
#       x <- X / (1 - Z)
#       y <- Y / (1 - Z)
#     } else if(latR<0) {
#       x <- X / (1 + Z)
#       y <- Y / (1 + Z)
#     }
#     if(!is.null(d)) {dim(x) <- d; dim(y) <- d; dim(visible) <- d} #; dim(Z) <- d}
#     output <- list(X=x, Y=y, visible=visible)
#   } else{
#     # Or return rotated Cartesian coordinates X, Y and Z
#     if(!is.null(d)) {dim(X) <- d; dim(Y) <- d; dim(Z) <- d; dim(visible) <- d}
#     visible <- Y > 1e-16  # Logical mask for points on the near side
#     output <- list(X = X, Y = Y, Z = Z, visible=visible)
#   }  
# 
#   return(output)
#   #return(list(X=X, Y=Y, Z=Z))
# }
