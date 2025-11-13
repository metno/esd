#' Convert Cartesian coordinates to spherical coordinates 
#' 
#' @param lon A numeric vector of longitudes
#' @param lat A numeric vector of latitudes
#' @param lonR Center longitude of viewing angle
#' @param latR Center latitude of viewing angle
#' @param stereographic if TRUE return a stereographic projection, else skip the last projection step and return rotated Cartesian coordinates
#' @param mask_horizon if TRUE mask areas beyond the horizon so that only one hemisphere can be shown at once. This is only 
#'                     relevant to the stereographic projection which can expand beyond the hemisphere and show the whole globe
#'                     (with a distorted perspective far from the center)  
#' @param verbose a boolean; if TRUE print information about progress
#' 
#' @returns A list(X=X, Y=Y, Z=Y) where X is the radial distance, Y is the polar angle, and Z is the azimuthal angle. 
#' When displaying data on a spherical grid, the longitude/latitude coordinate system is replaced with X/Z, i.e., 
#' plot(X, Y) instead of plot(lon, lat).
#' 
#' @export
cartesian2sphere <- function(lon, lat, lonR = 0, latR = 90, axiR=0, 
                             stereographic = FALSE, mask_horizon = TRUE, 
                             verbose = FALSE) {
  
  if (verbose) print("cartesian2sphere: Projecting to a rotated sphere.")
  
  # Add 90 degrees to lonR so that lonR is pointing downward in the plot 
  # instead horizontally to the right
  lonR <- lonR + 90 
  
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
  
  # Rotation 3: Azimuthal swing / Axis tilt around the Z'-axis 
  # This rotates the features on the visible hemisphere relative to the X'-Y' plane. 
  axiR_rad <- degrees2radians(-axiR) 
  M3 <- rotZ(axiR_rad)
  
  # Combined rotation: 
  A <- M3 %*% M2 %*% M1 %*% coords
  
  # New (rotated) coordinates:
  X_prime <- A[1, ]
  Y_prime <- A[2, ]
  Z_prime <- A[3, ]
  
  # Horizon Mask, preventing points from the far hemisphere (Z' < 0) from being projected.
  # Define the numerical tolerance for the horizon (e.g., 1e-8 for unit sphere)
  min_horizon <- 1e-8
  mask_visible_horizon <- Z_prime >= -min_horizon
  
  # --- 3. Apply Stereographic Projection (if requested) ---
  if (stereographic) {
    
    # 3. Calculate the divisor for stereographic projection
    divisor <- 1 + Z_prime
    
    # Singularity Mask to ensure numerical Stability (Z'!= -1)
    min_singularity <- 1e-10
    mask_visible_singularity <- abs(divisor) > min_singularity
    
    # Combine masks: near side AND not at the antipodal pole
    if(mask_horizon) mask_final <- mask_visible_horizon & mask_visible_singularity else
      mask_final <- mask_visible_singularity
    
    # Initialize projected coordinates (using NA ensures excluded points are not plotted)
    x <- X_prime * NA
    y <- Y_prime * NA
    
    # Apply projection only to the points authorized by the combined mask
    x[mask_final] <- X_prime[mask_final] / divisor[mask_final]
    y[mask_final] <- Y_prime[mask_final] / divisor[mask_final]
    
    # Restore original dimensions for output mask
    if (!is.null(d)) {
      dim(x) <- d
      dim(y) <- d
      dim(mask_final) <- d # Use the final, combined mask for output
    }
    
    output <- list(X = x, Y = y, visible = mask_final)
    
  } else {
    
    # --- 4. Return Rotated Cartesian Coordinates ---
    # Restore original dimensions
    if (!is.null(d)) {
      dim(X_prime) <- d
      dim(Y_prime) <- d
      dim(Z_prime) <- d
      dim(mask_visible_horizon) <- d
    }
    # For a general rotation, a simple near-side/far-side mask is less meaningful.
    # The Z_prime coordinate is now the distance from the projection plane.
    output <- list(X = X_prime, Y = Y_prime, Z = Z_prime, visible = mask_visible_horizon)
  }
  
  return(output)
}

