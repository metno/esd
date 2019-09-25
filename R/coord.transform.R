## Original Javascript by Chuck Taylor
## Port to C++ by Alex Hajnal (2)
## Port from C++ to R by Ketil Tunheim
##
## This script
##
## 1) http://home.hiwaay.net/~taylorc/toolbox/geography/geoutm.html
## 2) http://alephnull.net/software/gis/UTM_WGS84_C_plus_plus.shtml

# globals
sm_a <- 6378137.0
sm_b <- 6356752.314
UTMScaleFactor <- 0.9996

## Inputs are lat and lon vectors in degrees, and UTM zone to convert to
## Output is a vector with (easting, northing)
LatLon2UTM <- function(lat, lon, zone) {
  if ( (zone < 1) || (zone > 60) ) {
    zone <- floor((lon + 180.0) / 6) + 1
  }
  cmeridian = UTMCentralMeridian(zone)
  
  if (length(lat) != length(lon))
    stop("Lat and lon lists must be the same length")
  
  X <- c()
  Y <- c()
  for (i in 1:length(lat)) {
    if (is.null(lat[[i]]) || is.null(lon[[i]])) {
      X[[i]] <- NA
      Y[[i]] <- NA
    } else {
      XY <- MapLatLon2XY(deg2rad(lat[[i]]), deg2rad(lon[[i]]), cmeridian)
      XY[1] <- XY[1] * UTMScaleFactor + 500000
      XY[2] <- XY[2] * UTMScaleFactor
      if (XY[2] < 0) XY[2] <- XY[2] + 10000000
      
      X[[i]] <- round(XY[1], 3)
      Y[[i]] <- round(XY[2], 3)
    }
  }
  list(X, Y)
}

## Inputs are x (easting), y (northing), UTM zone they are defined in,
## and a boolean to distinguish the southern hemisphere
## Output is a vector with (latitude, longitude)
UTM2LatLon <- function(x, y, zone, southhemi=FALSE) {
  x <- (x - 500000) / UTMScaleFactor
  
  if (southhemi) y <- y - 10000000
  y <- y/UTMScaleFactor
  
  cmeridian = UTMCentralMeridian(zone)
  
  if (length(x) != length(x))
    stop("X and Y lists must be the same length")
  
  Lat <- c()
  Lon <- c()
  for (i in 1:length(x)) {
    if (is.null(x[[i]]) || is.null(y[[i]])) {
      Lat[[i]] <- NA
      Lon[[i]] <- NA
    } else {
      LatLon <- MapXY2LatLon(x[i], y[i], UTMCentralMeridian(zone))
    
      Lat[[i]] <- round(rad2deg(LatLon[1]), 4)
      Lon[[i]] <- round(rad2deg(LatLon[2]), 4)
    }
  }
  list(Lat, Lon)
}

## The math used by LatLon2UTM
MapLatLon2XY <- function(phi, lambda, lambda0) {
  ep2 <- (sm_a**2 - sm_b**2)/sm_b**2
  nu2 <- ep2 * cos(phi)**2
  N <- (sm_a**2) / (sm_b*sqrt(1+nu2))
  t <- tan(phi)
  l <- lambda - lambda0
  
  t2 <- t*t
  t4 <- t2*t2
  t6 <- t4*t2
  
  l3coef <- 1 - t2 + nu2
  l4coef <- 5 - t2 + 9*nu2 + 4*(nu2**2)
  l5coef <- 5 - 18*t2 + t4 + 14*nu2 - 58*t2*nu2
  l6coef <- 61 - 58*t2 + t4 + 270*nu2 - 330*t2*nu2
  l7coef <- 61 - 479*t2 + 179*t4 - t6
  l8coef <- 1385 - 3111*t2 + 543*t4 - t6
  
  x <- N*cos(phi)*l +
    N/6 * cos(phi)**3 * l3coef * (l**3) +
    N/120 * cos(phi)**5 * l5coef * (l**5) +
    N/5040 * cos(phi)**7 * l7coef * (l**7)
  
  y <- ArclengthOfMeridian(phi) +
    t/2*N * cos(phi)**2 * (l**2) +
    t/24*N * cos(phi)**4 * l4coef * (l**4) +
    t/720*N * cos(phi)**6 * l6coef * (l**6) +
    t/40320*N * cos(phi)**8 * l8coef * (l**8)
  
  c(x,y)
}

## The math used by UTM2LatLon
MapXY2LatLon <- function(x, y, lambda0) {
  phif <- FootpointLatitude(y)
  
  ep2 <- (sm_a**2 - sm_b**2)/sm_b**2
  cf <- cos(phif)
  nuf2 <- ep2 * (cf**2)
  Nf <- (sm_a**2) / (sm_b*sqrt(1+nuf2))
  tf <- tan(phif)
  tf2 <- tf*tf
  tf4 <- tf2*tf2
  
  Nfpow <- Nf
  x1frac <- 1/(Nfpow*cf)
  Nfpow <- Nfpow * Nf
  x2frac <- tf/(2*Nfpow)
  Nfpow <- Nfpow * Nf
  x3frac <- 1/(6*Nfpow*cf)
  Nfpow <- Nfpow * Nf
  x4frac <- tf/(24*Nfpow)
  Nfpow <- Nfpow * Nf
  x5frac <- 1/(120*Nfpow*cf)
  Nfpow <- Nfpow * Nf
  x6frac <- tf/(720*Nfpow)
  Nfpow <- Nfpow * Nf
  x7frac <- 1/(5040*Nfpow*cf)
  Nfpow <- Nfpow * Nf
  x8frac <- tf/(40320*Nfpow)
  
  x2poly <- -1 - nuf2
  x3poly <- -1 -2*tf2 - nuf2
  x4poly <- 5 + 3*tf2 + 6*nuf2 - 6*tf2*nuf2
  x5poly <- 5 + 28*tf2 + 24*tf4 + 6*nuf2 + 8*tf2*nuf2
  x6poly <- -61 - 90*tf2 - 45*tf4 - 107*nuf2 + 162*tf2*nuf2
  x7poly <- -61 - 662*tf2 - 1320*tf4 - 720*(tf4*tf2)
  x8poly <- 1385 + 3633*tf2 + 4095*tf4 + 1575*(tf4*tf2)   
  
  # latitude
  phi <- phif + x2frac*x2poly*(x*x) +
    x4frac*x4poly*(x**4) +
    x6frac*x6poly*(x**6) +
    x8frac*x8poly*(x**8)
  
  # longitude
  lambda <- lambda0 + x1frac*x +
    x3frac*x3poly*(x**3) +
    x5frac*x5poly*(x**5) +
    x7frac*x7poly*(x**7)
  
  c(phi,lambda)
}

## Helping function for MapLatLon2XY
ArclengthOfMeridian <- function(phi) {
  n <- (sm_a - sm_b)/(sm_a + sm_b)
  alpha <- (sm_a + sm_b)/2 * (1 + (n**2)/4 + (n**4)/64)
  beta  <- -3*n/2 + 9*(n**3)/16 - 3*(n**5)/32
  gamma <- 15*(n**2)/16 - 15*(n**4)/32
  delta <- -35*(n**3)/48 + 105*(n**5)/256
  epsilon <- 315*(n**4)/512
  
  result <- alpha * (phi +
                       beta*sin(2*phi) +
                       gamma*sin(4*phi) +
                       delta*sin(6*phi) +
                       epsilon*sin(8*phi)
  )
  
  result
}

## Helping function for MapXY2LatLon
FootpointLatitude <- function(y) {
  n <- (sm_a - sm_b)/(sm_a + sm_b)
  alpha_ <- (sm_a + sm_b)/2 * (1 + (n**2)/4 + (n**4)/64)
  beta_  <- 3*n/2 - 27*(n**3)/32 + 269*(n**5)/512
  gamma_ <- 21*(n**2)/16 - 55*(n**4)/32
  delta_ <- 151*(n**3)/96 - 417*(n**5)/128
  epsilon_ <- 1097*(n**4)/512
  
  y_ <- y / alpha_
  result <- y_ + 
    beta_*sin(2*y_) +
    gamma_*sin(4*y_) +
    delta_*sin(6*y_) +
    epsilon_*sin(8*y_)
  
  result
}

## Conversion between degrees and radians
deg2rad <- function(deg) {
  as.numeric(deg)/180*pi
}
rad2deg <- function(rad) {
  as.numeric(rad)/pi*180
}
## Convert from UTM zone to central meridian
UTMCentralMeridian <- function(zone) {
  deg2rad(-183 + zone*6)
}