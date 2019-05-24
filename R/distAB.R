#' Calculate distance between points on earth
#'
#' @param lon a longitude
#' @param lat a latitude
#' @param lons longitude or vector of longitudes 
#' @param lats latitude or vector of longitudes
#'
#' @return distance between [lon, lat] and [lons, lats] (unit: m)
#'
#' @export
distAB <- function(lon,lat,lons,lats,a=6.378e06) {
  good <- is.finite(lons) & is.finite(lats)
  lons <- lons[good]
  lats <- lats[good]
  if ( (length(lon) !=1) | (length(lat) !=1) |
       (length(lons)!=length(lats)) ) {
    print(paste("distAB [clim.pact]: length(lon)=",length(lon),
                "length(lat)=",length(lat),
                "length(lons)=",length(lons),"length(lats)=",length(lats)))
    print("length(lons) must equal length(lats) and lon and lat must have length=1")
    stop('Error in distAB - argument lengths do not match!')
  }
  theta <- pi*lon/180
  phi <- pi*lat/180
  dist <- rep(NA,length(lons))
  r1 <- c(cos(phi)*cos(theta),
          sin(phi),
          cos(phi)*sin(theta))
  dim(r1) <- c(3,length(lon))
  theta <- pi*lons/180
  phi <- pi*lats/180

  r2 <- cbind(cos(phi)*cos(theta),
              sin(phi),
              cos(phi)*sin(theta))
#  angle <- acos( sum(r1*r2) )
  angle <- acos( r2 %*% r1 )
#  if (sum(!is.finite(angle))>0) {
#    ibad <- !is.finite(angle)
#    print("MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM: distAB")
#    print("Detected a not-a-finite number...")
#    print(r2[ibad,])
#    print(r1)
#    print("MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM: distAB")
#  }
  dist <- rep(NA,length(lons))
  dist[good] <- a* angle
  dist
}
                   
