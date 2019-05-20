## Name		: spatial.avg.field
## Description	: Computes the spatial average of a field.
## Author 	: Abdelkader Mezghani, METNO
## contact 	: abdelkaderm@met.no
## Last Update	: 21-03-2013
## require	: zoo 
## input	: a zoo field object / 3 dimensional field with dimensions (time,lon,lat)



#' Spatial Average of a Field Object.
#' 
#' Computes the spatial average of a field object and return a zoo time series
#' object.
#' 
#' 
#' @aliases spatial
#' @param x A zoo field object with two (longitude, latitude) or three
#' dimensions (longitude, latitude, time)
#' @return A "zoo" "time series" object
#' @author A. Mezghani, MET Norway
#' @seealso \code{\link{retrieve.ncdf4}}
#' @keywords parameter , element
#' @examples
#' 
#' \dontrun{
#' # Consider the "gcm" object from the e.g. in \link{retrieve.ncdf4}
#' # Compute the spatial average along lon and lat in gcm2
#' gcm2 <- spatial.avg.field(gcm)
#' # keep all attributes in gcm2
#' gcm2 <- attrcp(gcm,gcm2)
#' # Compute the annual mean 
#' gcm.am <- as.annual(gcm2,FUN='mean',na.rm=TRUE)
#' # keep all attributes in gcm.am
#' gcm.am <- attrcp(gcm2,gcm.am)
#' year <- index(gcm.am)
#' # Compute the anomalies relative to the period 1986-2005
#' agcm.ave <- gcm.am-mean(gcm.am[is.element(as.numeric(year),c(1986:2005))])
#' agcm.ave <- attrcp(gcm.am,agcm.ave)
#' # plot anomaly time series of the global mean temperature 
#' frame.metno(agcm.ave,col="black",cex.lab=0.75,cex.axis=0.75) ! Should be
#'   updated !
#' # Add vertical margin text with y-label
#' mtext("Anomaly values relative to 1986-2005",side=2,line=2,cex=0.75,las=3)
#' }
#' 
#' @export spatial.avg.field
spatial.avg.field <- function(x) {
  ## Get dimensions
  d <- attr(x,"dimensions")
  ## Get latitude values
  lat <- attr(x,"lat")

  ## Get index attribute and initilaze z for the outputs
  index <- index(x)
  z <- zoo(NA,order.by=index)

  ## copy/update attributes
  z <- attrcp(x,z,ignore=c("longitude","latitude","dimensions"))
  attr(z,"type") <- "area"
  attr(z,"station_id") <- "999"
  attr(z,"longitude") <- mean(lon(x),na.rm=TRUE)
  attr(z,"latitude") <- mean(lat(x),na.rm=TRUE)
  attr(z,"loc") <- paste('area mean of',src(x))
  ## update the class of the object
  ## class(x)[1] <- "station"
  ## class(z) <- class(x)

  ## transpose x to have time in last dimension
  x <- t(x)
 
  ## reorganize in 3D matrix
  dim(x) <- c(d[1],d[2],d[3])
  ## reorganize in 2D
  dim(x) <- c(d[1]*d[2],d[3])

  ## compute the weights for all lon and lat combination
  w <- cos(pi*lat/180) ; sw <- sum(w)
  wl <- rep(w,d[1])

  for (t in 1:d[3]) 
    z[t] <- sum(wl*x[,t],na.rm=TRUE)/sum(wl)
  z <- as.station(z)
  index(z) <- index
  ## return(z)
  invisible(z)
}
