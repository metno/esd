## Name		: spatial.avg.field
## Description	: Computes the spatial average of a field.
## Author 	: Abdelkader Mezghani, METNO
## contact 	: abdelkaderm@met.no
## Last Update	: 21-03-2013
## require	: zoo 
## input	: a zoo field object / 3 dimensional field with dimensions (time,lon,lat)

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
  ## return(z)
  invisible(z)
}
