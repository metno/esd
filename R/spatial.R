# How is spatial.avg.field different from aggregate.area?
# I replaced spatial.avg.field with aggregate.area in DSensemble and iid.test
# and will not export it.

# Spatial Average of a Field Object.
# 
# Computes the spatial average of a field object and return a zoo time series
# object.
# 
# @aliases spatial
#
# @param x A zoo field object with two (longitude, latitude) or three
# dimensions (longitude, latitude, time)
#
# @return A "zoo" "time series" object
# 
# @seealso \code{\link{retrieve.ncdf4}}
#
# @examples
# 
# slp <- slp.DNMI(lon=c(-60,60),lat=c(30,65))
# slp.avg <- spatial.avg.field(slp)
# 
# @export spatial.avg.field
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
