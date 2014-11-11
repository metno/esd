# Function that reads data stored on an irregular grid. The data is returned as a 'station' object.
retrieve.rcm <- function(ncfile,param=NULL,is=NULL,it=NULL) {
  ncold <- nc_open(ncfile)
  # Extract the time information: unit and time origin
  tunit <- ncatt_get( ncold, varid='time', attname='units')
  a <- regexpr("since",tunit)
  torg <- substr(tunit,a + attr(a,'match.length')+1,a + attr(a,'match.length')+10)
  tunit <- tolower(substr(tunit,1,a-2))
  xunit <- ncatt_get( ncold, varid=param, attname='units')
  longname <- ncatt_get( ncold, varid=param, attname='long_name')
  lat <- c(ncvar_get(ncold,varid='lat'))
  lon <- c(ncvar_get(ncold,varid='lon'))
  if (!is.null(is)) {
    nms <- names(is)
    ix <- grep("lon", tolower(substr(nms, 1, 3)))
    iy <- grep("lat", tolower(substr(nms, 1, 3)))
  } else {ix <- NA; iy <- NA}
  if (!is.null(it)) {
    if (is.character(it))
    
  } else it <- NA
  
  time <- ncvar_get(ncold,varid='time')
  time <- switch(str(tunit,1,4),'day'=as.Date(time+julian(as.Date(torg))))
  
  rcm <- ncvar_get(ncold,varid=param,start=start, count=count)
  nc_close( ncold )

  d <- dim(rcm)
  dim(rcm) <- c(d[1]*d[2],d[3])
  RCM <- zoo(t(rcm)*scal,order.by=time)
  attr(RCM,'longitude') <- lon
  attr(RCM,'latitude') <- lat
  attr(RCM,'altitude) <- rep(NA,length(lon))
  attr(RCM,'variable') <- param
  attr(RCM,'unit') <- xunit
  attr(RCM,'source') <- fname
  attr(RCM,'location') <- rep(NA,length(lon))
  attr(RCM,'longname') <- longname
  attr(RCM,'history') <- history.stamp(x)
  rm('rcm')
  class(RCM) <- c('station','day','zoo')
  return(RCM)
}
