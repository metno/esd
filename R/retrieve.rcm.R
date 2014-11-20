# Function that reads data stored on an irregular grid. The data is returned as a 'station' object.
retrieve.rcm <- function(ncfile,param=NULL,is=NULL,it=NULL,verbose=FALSE) {
  if (verbose) print(paste('retrieve.rcm',ncfile))
  ncold <- nc_open(ncfile)
  
  # Extract the time information: unit and time origin
  tatt <- ncatt_get( ncold, varid='time' )
  #if (verbose) print(names(tatt))
  itunit <- (1:length(names(tatt)))[is.element(substr(names(tatt),1,4),'unit')]
  tunit <- tatt[[itunit]]
  a <- regexpr("since",tunit)
  torg <- substr(tunit,a + attr(a,'match.length')+1,a + attr(a,'match.length')+10)
  tunit <- tolower(substr(tunit,1,a-2))
  #browser()
  
  # Extract unit etc for the parameter
  vatt <- ncatt_get( ncold, varid=param )
  if (verbose) print(names(vatt))
  ivunit <- (1:length(names(vatt)))[is.element(substr(names(vatt),1,4),'unit')]
  vunit <- vatt[[ivunit]]
  if (verbose) print(paste('unit: ',vunit,'; time unit: ',tunit,'; time origin: ',torg,sep=''))
  longname <- ncatt_get( ncold, varid=param, attname='long_name')
  if (is.null(longname)) {
    longname <- switch(param,'t2m'='temperature','tmax'='maximum temperature','tmin'='minimum temperature',
                             'precip'='precipitation','slp'='mean sea level pressure','pt'='precipitation',
                             'pr'='precipitation')
  }
  # Extract the spatial coordinates:
  vnames <- names(ncold$var)
  latid <- vnames[is.element(tolower(substr(vnames,1,3)),'lat')]
  lonid <- vnames[is.element(tolower(substr(vnames,1,3)),'lon')]
  lat <- ncvar_get(ncold,varid=latid)
  lon <- ncvar_get(ncold,varid=lonid)
  time <- ncvar_get(ncold,varid='time')
  d <- c(dim(lat),length(time))
  #str(lat); str(lon)
  if (verbose) print(paste('region: ',min(lon),'-',max(lon),'E /',min(lat),'-',max(lat)))
  
  # Extract only the region of interest: only read the needed data
  #browser()
  if (!is.null(is)) {
    if (inherits(is,c('field','station'))) {
      y <- is
      if (verbose) print(paste('Use spatial coverage from an object:',floor(min(c(lon(y)))),'-',
                          ceiling(max(c(lon(y)))),'E /',floor(min(c(lat(y)))),'-',ceiling(max(c(lat(y)))),'N'))
      if (!is.null(attr(y,'lon_ref')) & !is.null(attr(y,'lat_ref'))) 
        is <- list( lon=attr(y,'lon_ref'),
                    lat=attr(y,'lat_ref') ) else
        is <- list(lon=c(floor(min(c(lon(y)))),ceiling(max(c(lon(y))))),
                   lat=c(floor(min(c(lat(y)))),ceiling(max(c(lat(y))))))
      rm('y')
    } 
    if (is.list(is)) {
      nms <- names(is)    
      iy <- grep("lat", tolower(substr(nms, 1, 3)))
      if (length(iy)>0) {
        lat.rng <- range(is[[iy]])
        mx <- trunc(d[1]/2)  # use the middle of the region for defining latitude range
        suby <- (lat.rng[1] <= lat[mx,]) & (lat.rng[2] >= lat[mx,])
        if (sum(suby)==0) stop(paste('retrieve.rcm: problems, the requested latitude range (',
                                      lat.rng[1],'-',lat.rng[2],') is not within present data (',
                                      min(lat[mx,]),'-',max(lat[mx,]),')'))
        #print(lat[mx,suby])
        starty <- min( (1:length(lat[1,]))[suby] )
        county <- sum(suby)
        if (verbose) print(paste('latitudes:',min(is[[iy]]),'-',max(is[[iy]]),
                                  'extracted:',min(lat[,suby]),'-',max(lat[,suby]),
                                  'start=',starty,'count=',county))
      } else if (length(grep("iy", tolower(substr(nms, 1, 2))))>0) {
      # Select single columns or rows in the spatial maxtix:
        iy <- grep("iy", tolower(substr(nms, 1, 2)))
        starty <- min(is[[iy]]); county <- max(is[[iy]])
      } else {starty <- 1; county <- d[2]}
  
      ix <- grep("lon", tolower(substr(nms, 1, 3)))
      if (length(ix)>0) {
      # The coordinates lon and lat are [X,Y] maxtrices:
        my <- trunc(d[2]/2) # use the middle of the region for defining longitude range
        lon.rng <- range(is[[ix]]) 
        subx <- (lon.rng[1] <= lon[,my]) & (lon.rng[2] >= lon[,my])
        if (sum(subx)==0) stop(paste('retrieve.rcm: problems, the requested longitude range (',
                                      lon.rng[1],'-',lon.rng[2],') is not within present data (',
                                      min(lon[,my]),'-',max(lon[,my]),')'))
        startx <- min( (1:length(lon[,my]))[subx] )
        countx <- sum(subx)
        if (verbose) print(paste('longitudes:',min(is[[ix]]),'-',max(is[[ix]]),
                                  'extracted:',min(lon[subx,]),'-',max(lon[subx,]),
                                  'start=',startx,'count=',countx))
      } else if (length(grep("ix", tolower(substr(nms, 1, 2))))>0) {
      # Select single columns or rows in the spatial maxtix:
        ix <- grep("ix", tolower(substr(nms, 1, 2)))
        startx <- min(is[[ix]]); countx <- max(is[[ix]])
      } else {startx <- 1; countx <- d[1]}
    
    } else if (is.numeric(is) | is.integer(is)) {
    # Select 
      if (verbose) print('Select by spatial index')
      startx <- 1
      starty <- min(it) %/% d[1] + 1
      countx <- d[1]
      county <- max(it) %/% d[1] + 1
      if (verbose) print(paste('selecting is: ',min(is),'-',max(is), 'reads rows',starty,'to',county))
    }
  } else {
    startx <- 1; countx <- d[1]; 
    starty <- 1; county <- d[2];
    subx <- rep(TRUE,d[1]); suby <- rep(TRUE,d[2])
  }
  
  # Extract only the time of interest: assume only an interval
  #print(tunit); browser()
  time <- switch(substr(tunit,1,3),'day'=as.Date(time+julian(as.Date(torg))),
       'mon'=as.Date(julian(as.Date(paste(time%/%12,time%%12+1,'01',sep='-'))) + julian(as.Date(torg))))
  if (verbose) print(paste(start(time),end(time),sep=' - '))
  if (!is.null(it)) {
    if (inherits(it,c('field','station'))) {
      if (verbose) print('Use time interval from an object')
      y <- it
      it <- c(start(y),end(y))
      rm('y')
    }
    if (inherits(it,'Date')) {
      startt <- min( (1:length(time))[it >= time] )
      stoptt <- max( (1:length(time))[it <= time] )
      countt <- max(it) - startt + 1
    } else if (sum(is.element(it,1600:2200)) > 0) {
    # Assume the years:
      startt <- min( (1:length(time))[it >= year(time)] )
      stoptt <- max( (1:length(time))[it <= year(time)] )
      countt <- max(it) - startt + 1
    } else if ( (max(it) <= length(time)) & min(it >= 1) ) {
      startt <- min(it); countt <- max(it) - startt + 1
    } 
  } else {startt <- 1; countt <- length(time); it <- NA}
  
  # This information is used when retrieve.rcm is used again to extract similar region
  mx <- trunc(d[1]/2); my <- trunc(d[2]/2)
  lon.ref <- range(lon[subx,my])
  lat.ref <- range(lat[mx,suby])
  
  # Test the dimensions so that the count does not exceed the array:
  if (startx + countx - 1 > d[1]) {
    countx <- d[1] - startx + 1
    warning("retrieve.rcm: number of points along the longitude exceeds data dimensions")
  }
  if (starty + county - 1 > d[2]) {
    county <- d[2] - starty + 1
    warning("retrieve.rcm: number of points along the latitude exceeds data dimensions")
  }
  if (startt + countt - 1> d[3]) {
    countt <- d[3] - startt + 1
    warning("retrieve.rcm: number of points in time exceeds data dimensions")
  }
  
  lon <- lon[startx:(startx+countx-1),starty:(starty+county-1)]
  lat <- lat[startx:(startx+countx-1),starty:(starty+county-1)]
  time <- time[startt:(startt+countt-1)]
  start <- c(startx,starty,startt)
  count <- c(countx,county,countt)
  
  if (verbose) {print(start); print(count)}
  rcm <- ncvar_get(ncold,varid=param,start=start, count=count)
  nc_close( ncold )
  
  d <- dim(rcm)
  dim(rcm) <- c(d[1]*d[2],d[3])
  
  if (is.numeric(is) | is.integer(is)) {
  # If only reading a set index, then remove the ones before and after, i.e. read the first 1000 grid points or
  # the next 1000 grid points. Useful for processing the data chunck-wise.
    i1 <- min(it) %% d[1]
    i2 <- (max(it) %/% d[1])*d[1] + max(it) %% d[1]
    rcm <- rcm[i1:i2,]
    lon <- lon[i1:i2]; lat <- lat[i1:i2]
  }
  RCM <- zoo(t(rcm),order.by=time)
  attr(RCM,'longitude') <- c(lon)
  attr(RCM,'latitude') <- c(lat)
  attr(RCM,'lat_ref') <- lat.ref
  attr(RCM,'lon_ref') <- lon.ref
  attr(RCM,'altitude') <- rep(NA,length(lon))
  attr(RCM,'variable') <- param
  attr(RCM,'unit') <- vunit
  attr(RCM,'source') <- ncfile
  attr(RCM,'location') <- rep(NA,length(lon))
  attr(RCM,'longname') <- longname
  attr(RCM,'history') <- history.stamp()
  rm('rcm')
  class(RCM) <- c('station','day','zoo')
  return(RCM)
}
