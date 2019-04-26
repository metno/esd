# Function that reads data stored on an irregular grid. The data is returned as a 'station' object.
retrieve.rcm <- function(ncfile,path=NULL,param=NULL,is=NULL,it=NULL,verbose=FALSE) {
    if(verbose) print("retrieve.rcm")
    if (!is.null(path))
        ncfile <- file.path(path,ncfile,fsep = .Platform$file.sep)
    if (verbose) print(paste('retrieve ',ncfile))
    ncold <- nc_open(ncfile)

    # Extract the time information: unit and time origin
    tatt <- ncatt_get( ncold, varid='time' )
                                        #if (verbose) print(names(tatt))
    itunit <- (1:length(names(tatt)))[is.element(substr(names(tatt),1,4),'unit')]
    tunit <- tatt[[itunit]]
    tcal <-""
    if (sum(is.element(substr(names(tatt),1,4),'cale'))!=0) {
      if (verbose) print("Calender found")
      tcal <- tatt$cale
    }
    a <- regexpr("since",tunit)
    torg <- substr(tunit,a + attr(a,'match.length')+1,a + attr(a,'match.length')+10)
    torig <- paste(unlist(strsplit(tunit," "))[3:4],collapse=" ")
    tunit <- tolower(substr(tunit,1,a-2))

    # Extract unit etc for the parameter
    vatt <- ncatt_get( ncold, varid=param )
    # 16.11.2017 hbe added option names(vatt below)
    if(any(grepl('unit',vatt)) | any(grepl('unit',names(vatt)))) {
      ivunit <- (1:length(names(vatt)))[is.element(substr(names(vatt),1,4),'unit')]
      vunit <- vatt[[ivunit]]
    } else {
      if(verbose) print("Unkown unit")
      vunit <- ""
    }
    if (verbose) print(paste('unit: ',vunit,'; time unit: ',tunit,'; time origin: ',torg,sep=''))
    #HBE added value 12.11.2017
    longname <- ncatt_get( ncold, varid=param, attname='long_name')$value
    if (is.null(longname)) {
        longname <- switch(param,'t2m'='temperature','tmax'='maximum temperature','tmin'='minimum temperature',
                           'precip'='precipitation','slp'='mean sea level pressure','pt'='precipitation',
                           'pr'='precipitation')
    }
    # Extract the spatial coordinates:
    vnames <- names(ncold$var)
    ##latid <- vnames[is.element(tolower(substr(vnames,1,3)),'lat')]
    ##lonid <- vnames[is.element(tolower(substr(vnames,1,3)),'lon')]
    #latid <- vnames[grep('lat',tolower(vnames))]
    #lonid <- vnames[grep('lon',tolower(vnames))]
    ## KMP 2016-12-20: grep('lat',...) sometimes finds more than 1 match
    latid <- vnames[tolower(vnames) %in% c("lat","latitude")]
    lonid <- vnames[tolower(vnames) %in% c("lon","longitude")]
    lat <- ncvar_get(ncold,varid=latid)
    lon <- ncvar_get(ncold,varid=lonid)
    time <- ncvar_get(ncold,varid='time')
    d <- c(dim(lat),length(time))
                                        #str(lat); str(lon)
    if (verbose) print(paste('region: ',round(min(lon),digits=2),'-',round(max(lon,digits=2)),'E /',round(min(lat),digits=2),'-',round(max(lat),digits=2),'N'))

                                        # Extract only the region of interest: only read the needed data

    if (!is.null(is)) {
        if (inherits(is,c('field','station'))) {
            y <- is
            if (verbose) print(paste('Use spatial coverage from an object:',floor(min(c(lon(y)))),'-',
                                     ceiling(max(c(lon(y)))),'E /',floor(min(c(lat(y)))),'-',ceiling(max(c(lat(y)))),'N'))
                                        #      if (!is.null(attr(y,'lon_ref')) & !is.null(attr(y,'lat_ref')))
                                        #        is <- list( lon=attr(y,'lon_ref'),
                                        #                    lat=attr(y,'lat_ref') ) else
            is <- list(lon=c(floor(min(c(lon(y)))),ceiling(max(c(lon(y))))),
                       lat=c(floor(min(c(lat(y)))),ceiling(max(c(lat(y))))))
            rm('y')
        }
        if (is.list(is)) {
            nms <- names(is)
            iy <- grep("lat", tolower(substr(nms, 1, 3)))
            if (length(iy)>0) {
                lat.rng <- range(is[[iy]])
                                        # use the smallest and largest
                if(length(dim(lat))==2) {
                  latn <- apply(lat,2,min); latx <- apply(lat,2,min)
                } else {
                  latn <- lat; latx <- lat
                }
                suby <- (lat.rng[1] <= latn) & (lat.rng[2] >= latx)
                starty <- min( (1:length(latx))[suby] )
                county <- sum(suby)
                if (sum(suby)==0) stop(paste('retrieve.rcm: problems, the requested latitude range (',
                           lat.rng[1],'-',lat.rng[2],') is not within present data (',
                           min(latn),'-',max(latx),')'))
                                        #print(lat[mx,suby])

                if (verbose) print(paste('latitudes:',min(is[[iy]]),'-',max(is[[iy]]),
                                         'extracted:',min(latn[suby]),'-',max(latn[suby]),
                                         'start=',starty,'count=',county))
            } else if (length(grep("iy", tolower(substr(nms, 1, 2))))>0) {
                                        # Select single columns or rows in the spatial maxtix:
                iy <- grep("iy", tolower(substr(nms, 1, 2)))
                starty <- min(is[[iy]]); county <- length(is[[iy]])
                suby <- is.element(1:d[2],iy)
            } else {starty <- 1; county <- d[2]; suby <- is.finite(1:d[2])}

            ix <- grep("lon", tolower(substr(nms, 1, 3)))
            if (length(ix)>0) {
              lon.rng <- range(is[[ix]])
              if(length(dim(lat))==2) {
                # The coordinates lon and lat are [X,Y] maxtrices:
                lonn <- apply(lon,1,min); lonx <- apply(lon,1,min)
              } else {
                lonn <- lon; lonx <- lon
              }
              subx <- (lon.rng[1] <= lonn) & (lon.rng[2] >= lonx)
              startx <- min( (1:length(lonx))[subx] )
              countx <- sum(subx)
              if (sum(subx)==0) stop(paste('retrieve.rcm: problems, the requested longitude range (',
                         lon.rng[1],'-',lon.rng[2],') is not within present data (',
                         min(lonn),'-',max(lonx),')'))

              if (verbose) print(paste('longitudes:',min(is[[ix]]),'-',max(is[[ix]]),
                                         'extracted:',min(lonn[subx]),'-',max(lonn[subx]),
                                         'start=',startx,'count=',countx))
            } else if (length(grep("ix", tolower(substr(nms, 1, 2))))>0) {
                                        # Select single columns or rows in the spatial maxtix:
              ix <- grep("ix", tolower(substr(nms, 1, 2)))
              startx <- min(is[[ix]]); countx <- length(is[[ix]])
              subx <- is.element(1:d[1],ix)
            } else {startx <- 1; countx <- d[1]; subx <- is.finite(1:d[1])}

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
    time <- switch(substr(tunit,1,3),
                   'day'=as.Date(time+julian(as.Date(torg))),
                   'mon'=as.Date(julian(as.Date(paste(time%/%12,time%%12+1,'01',sep='-'))) + julian(as.Date(torg))),
                   'hou'=strptime(torig,format="%Y-%m-%d %H") + time*3600,
                   'sec'=strptime(torig,format="%Y-%m-%d %H") + time)
    # next save for later if adding no_leap func
    #if ((tcal %in% c("365_day", "365day", "no_leap", "no leap")) && (any(grepl('hou',tunit))) && ((diff(ttest)>29) && (diff(ttest) <= 31 )) )
    #HBE 2018/1/17 saving POSIX with monthly freq as Date at month start
    if (((diff(time)>=28) && (diff(time) <= 31 )) | (length(time) == 1)) {
      time <- as.Date(strftime(time, format="%Y-%m-01"))
      if (verbose) print("monthly frequency, saving as Date Y-m-01")
      }
    if (verbose) print(paste(start(time),end(time),sep=' - '))
    if (!is.null(it)) {
        if (inherits(it,c('field','station'))) {
            if (verbose) print('Use time interval from an object')
            y <- it
            it <- c(start(y),end(y))
            rm('y')
        }

        if (is.character(it)) {
          if ((levels(factor(nchar(it)))==10)) {
            it <- as.Date(it,format="%Y-%m-%d")
          } else if ((levels(factor(nchar(it)))==8)) {
            it <- as.Date(it,format="%Y%m%d")
          } else if (levels(factor(nchar(it)))==4) {
            it <- as.Date(c(paste(it[1],'-01-01',sep=''),
                            paste(it[2],'-12-31',sep='')))
          }
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
        } else {
            print(paste("unkown format of input it:",it))
        }
    } else {startt <- 1; countt <- length(time); it <- NA}

                                        # This information is used when retrieve.rcm is used again to extract similar region
                                        #mx <- trunc(d[1]/2); my <- trunc(d[2]/2)
                                        #lon.ref <- range(lon[subx,my])
                                        #lat.ref <- range(lat[mx,suby])

                                        # Test the dimensions so that the count does not exceed the array:

      d1 <- d[1]
      if(length(d)==3) {
        d2 <- d[2]; d3 <- d[3]
      } else {
        d2 <- d[1]; d3 <- d[2]
      }
      if (startx > d1) {
          startx <- d1
          warning("retrieve.rcm: points along the longitude exceed data dimensions")
      }
      if (starty > d2) {
          starty <- d2
          warning("retrieve.rcm: points along the latitude exceed data dimensions")
      }
      if (startt > d3) {
          startt <- d3
          warning("retrieve.rcm: points in time exceed data dimensions")
      }
      if (startx < 1) {
          startx <- 1
          warning("retrieve.rcm: points along the longitude exceed data dimensions")
      }
      if (starty < 1) {
          starty <- 1
          warning("retrieve.rcm: points along the latitude exceed data dimensions")
      }
      if (startt < 1) {
          startt <- 1
          warning("retrieve.rcm: points in time exceed data dimensions")
      }
                                        # Test the dimensions so that the count does not exceed the array:
    if (startx + countx - 1 > d1) {
        countx <- d[1] - startx + 1
        warning("retrieve.rcm: number of points along the longitude exceeds data dimensions")
    }
    if (starty + county - 1 > d2) {
        county <- d[2] - starty + 1
        warning("retrieve.rcm: number of points along the latitude exceeds data dimensions")
    }
    if (startt + countt - 1> d3) {
        countt <- d[3] - startt + 1
        warning("retrieve.rcm: number of points in time exceeds data dimensions")
    }
      #HBE added y-dimension to lon as well as lat (before only lat had)
    if(length(d)==3) {
      lon <- lon[startx:(startx+countx-1),starty:(starty+county-1)]
      lat <- lat[startx:(startx+countx-1),starty:(starty+county-1)]
      time <- time[startt:(startt+countt-1)]
      start <- c(startx,starty,startt)
      count <- c(countx,county,countt)
    } else {
      startxy <- max(startx,starty)
      countxy <- min(countx,county)
      lon <- lon[startxy:(startxy+countxy-1)]
      lat <- lat[startxy:(startxy+countxy-1)]
      time <- time[startt:(startt+countt-1)]
      start <- c(startxy,startt)
      count <- c(countxy,countt)
    }
    if (verbose) {print(start); print(count)}
    rcm <- ncvar_get(ncold,varid=param,start=start, count=count)
    nc_close( ncold )
    if(length(dim(rcm))==3) {
      d <- dim(rcm)
    } else if (length(dim(rcm))!=3 & length(d)==3) {
      # If there are less than 3 dimensions, add one dummy dimension. To avoid crashes...
      n1 <- (1:3)[count==1]; nm <- (1:3)[count>1]
      D <- rep(1,3); D[nm] <- d; D[n1] <- 1
      d <- D; rm('D')
    } else {
      d <- c(1,dim(rcm))
    }
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
    attr(RCM,'count') <- count
    attr(RCM,'start') <- start
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
