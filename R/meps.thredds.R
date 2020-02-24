#' Function to read weather forecast data for the Nordic domain from thredds.met.no
#' 
#' The data are post processed forecasts or analysis based on the MetCoOp Ensemble Prediction System (MEPS).
#' As the resolution is high both in space (1x1km) and time (1 hours), it may be necessary
#' to limit the spatial and/or temporal scope (see input options 'it', 'lon', 'lat')
#' as you will otherwise run out of memory.
#' 
#' @aliases meps
#'
#' @param url URL for the data on thredds.met.no
#' @param lon Longitude selection (=NULL reads all)
#' @param lat Latitude selection (=NULL reads all)
#' @param param Variable name (rr = precipitation, tg = temperature)
#' @param FUN Function for daily aggregation. =NULL gives raw data
#' @param it Index time - a range of dates or years to select (e.g., it=as.Date(c("2010-01-01","2010-01-31")))
#' @param dt Number of time steps to access at a time (default: 50). 
#' A smaller value can be set when requesting a large spatial domain as the server 
#' doesn't like to deal with large data amounts.
#' @param verbose write out diagnostics
#' @param plot plot the results while reading. 
#' @param ncfile netCDF file to access (path to a file or url)
#' @param path path to ncfile (use only if ncfile doesn't contain the full path)
#' @return A "zoo" "station" object with additional attributes used for further
#' processing.
#'
#' @seealso senorge, station.thredds
#' 
#' @examples
#' it <- "latest"
#' lon <- c(-2,15)
#' lat <- c(55,63)
#' slp <- meps(param="slp", lon=lon, lat=lat, it=it, verbose=TRUE)
#' map(slp, FUN="mean")
#' 
#' @export
meps <- function(url='https://thredds.met.no/thredds/catalog/metpplatest',
                 type="forecast", param='slp', lon=c(9.5,11.5), lat=c(59,61), 
                 it="latest", dt=50, verbose=FALSE, plot=FALSE) {
  if (verbose) print('esd::meps')
  path <- sub("/catalog/","/dodsC/",url)
  contents <- readLines(paste0(url,'/catalog.html'))
  files <- contents[grepl(".nc",contents,fixed=TRUE)]
  files <- sort(gsub(".*<tt>|</tt>.*","",files,fixed=FALSE))
  files <- files[grepl(type,files)]
  if(is.character(it)) {
    files <- files[grepl(it,files)]
  } else if(is.numeric(it)) {
    files <- files[grepl(paste0(as.character(it),collapse="|"),files)]
  } else if(inherits(it,c("Date","POSIXt"))) {
    files <- files[!grepl("latest",dates)]
    dates <- gsub(".*_|.nc","",files)
    dates <- gsub("[A-Z]|[a-z]","",dates)
    if(inherits(it,"Date")) {
      dates <- as.Date(dates,format="%Y%m%d%H")
    } else {
      dates <- strptime(dates,format="%Y%m%d%H")
    }
    ok <- dates>=min(it) & dates<=max(it)
    files <- files[ok] 
  }
  if(length(files)>0) {
    ncid <- nc_open(file.path(path,files[[1]]))
    nc_close(ncid)
    varnames <- names(ncid$var)
    if(!param %in% varnames) {
      if(grepl("t2m|tas|temp",tolower(param))) {
        param <- "air_temperature_2m"
      } else if(grepl("pr|rr",tolower(param))) {
        param <- "precipitation_amount"
      } else if(grepl("slp|psl|pressure",tolower(param))) {
        param <- "air_pressure_at_sea_level"
      } else if(grepl("cc|cloud",tolower(param))) {
        param <- "cloud_area_fraction"
      } else if(grepl("ws",tolower(param))) {
        param <- "wind_speed"
      } else if(grepl("wd",tolower(param))) {
        param <- "wind_direction"
      } else if(grepl("rh",tolower(param))) {
        param <- "relative_humidity"
      } else if(grepl("sw|ssr",tolower(param))) {
        param <- "integral_of_surface_downwelling_shortwave_flux_in_air_wrt_time"
      }
    }
  }
  Z <- NULL
  index.Z <- NULL
  i <- 1
  for (f in files) {
    if (verbose) print(f)
    X <- retrieve.meps(ncfile=file.path(path,f), param=param,
                       lon=lon, lat=lat, #it=it,
                       dt=dt, plot=plot, verbose=verbose)
    #Y[[f]] <- X
    if(is.null(Z)) Z <- array(NA, dim=c(length(files)*366, dim(X)[2]))
    Z[i:(i+dim(X)[1]-1),] <- X
    index.Z <- c(index.Z, index(X))
    i <- i + dim(X)[1]
  }
  ok <- apply(Z, 1, function(x) !all(is.na(x)))
  Z <- Z[ok,]
  Z <- zoo(Z, order.by=as.Date(index.Z))
  Z <- attrcp(X, Z, ignore = c("start","count"))
  attr(Z,"history") <- c(history.stamp(),attr(X,"history"))
  class(Z) <- class(X)
  invisible(Z)
}

retrieve.meps <- function(ncfile, path=NULL, param='rr',
                          lon=c(9.5,11.5), lat=c(59,61), it=NULL, 
                          dt=50, verbose=FALSE, plot=FALSE) {
  if (verbose) print('retrieve.meps')
  if (!is.null(path)) ncfile <- file.path(path,ncfile,fsep = .Platform$file.sep)
  X <- NULL
  qf <- NULL
  test <- NULL
  
  ## check if file exists and type of ncfile object
  if (is.character(ncfile)) {
    if(grepl("https://|http://",ncfile)) {
      if(verbose) print(paste("Input",ncfile,"is a url."))
      ncid <- try(nc_open(ncfile))
      if(inherits(ncid,"try-error")) {
        stop(paste0("Sorry, the url '", ncfile,
                    "' could not be opened!"))
      }
    } else if(!file.exists(ncfile)) {
      stop(paste("Sorry, the netcdf file '", ncfile,
                 "' does not exist or the path has not been set correctly!",sep =""))
    } else {
      ncid <- nc_open(ncfile)
    }
  } else if (class(ncfile) == "ncdf4") {
    ncid <- ncfile
  } else {
    stop("ncfile format should be a valid netcdf filename or a netcdf id of class 'ncdf4'")
  }
  ## Read and put attributes in model
  model <- ncatt_get(ncid,0)
  dimnames <- names(ncid$dim)
  varnames <- names(ncid$var)
  latid <- varnames[tolower(varnames) %in% c("lat","latitude")]
  lonid <- varnames[tolower(varnames) %in% c("lon","longitude")]
  altid <- varnames[tolower(varnames) %in% c("lon","longitude")]
  timeid <- dimnames[grep("time",dimnames)]
  vlat <- ncvar_get(ncid, varid=latid)
  vlon <- ncvar_get(ncid, varid=lonid)
  valt <- ncvar_get(ncid, varid=altid)
  time <- ncvar_get(ncid, varid=timeid)
  d <- c(dim(vlat),length(time))
  if(verbose) print(paste0('region: ',
               paste(round(range(vlon),digits=2),collapse='-'),'E/',
               paste(round(range(vlat),digits=2),collapse='-'),'N'))
  vunit <- ncatt_get(ncid, varid=param, attname="units")$value
  longname <- ncatt_get(ncid, varid=param, attname='standard_name')$value
  projid <- varnames[grep("projection",tolower(varnames))]
  proj <- ncatt_get(ncid, varid=projid)
  if(length(lat)>0) {
    lat.rng <- range(lat)
  } else {
    lat.rng <- range(vlat)
  }
  if(length(dim(vlat))==2) {
    latn <- apply(vlat,2,min)
    latx <- apply(vlat,2,min)
  } else {
    latn <- vlat
    latx <- vlat
  }
  suby <- (lat.rng[1] <= latn) & (lat.rng[2] >= latx)
  starty <- min( (1:length(latx))[suby] )
  county <- sum(suby)
  if (sum(suby)==0) stop(paste('retrieve.rcm: problems, the requested latitude range (',
                               lat.rng[1],'-',lat.rng[2],') is not within present data (',
                               min(latn),'-',max(latx),')'))
  if (verbose) print(paste('latitudes:',min(lat.rng),'-',max(lat.rng),
                           'extracted:',min(latn[suby]),'-',max(latn[suby]),
                           'start=',starty,'count=',county))

  if(length(lon)>0) {
    lon.rng <- range(lon)
  } else {
    lon.rng <- range(vlon)
  }
  if(length(dim(vlat))==2) {
    lonn <- apply(vlon,1,min)
    lonx <- apply(vlon,1,min)
  } else {
    lonn <- vlon
    lonx <- vlon
  }
  subx <- (lon.rng[1] <= lonn) & (lon.rng[2] >= lonx)
  startx <- min( (1:length(lonx))[subx] )
  countx <- sum(subx)
  if (sum(subx)==0) stop(paste('retrieve.rcm: problems, the requested longitude range (',
                                lon.rng[1],'-',lon.rng[2],') is not within present data (',
                                min(lonn),'-',max(lonx),')'))
        
  if (verbose) print(paste('longitudes:',min(lon.rng),'-',max(lon.rng),
                           'extracted:',min(lonn[subx]),'-',max(lonn[subx]),
                           'start=',startx,'count=',countx))
  
  # Extract the time information: unit and time origin
  tatt <- ncatt_get( ncid, varid='time' )
  tunit <- tatt[[grep("unit",names(tatt))]]
  tcal <-""
  if(any(grepl("^cale", names(tatt)))) {
    if (verbose) print("Calender found")
    tcal <- tatt$cale
  }
  a <- regexpr("since",tunit)
  torg <- substr(tunit,a + attr(a,'match.length')+1,a + attr(a,'match.length')+10)
  torig <- paste(unlist(strsplit(tunit," "))[3:4],collapse=" ")
  tunit <- tolower(substr(tunit,1,a-2))
  time <- switch(substr(tunit,1,3),
                 'day'=as.Date(time, origin=torg),
                 'mon'=as.Date(julian(as.Date(paste(time%/%12,time%%12+1,'01',sep='-'))), origin=torg),
                 'hou'=strptime(torig,format="%Y-%m-%d %H") + time*3600,
                 'sec'=strptime(torig,format="%Y-%m-%d %H") + time)
  if (((diff(time)>=28) && (diff(time) <= 31 )) | (length(time) == 1)) {
    time <- as.Date(strftime(time, format="%Y-%m-01"))
    if (verbose) print("monthly frequency, saving as Date Y-m-01")
  }
  if (verbose) print(paste(start(time),end(time),sep=' - '))

  startt <- 1
  countt <- length(time)
  if(!is.null(it)) {
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
      startt <- min(which(round(as.Date(time))>=min(it)))
      stoptt <- max(which(round(as.Date(time))<=max(it)))
      countt <- stoptt - startt + 1
    } else if (inherits(it,'POSIXt')) {
      startt <- min(which(time>=min(it)))
      stoptt <- max(which(time<=max(it)))
      countt <- stoptt - startt + 1    
    } else if (sum(is.element(it,1600:2200)) > 0) {
      startt <- min(which(year(time)>=min(it)))
      stoptt <- max(which(year(time)<=max(it)))
      countt <- stoptt - startt + 1
    } else if ( (max(it) <= length(time)) & min(it >= 1) ) {
      startt <- min(it)
      stoptt <- max(it)
      countt <- stoptt - startt + 1
    } else {
      print(paste("unkown format of input it:",it))
    }
  }
  
  # Test the dimensions so that the count does not exceed the array:
  d1 <- d[1]
  if(length(d)==3) {
    d2 <- d[2]; d3 <- d[3]
  } else {
    d2 <- d[1]; d3 <- d[2]
  }
  if (startx > d1) {
    startx <- d1
    warning("retrieve.senorge: points along the longitude exceed data dimensions")
  }
  if (starty > d2) {
    starty <- d2
    warning("retrieve.senorge: points along the latitude exceed data dimensions")
  }
  if (startt > d3) {
    startt <- d3
    warning("retrieve.senorge: points in time exceed data dimensions")
  }
  if (startx < 1) {
    startx <- 1
    warning("retrieve.senorge: points along the longitude exceed data dimensions")
  }
  if (starty < 1) {
    starty <- 1
    warning("retrieve.senorge: points along the latitude exceed data dimensions")
  }
  if (startt < 1) {
    startt <- 1
    warning("retrieve.senorge: points in time exceed data dimensions")
  }
  # Test the dimensions so that the count does not exceed the array:
  if (startx + countx - 1 > d1) {
    countx <- d[1] - startx + 1
    warning("retrieve.senorge: number of points along the longitude exceeds data dimensions")
  }
  if (starty + county - 1 > d2) {
    county <- d[2] - starty + 1
    warning("retrieve.senorge: number of points along the latitude exceeds data dimensions")
  }
  if (startt + countt - 1> d3) {
    countt <- d[3] - startt + 1
    warning("retrieve.senorge: number of points in time exceeds data dimensions")
  }
  
  if(length(d)==3) {
    vlon <- vlon[startx:(startx+countx-1),starty:(starty+county-1)]
    vlat <- vlat[startx:(startx+countx-1),starty:(starty+county-1)]
    time <- time[startt:(startt+countt-1)]
    startxy <- c(startx,starty)
    countxy <- c(countx,county)
  } else {
    startxy <- max(startx,starty)
    countxy <- min(countx,county)
    vlon <- vlon[startxy:(startxy+countxy-1)]
    vlat <- vlat[startxy:(startxy+countxy-1)]
    time <- time[startt:(startt+countt-1)]
  }
  if (verbose) {print(c(startxy,startt)); print(c(countxy,countt))}
  
  ## Workaround because the server won't give me all the time steps at once:
  if(countt<=dt) {
    X <- ncvar_get(ncid, varid=param, start=c(startxy,startt), count=c(countxy,countt))
  } else {
    X <- array(NA, dim=c(countxy, countt))
    for(t in seq(startt,startt+countt-1,dt)) {
      dtt <- min(dt, startt+countt-t)
      X[,,t:(t+dtt-1)] <- ncvar_get(ncid, varid=param, 
          start=c(startxy,t), count=c(countxy,dtt))
    }
  }
  nc_close(ncid)
  
  # If there are less than 3 dimensions, add one dummy dimension.
  if(length(dim(X))==3) {
    d <- dim(X)
  } else if (length(dim(X))!=3 & length(d)==3) {
    d <- c(countxy, countt)
  } else {
    d <- c(1,dim(X))
  }
  dim(X) <- c(d[1]*d[2],d[3])
  
  Y <- zoo(t(X),order.by=time)
  attr(Y,'longitude') <- c(vlon)
  attr(Y,'latitude') <- c(vlat)
  attr(Y,'altitude') <- c(valt)
  attr(Y,'count') <- countt
  attr(Y,'start') <- startt
  attr(Y,'variable') <- param
  attr(Y,'unit') <- vunit
  attr(Y,'source') <- ncfile
  attr(Y,'location') <- rep(NA,length(vlon))
  attr(Y,'longname') <- longname
  attr(Y,'history') <- history.stamp()
  rm('X')
  class(Y) <- c('station','day','zoo')
  
  ## plot the results
  if (plot) map(Y,...)
  invisible(Y)
}


