#' Function to read seNorge daily temperature and precipitation data from thredds.met.no
#' 
#' The seNorge data set offers gridded fields of meteorological variables at a resolution of 1x1km. 
#' The 'senorge' function reads the seNorge data directly from the metno thredds server. 
#' As the resolution is very high, it may be necessary to limit the spatial and/or temporal scope 
#' (see input options 'it', 'lon', and 'lat') as you will otherwise run out of memory.
#' 
#' @aliases senorge, retrieve.senorge
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
#' @seealso meta.thredds, station.thredds
#' 
#' @examples
#' it <- as.Date(c("1962-12-01","1962-12-31"))
#' lon <- c(-2,15)
#' lat <- c(55,63)
#' rr <- senorge(param="rr", lon=lon, lat=lat, it=it, verbose=TRUE)
#' map(rr, FUN="mean")
#' 
#' @export
senorge <- function(url='https://thredds.met.no/thredds/catalog/senorge/seNorge_2018/Archive',
                    param='rr', lon=c(9.5,11.5), lat=c(59,61), 
                    it=2010:2019, dt=50, verbose=FALSE, plot=FALSE) {
  if (verbose) print('esd::senorge')
  path <- sub("/catalog/","/dodsC/",url)
  contents <- readLines(paste0(url,'/catalog.html'))
  files <- contents[grepl("_[0-9]{4}.nc",contents,fixed=FALSE)]
  files <- sort(gsub(".*<tt>|</tt>.*","",files,fixed=FALSE))
  if(!is.null(it)) {
    years <- as.numeric(gsub(".*_|.nc","",files))
    if (inherits(it,c('field','station'))) {
      if (verbose) print('Use time interval from an object')
      y <- it
      it <- c(start(y),end(y))
      rm('y')
    }
    if(is.character(it)) {
      if(levels(factor(nchar(it)))==10) {
        it <- as.Date(it,format="%Y-%m-%d")
      } else if((levels(factor(nchar(it)))==8)) {
        it <- as.Date(it,format="%Y%m%d")
      } else if(levels(factor(nchar(it)))==4) {
        it <- as.Date(c(paste(min(it),'-01-01',sep=''),
                        paste(max(it),'-12-31',sep='')))
      }
    }
    if(is.numeric(it)) {
      files <- files[years %in% it]
    } else if(inherits(it,c("Date","POSIXt"))) {
      files <- files[years %in% as.numeric(strftime(it,format="%Y"))]
    }
  }
  if(length(files)>0) {
    ncid <- nc_open(file.path(path,files[[1]]))
    nc_close(ncid)
    varnames <- names(ncid$var)
    if(!param %in% varnames) {
      if(grepl("t2m|tas|temp",param)) {
        param <- "tg"
      } else if(grepl("pr",param)) {
        param <- "rr"
      }
    }
  }
  Z <- NULL
  #Y <- list()
  index.Z <- NULL
  i <- 1
  for (f in files) {
    if (verbose) print(f)
    X <- retrieve.senorge(ncfile=file.path(path,f), param=param,
                          lon=lon, lat=lat, it=it,
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

retrieve.senorge <- function(ncfile, path=NULL, param='rr',
                             lon=c(9.5,11.5), lat=c(59,61), it=NULL, 
                             dt=50, verbose=FALSE, plot=FALSE) {
  if (verbose) print('retrieve.senorge')
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
  timeid <- dimnames[grep("time",dimnames)]
  vlat <- ncvar_get(ncid, varid=latid)
  vlon <- ncvar_get(ncid, varid=lonid)
  time <- ncvar_get(ncid, varid=timeid)
  d <- c(dim(vlat),length(time))
  if(verbose) print(paste0('region: ',paste(round(range(vlon),digits=2),collapse='-'),'E/',
                paste(round(range(vlat),digits=2),collapse='-'),'N'))
  vunit <- ncatt_get(ncid, varid=param, attname="units")$value
  longname <- ncatt_get(ncid, varid=param, attname='long_name')$value
  utmid <- varnames[grep("utm",tolower(varnames))]
  utm <- ncatt_get(ncid, varid=utmid)
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
  itunit <- (1:length(names(time)))[is.element(substr(names(tatt),1,4),'unit')]
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
      startt <- min(which(round(time)>=min(it)))
      stoptt <- max(which(round(time)<=max(it)))
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
  } else {
    startt <- 1
    countt <- length(time)
    it <- NA
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
  if(countt<=100) {
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
  attr(Y,'count') <- countt
  attr(Y,'start') <- startt
  attr(Y,'altitude') <- rep(NA,length(vlon))
  attr(Y,'variable') <- param
  attr(Y,'unit') <- vunit
  attr(Y,'source') <- ncfile
  attr(Y,'location') <- rep(NA,length(vlon))
  attr(Y,'longname') <- longname
  attr(Y,'history') <- history.stamp()
  rm('X')
  class(Y) <- c('station','day','zoo')
  
  ## plot the results
  if (plot) map(Y)
  invisible(Y)
}


