#' Retrieve field data from a netcdf file.
#' 
#' Retrieve data from a netcdf file and return a zoo field object with
#' attributes.  \code{retrieve} assumes data on a regular lon-lat grid and
#' \code{retrieve.rcm} reads data on irregular (rotated) grid (typically output
#' from RCMs).
#' 
#' @aliases retrieve retrieve.default retrieve.ncdf4 retrieve.rcm
#' retrieve.station retrieve.stationsummary
#'
#' @import ncdf4
#'
#' @seealso summary.ncdf4 check.ncdf4 file.class
#'
#' @param file Name of the existing netCDF file to be opened or an object of class 'ncdf4'. 
#' The full path to the netCDF file can either be included in 'ncfile' or entered as a separate input ('path').
#' @param path Path to netcdf file
#' @param ncid An object of class 'ncdf4'
#' @param stid Station IDs to read with retrieve.station
#' @param loc locations to read with retrieve.station
#' @param lon Numeric value of longitude for the reference point (in decimal
#' degrees East) or a vector containing the range of longitude values in the
#' form of c(lon.min,lon.max)
#' @param lat Numeric value of latitude for the reference point (in decimal
#' degrees North) or a vector containing the range of latitude values in the
#' form of c(lat.min,lat.max)
#' @param lev Numeric value of pressure levels or a vector containing the range
#' of pressure level values in the form of c(lev.min,lev.max)
#' @param alt Altititude for stations to read with retrieve.station. Negative
#' values for reading stations below the altitude. For a range use
#' c(alt.min,alt.max)
#' @param cntr Countries of stations to read with retrieve.station
#' @param it Numerical or date values of time or a vector containing the range
#' of values in the form of c(start,end). Date format should be in the form of
#' "YYYY-MM-DD".
#' @param is Numerical or logical values of spatial indexing for reading
#' station data (retrieve.station).
#' @param param Parameter or element type. There are several core parameters or
#' elements as well as a number of additional parameters. The parameters or
#' elements are: auto = automatic selection.  precip, prcp, pr = Precipitation
#' (mm) tas, tavg = 2m-surface temperature (in degrees Celcius) tmax, tasmax =
#' Maximum temperature (in degrees Celcius) tmin, tasmin = Minimum temperature
#' (in degrees Celcius)
#' @param plot Logical value. if, TRUE provides a map.
#' @param greenwich Logical value. If FALSE, convert longitudes to -180E/180E
#' or centre maps on Greenwich meridian (0 deg E). In other words, when
#' Greenwich == TRUE, the left boundary of a global field is set to Greenwich
#' and not the dateline.
#' @param ncdf.check Logical value. If TRUE, performs a quick check of the
#' ncfile contents
#' @param miss2na Logical value. If TRUE missing values are converted to "NA"
#' @param verbose Logical value defaulting to FALSE. If FALSE, do not display
#' comments (silent mode). If TRUE, displays extra information on progress.
#' @param onebyone Logical value. If TRUE, retrieve.station reads one station
#' at the time rather than reading a block of data which can be demaning if the
#' stations are stored in widely different parts of the netCDF file.
#' @param sort - if TRUE, sort the metadata according to location name
#' @return A "zoo" "field" object with additional attributes used for further
#' processing.
#'
#' @examples
#' 
#' \dontrun{
#'   # Download air surface temperature (tas) for the 'NorESM1-ME' model
#'   # output prepared for 'CMIP5 RCP4.5' and for run 'r1i1p1' from the climate
#'   # explorer web portal (http://climexp.knmi.nl) and store the file into the 
#'   # local machine, e.g. temporary folder '/tmp' (Size ~96Mb) using the following
#'   # command if needed. Otherwise, specify a netcdf file to retrieve data from. 
#'   url <- "http://climexp.knmi.nl/CMIP5/monthly/tas"
#'   noresm <- "tas_Amon_NorESM1-ME_rcp45_000.nc"
#'   download.file(url=file.path(url,noresm), destfile=noresm,
#'                method="auto", quiet=FALSE, mode="w",
#'                cacheOK = TRUE)
#'   # Retrieve the data into "gcm" object
#'   gcm <- retrieve(file=file.path(~,noresm),param="tas",
#'                 lon=c(-20,30),lat=c(40,90),plot=TRUE)
#'   # Download the air surface temperature (tas) for RCP 4.5 scenarios and
#'   # NorESM1-ME model from the climate explorer and store it in destfile. 
#'   # Compute the anomalies
#'   gcm.a <- as.anomaly(gcm,ref=c(1960:2001))
#'   map(gcm.a,projection="sphere")
#' }
#' 
#' @export retrieve
retrieve <- function(file=NULL,...) UseMethod("retrieve")

## Default function
#' @exportS3Method 
#' @export retrieve.default
retrieve.default <- function(file,param="auto",
                             path=NULL,verbose=FALSE,...) {
  ncfile <- file
  if (verbose) print(paste('retrieve.default - param=',param,'in',ncfile))
  
  if (!is.null(path)) ncfile <- file.path(path,ncfile,fsep = .Platform$file.sep)
  X <- NULL
  qf <- NULL
  test <- NULL
  nc <- nc_open(ncfile)
  dimnames <- names(nc$dim)
  varnames <- names(nc$var)
  if (verbose) {print(dimnames); print(varnames)}
  ## REB 2021-04-16: Check if the file contains station data - if it does, use the retrieve.station method
  if (sum(tolower(dimnames) %in% c("stid"))>0) {
    if (verbose) print('Detected station netCDF')
    nc_close(nc)
    Y <- retrieve.station(file=file,param=param,path=path,verbose=verbose,...)
    return(Y)
  } else if (sum(tolower(varnames) %in% c("longitude","latitude"))>1) {
    if (verbose) print('Detected rotated-grid netCDF (RCMs)')
    nc_close(nc)
    Y <- retrieve.rcm(file=file,param=param,path=path,verbose=verbose,...)
    return(Y)
  } else if (verbose) print('Detected ordinary netCDF')
  
  ilon <- tolower(dimnames) %in% c("x","i") | grepl("lon",tolower(dimnames))
  ilat <- tolower(dimnames) %in% c("y","j") | grepl("lat",tolower(dimnames))
  if(any(ilon) & any(ilat)) {
    if (verbose) print(c(dimnames[ilon],dimnames[ilat]))
    lons <- ncvar_get(nc,dimnames[ilon])
    lats <- ncvar_get(nc,dimnames[ilat])
  } else {
    lons <- NULL
    lats <- NULL
  }
  
  if ( (length(dim(lons))==1) & (length(dim(lats))==1) )  {
    if (verbose) print(paste('Regular grid field found',ncfile))
    #nc_close(nc)
    #X <- retrieve.ncdf4(ncfile,param=param,verbose=verbose,...)
    X <- retrieve.ncdf4(nc,param=param,verbose=verbose,...)
  } else {
    if (verbose) print('Irregular grid field found')
    nc_close(nc)
    class.x <- file.class(ncfile)
    if (tolower(class.x$value[1])=='station' | sum(is.element(class.x$dimnames,'stid')) > 0) {
      X <- retrieve.station(ncfile,param=param,verbose=verbose,...)
    } else {
      X <- retrieve.rcm(ncfile,param=param,verbose=verbose,...)
    }
  }
}

## Set retrieve for ncdf4 object
#' @exportS3Method
#' @export retrieve.ncdf4
retrieve.ncdf4 <- function (file, path=NULL , param="auto",
                            lon=NULL, lat=NULL, lev=NULL, is=NULL, it=NULL,
                            miss2na=TRUE, greenwich=FALSE,
                            plot=FALSE, verbose=FALSE, ...)  {
  ncfile <- file
  if(verbose) print("retrieve.ncdf4")
  if (!is.null(path)) ncfile <- file.path(path,ncfile,fsep = .Platform$file.sep)
  
  ## check if file exists and type of ncfile object
  if (is.character(ncfile)) {
    if(grepl("https://|http://",ncfile)) {
      if(verbose) print(paste("Input",ncfile,"is a url."))
      ncid <- try(nc_open(ncfile))
      if(inherits(ncid,"try-error")) {
        stop(paste0("Sorry, the url '", ncfile,"' could not be opened!"))
      }
    } else if(!file.exists(ncfile)) {
      stop(paste("Sorry, the netcdf file '", ncfile,
                 "' does not exist or the path has not been set correctly!",sep =""))
    } else {
      ncid <- nc_open(ncfile, verbose=verbose)
    }
  } else if (class(ncfile) == "ncdf4") {
    ncid <- ncfile
  } else {
    stop("ncfile format should be a valid netcdf filename or a netcdf id of class 'ncdf4'")
  }
  
  #class.x <- file.class(ncfile)
  #if (verbose) {print('Check class'); print(class.x$value)}
  lon.rng  <- lon
  lat.rng  <- lat
  ## KMP 2021-08-24: For consistency, adding the input argument 'is' 
  ## which is standard for the esd package and used in retrieve.rcm
  if (is.list(is)) {
    nms <- names(is)
    iy <- grep("lat", tolower(substr(nms, 1, 3)))
    if (length(iy)>0) lat.rng <- range(is[[iy]])
    ix <- grep("lon", tolower(substr(nms, 1, 3)))
    if (length(ix)>0) lon.rng <- range(is[[ix]])
  }
  lev.rng  <- lev
  time.rng <- it
  
  ## Read and put attributes in model
  model <- ncatt_get(ncid,0)
  ## Get variable attributes in v1
  namevars <- names(ncid$var)
  ## REB 2023-01-18: If unspecified, select the variable that has the most dimensions 
  if (tolower(param) == "auto") {
    if (verbose) print('<<< varid not provided - selecting variable with the most dimensions >>>')
    # if (ncid$nvars > 1) {
    #   i <- length(namevars)
    # } else {
    #   i <- 1
    # }
    i <- 0; dmax <- 0
    for (iv in 1:length(namevars)) {
      if (ncid$var[[iv]]$ndims > dmax) {
        i <- iv; dmax <- ncid$var[[iv]]$ndims
      }
      if (verbose) print(c(iv,ncid$var[[iv]]$ndims,i,dmax))
    }
    param <- names(ncid$var)[i]
    if (verbose) print(paste('selected',param))
    v1 <- ncid$var[[i]] 
  } else {
    v1 <- NULL
    i <- grep(param,namevars)
    v1 <- eval(parse(text=paste("ncid$var[[",i,"]]",sep="")))
    if (is.null(v1)) {
      stop(paste("Variable ",param," could not be found!",sep=""))
    }
  }
  
  ## Get dimensions and dimension names
  dimnames <- rep(NA,v1$ndims)
  for (i in 1:v1$ndim) {
    dimnames[i] <- tolower(v1$dim[[i]]$name)
  }
  ## Get lon, lat, lev, time attr and values and update values if necessary
  ## Longitudes
  ilon <- which(tolower(dimnames) %in% c("x","i") | grepl("lon|ncells",tolower(dimnames)))
  if (length(ilon)==0) {
    ilon <- NULL
  } else if (length(ilon)>1) {
    stop("Error in dim lon")
  }
  if (!is.null(ilon)) {
    lon <- eval(parse(text=paste("v1$dim[[",as.character(ilon),"]]",sep="")))
  } else {
    lon <- NULL
  }
  if (!is.null(ilon)) {
    ilonunit <- grep("unit",names(lon))
    if (length(ilonunit>1)) {
      if (verbose) print(paste("Longitude unit is :",lon$unit,sep=" "))
      lonunit <- eval(parse(text = paste("lon$",names(lon)[ilonunit],sep="")))
      #if (length(grep("degree.*.east",lonunit))<1) {
      if (length(grep("degree.*east|degree.*E",lonunit))<1) {
        stop("'retrieve.ncdf4' is not suited to extract longitude units different from 'degrees_east'")
      }
    }
  }
  
  ## Update longitude values if greenwich is FALSE
  if (!greenwich) {
    id <- lon$vals > 180
    if (sum(id) > 0) {
      if (verbose) print("Convert to non-Greenwich")
      lon$vals[id] <- lon$vals[id] - 360
    }
  } else {
    id <- lon$vals < 0
    if (sum(id) > 0) {
      if (verbose) print("Convert to Greenwich")
      lon$vals[id] <- lon$vals[id] + 360
    }
  }
  
  ## Latitudes
  ilat <- which(tolower(dimnames) %in% c("y","j") | grepl("lat",tolower(dimnames)))
  if (length(ilat) ==0) {
    ilat <- NULL
  } else if (length(ilat) > 1) {
    stop("Error in dim lat")
  }
  if (!is.null(ilat)) {
    lat <- eval(parse(text=paste("v1$dim[[",as.character(ilat),"]]",sep="")))
  } else {
    lat <- NULL
  }
  
  ## Pressure Level if pressure variable / not used for surface variables
  ilev <- grep("lev|hei|expver", dimnames)
  if (length(ilev) ==0) {
    ilev <- NULL
  } else if (length(ilev)>1) {
    stop("Error in dim lev")
  }
  if (!is.null(ilev)) {
    lev <- eval(parse(text=paste("v1$dim[[",as.character(ilev),"]]",sep="")))
  } else {
    lev <- NULL
  }
  
  ## Time
  itime <- grep("tim", dimnames)
  if (length(itime) ==0) {
    itime <- NULL
  } else if (length(itime)>1) {
    stop("Error in dim time")
  }
  if (!is.null(itime)) {
    time <- eval(parse(text=paste("v1$dim[[",as.character(itime),"]]",sep="")))
  } else {
    time <- NULL
  }
  ## Check & update meta data from the data itself
  ncid2 <- check.ncdf4(ncid,param=param,verbose=verbose) 
  if (length(grep("model",ls())) > 0) model <- ncid2$model 
  if (!is.null(itime)) time <- ncid2$time
  rm(ncid2)
  if (verbose) print(model$frequency)
  ## Subselect a spatial and a temporal domain
  ## Single point extraction
  one.cell <- FALSE
  if ((length(lon.rng) == 1) & (length(lat.rng)==1)) {
    lons <- rep(as.vector(lon$vals),length(lat$vals))
    lats <- rep(as.vector(lat$vals),length(lon$vals))
    
    dmin <- distAB(lon=lon.rng,lat=lat.rng,lons=lons,lats=lats) 
    id <- which(dmin==min(dmin,na.rm=TRUE))
    lon.w <- which(lon$vals==lons[id])
    lat.w <- which(lat$vals==lats[id])
    if (verbose) {
      print(paste("Single point extraction"))
      print(paste("Selected nearest grid cell lon :",
                  as.character(lon$vals[lon.w]),lon$unit,sep=" "))
    }  
    one.cell <- TRUE
  }
  
  ## longitude extract range
  if (!is.null(ilon)) {
    if (!is.null(lon.rng)) {
      if (length(lon.rng) > 2) {
        stop("lon.rng should be in the form of c(x1,x2)")
      } else if (length(lon.rng) == 1) {
        lon.w <- which((lon$vals-lon.rng) == min(abs(lon$vals-lon.rng)))
        if (verbose) print(paste("Single point extraction / Selected nearest grid cell lon :",
                                 as.character(lon$vals[lon.w]),lon$unit,sep=" "))
      } else if (length(lon.rng)==2) {
        lon.w <- which((lon$vals >= lon.rng[1]) &
                         (lon$vals <= lon.rng[length(lon.rng)]))
        if (verbose) print(paste("Selected longitudes:",paste(as.character(sort(lon$vals[lon.w])),
                                                              collapse="/"),lon$units,sep=" "))
      }
    } else {
      lon.w <- seq(1,length(lon$vals),1)
    }
    lon$len <- length(lon.w)
  }
  
  ## latitude extract range
  if (!is.null(ilat)) {
    if (!is.null(lat.rng)) {
      if (length(lat.rng) > 2) {
        stop("lat.rng should be in the form of c(y1,y2)")
      } else if (length(lat.rng) == 1) {
        lat.w <- which((lat$vals-lat.rng) == min(abs(lat$vals-lat.rng)))
        if (verbose) print(paste("Single point extraction / Selected nearest grid cell lat :",
                                 as.character(lat$vals[lat.w]),lat$unit,sep=" "))
      } else if (length(lat.rng) == 2) { 
        lat.w <- which((lat$vals >= lat.rng[1]) &
                         (lat$vals <= lat.rng[length(lat.rng)]))
        if (verbose) print(paste("Selected Latitudes:",paste(as.character(lat$vals[lat.w]),
                                                             collapse="/"),lat$units,sep=" "))
      }
    } else {
      lat.w <- seq(1,length(lat$vals),1)
    }
    lat$len <- length(lat.w)
  }
  ## time extract range
  if (!is.null(itime)) {
    if (verbose) print('Select time sequence')
    if (!is.null(time.rng)) {
      if (verbose) print(time.rng)
      if (length(time.rng) > 2) {
        stop("time.rng should be in the form of c(year1,year2)/c(date1,date2)/c(i1,i2)")
      } else if (length(time.rng) == 1) {
        if ( (max(time$vals) >= max(time.rng)) & (min(time$vals) >= min(time.rng)) )
          time.w <- which((time$vals-time.rng) == min(abs(time$vals-time.rng))) else 
            if ( (max(time.rng) <= length(time$vals)) & (min(time.rng) > 1) )
              time.w <- time.rng
            if (verbose) print(paste("Single time extraction:",as.character(time$vals[time.w]),
                                     time$unit,sep=" "))
      } else if (length(time.rng) == 2) {
        if(is.years(time.rng)) {
          if (sum(is.element(time.rng,format.Date(time$vdate,"%Y"))) < 1) {
            stop("Selected time interval is outside the range of the data")
          }
          time.w <- which((format.Date(time$vdate,"%Y") >= time.rng[1]) &
                            (format.Date(time$vdate,"%Y") <= time.rng[length(time.rng)]))
        } else if(is.dates(time.rng)) {
          if(inherits(time.rng, c("POSIXt"))) {
            time.w <- which((format.Date(time$vdate,"%Y%m%d %H") >= format.Date(time.rng[1], "%Y%m%d %H") &
                               (format.Date(time$vdate,"%Y%m%d %H") <= format.Date(time.rng[length(time.rng)], "%Y%m%d %H"))))
            
          } else {
            time.w <- which((format.Date(time$vdate,"%Y%m%d") >= format.Date(time.rng[1], "%Y%m%d") &
                               (format.Date(time$vdate,"%Y%m%d") <= format.Date(time.rng[length(time.rng)], "%Y%m%d"))))
          }
          if(length(time.w)==0) {
            stop("Selected time interval is outside the range of the data")
          }
        } else if ( (max(time.rng) <= length(time$vals)) & (min(time.rng) >= 1) ) { 
          time.w <- time.rng[1]:time.rng[2]
        } else {
          print(length(time$vals))
          stop("Unknown format of time.rng")
        }
        if (verbose) {
          if (model$frequency == "mon") {
            print(paste("Selected time values:",
                        paste(as.character(format.Date(time$vdate[time.w],"%Y-%m")),
                              collapse="/"),model$frequency,sep=" "))
          } else {
            print(paste("Selected time values:",
                        paste(as.character(time$vdate[time.w]),collapse="/"),
                        model$frequency,sep=" "))
          }
        }
        if ((length(grep("time.w",ls())) < 1) | (length(time.w)<1)) {
          stop("No time overlapping with selected time interval")
        }
      }
    } else {
      time.w <- seq(1,length(time$vals),1)
    }
    time$vdate <- time$vdate[time.w]
    time$len <- length(time.w)
  } 
  
  ## level extract range
  if (!is.null(ilev)) {
    if (verbose) { 
      print(dimnames[grep("lev|hei|expver", dimnames)])
      print(lev)
    }
    levid <- dimnames[grep("lev|hei|expver", dimnames)]
    
    ## REB: 2020-07-21: If the extra dimension is 'expver' pick the most recent one by default 
    if ((length(grep('expver',dimnames))>0) & is.null(lev.rng)) {
      lev.rng=lev$vals[1]
      if (verbose) print(paste('exver: ilev=',ilev,'lev.rng=',lev.rng))
    }
    
    if (is.null(lev.rng)) {
      lev.rng <- as.integer(readline(paste("Warning: 'esd-package' cannot handle more than one pressure level,",
                                           " specify one level from the list and type 'Enter'",
                                           paste(param,"(",paste(lev$val,collapse="/"),lev$levelUnit,")",sep=""))))
    } else if (length(lev.rng)>1) {
      lev.rng <- as.integer(readline(paste("Warning: 'esd-package' cannot handle more than one pressure level,",
                                           " enter a single level value from the list and type 'Enter'",
                                           paste(param,"(",paste(lev$val,collapse="/"),lev$levelUnit,")",sep=""))))
    } 
    if (length(lev.rng) == 1) {
      lev.w <- which((lev$vals-lev.rng) == min(abs(lev$vals-lev.rng)))
      if (verbose) print(paste("Single level extraction:",
                               as.character(lev$vals[lev.w]),
                               lev$unit,sep=" "))
    } else {
      lev.w <- seq(1,length(lev$vals),1)
    }
    lev$len <- length(lev.w)
  } else levid <- NULL
  ## KMP 2020-08-25: Check the order of dimensions and use idim and idim2 
  ## to rearrange start, count and val in case the they are not in standard order 
  ## (lon,lat,time)
  ## REB 2020-09-31: Fixed some minor problems reading ERA5-data with the fourth dimension 'expver'
  ##
  dimnames <- names(ncid$dim)
  lonid <- dimnames[tolower(dimnames) %in% c("lon","longitude","nlon")]
  latid <- dimnames[tolower(dimnames) %in% c("lat","latitude","nlat")]
  timeid <- dimnames[grep("time", tolower(dimnames))]
  if(length(timeid)>1 & "time" %in% dimnames) timeid <- "time"  
  idlon <- ncid$dim[[lonid]]$id
  idlat <- ncid$dim[[latid]]$id
  if (!is.null(levid)) idlev <- ncid$dim[[levid]]$id else idlev <- NULL
  idtime <- ncid$dim[[timeid]]$id
  dimids <- ncid$var[[param]]$dimids
  if (verbose) {print(dimnames); print(dimids); print(idlev)}
  if (is.null(ilev)) idim <- sapply(c(idlon,idlat,idtime), function(x) grep(x, dimids)) else
                     idim <- sapply(c(idlon,idlat,idlev,idtime), function(x) grep(x, dimids))
  if (verbose) {print(c(idlon,idlat,idlev,idtime)); print(idim); print('-------')}
  idim2 <- sapply(idim, function(x) grep(x, seq(length(idim))))
  ## Extract values and add Scale Factor and offset if any
  if (verbose) print(paste("Reading data for ",v1$longname,sep=""))
  if ((one.cell) & (!is.null(itime))) {
    if (verbose) print('one.cell & !NULL(itime)')
    if (!is.null(ilev)) {
      if (verbose) print('Found four dimensions')
      start <- c(lon.w,lat.w,lev.w[1],time.w[1])
      count <- c(1,1,length(lev.w),length(time.w))
      val <- ncvar_get(ncid,param,start[idim],count[idim])
      dim(val) <- count[idim]
      val <- aperm(val, idim2)
    } else {
      if (verbose) print('Only three dimensions')
      start <- c(lon.w,lat.w,time.w[1])
      count <- c(1,1,length(time.w))
      val <- ncvar_get(ncid,param,start[idim],count[idim])
      dim(val) <- count[idim]
      val <- aperm(val, idim2)
    }
    lon$vals <- lon$vals[lon.w]
    lat$vals <- lat$vals[lat.w]
  } else if ((!is.null(ilon)) & (!is.null(itime))) {
    if (verbose) print('!NULL(ilon) & !NULL(itime)')
    diff.lon.w <- diff(rank(lon$vals[lon.w]))
    if(any(diff.lon.w!=1)) {
      id2 <- which(diff.lon.w!=1)
    } else id2 <- 1
    if (!is.null(ilev)) {
      if (verbose) print('!is.null(ilev)')
      if ((sum(id) > 0) & (sum(id2)!=0)) { ## & !greenwich   
        if (verbose) print('((sum(id) > 0) & (sum(id2)!=0))')
        count <- c(length(lon.w),length(lat.w),length(lev.w),length(time.w))
        lon.w1 <-lon.w[1:id2]
        lon.w2 <- lon.w[(id2+1):length(lon.w)]
        start1 <- c(lon.w1[1],lat.w[1],lev.w[1],time.w[1])
        count1 <- c(length(lon.w1),length(lat.w),length(lev.w),length(time.w))
        if (verbose) print(rbind(start1,count1))
        if (verbose) {print(idim); print(idim2)}
        val1 <- ncvar_get(ncid,param,start1[idim],count1[idim],collapse_degen=FALSE)
        val1 <- aperm(val1, idim2)
        d1 <- dim(val1)
        dim(val1) <- c(d1[1],prod(d1[2:length(d1)]))
        start2 <- c(lon.w2[1],lat.w[1],lev.w[1],time.w[1])
        count2 <- c(length(lon.w2),length(lat.w),length(lev.w),length(time.w))
        if (verbose) print(rbind(start2,count2))
        val2 <- ncvar_get(ncid,param,start2[idim],count2[idim],collapse_degen=FALSE)
        val2 <- aperm(val2, idim2)
        d2 <- dim(val2)
        dim(val2) <- c(d2[1],prod(d2[2:length(d2)]))
        val <- rbind(val1,val2)
      } else {
        if (verbose) print('!((sum(id) > 0) & (sum(id2)!=0))')
        start <- c(lon.w[1],lat.w[1],lev.w[1],time.w[1])
        count <- c(length(lon.w),length(lat.w),length(lev.w),length(time.w))
        if (verbose) {print(rbind(start,count)); print(ncid$dim); print(idim)}
        check.d <- rep(NA,ncid$ndims)  
        for (jj in 1:ncid$ndims) check.d[jj] <- ncid$dim[[jj]]$len  
        ## REB 2021-02-08: These lines causes problems with ERA5. 'idim' has repeated element 1 & 3 for some 
        ## strange reason. The original line does not seem consistent with the lines within this if-block.
        # val <- ncvar_get(ncid,param,start[idim],count[idim],collapse_degen=FALSE)
        #dim(val) <- count[idim]
        # val <- aperm(val, idim2)
        val <- ncvar_get(ncid,param,start,count,collapse_degen=FALSE)
        dim(val) <- count
      }
      #dim(val) <- count
      lon$vals <- lon$vals[lon.w]
      lon.srt <- order(lon$vals)
      if (sum(diff(lon.srt)!=1)) {
        if (verbose) print("Sort Longitudes") 
        lon$vals <- lon$vals[lon.srt]
      }
      lat$vals <- lat$vals[lat.w]
      lat.srt <- order(lat$vals)
      if (sum(diff(lat.srt)!=1)) {
        if (verbose) print("Sort Latitudes") 
        lat$vals <- lat$vals[lat.srt]
      }
      dim(val) <- count
      val <- val[lon.srt,lat.srt,,]
      dim(val) <- count
    } else {
      if (verbose) print('otherwise')
      if (!is.null(ilev)) {print('HERE IS A PROBLEM...')}
      if ((sum(id) > 0) & (sum(id2)!=0)) { ## & !greenwich
        if (verbose) print('((sum(id) > 0) & (sum(id2)!=0))')
        count <- c(length(lon.w),length(lat.w),length(time.w))
        lon.w1 <-lon.w[1:id2]
        lon.w2 <- lon.w[(id2+1):lon$len]
        start1 <- c(lon.w1[1],lat.w[1],time.w[1])
        count1 <- c(length(lon.w1),length(lat.w),length(time.w))
        val1 <- ncvar_get(ncid,param,start1[idim],count1[idim],collapse_degen=FALSE)
        val1 <- aperm(val1, idim2)
        d1 <- dim(val1)
        dim(val1) <- c(d1[1],prod(d1[2:length(d1)]))
        start2 <- c(lon.w2[1],lat.w[1],time.w[1])
        count2 <- c(length(lon.w2),length(lat.w),length(time.w))
        val2 <- ncvar_get(ncid,param,start2[idim],count2[idim],collapse_degen=FALSE)
        val2 <- aperm(val2, idim2)
        d2 <- dim(val2)
        dim(val2) <- c(d2[1],prod(d2[2:length(d2)]))
        val <- rbind(val1,val2)
        stopifnot((d1[2]==d2[2]) | (d1[3]==d2[3]))
        dim(val) <- count
      } else {
        if (verbose) print('!((sum(id) > 0) & (sum(id2)!=0))')
        start <- c(lon.w[1],lat.w[1],time.w[1])
        count <- c(length(lon.w),length(lat.w),length(time.w))
        if (verbose) print(rbind(start,count))
        val <- ncvar_get(ncid,param,start[idim],count[idim])
        dim(val) <- count[idim]
        val <- aperm(val, idim2)
      }
      #dim(val) <- count
      lon$vals <- lon$vals[lon.w]
      lon.srt <- order(lon$vals)
      if (sum(diff(lon.srt)!=1)!=0) {
        if (verbose) print("Sort Longitudes") 
        lon$vals <- lon$vals[lon.srt]
      } 
      lat$vals <- lat$vals[lat.w]
      lat.srt <- order(lat$vals)
      if (sum(diff(lat.srt)!=1)!=0) {
        if (verbose) print("Sort Latitudes") 
        lat$vals <- lat$vals[lat.srt]
      }
      val <- val[lon.srt,lat.srt,]
    }
  }
  
  
  ## Convert units
  if(verbose) print("Check and convert units")
  iunit <- grep("unit",names(v1))
  if (length(iunit)>0) {
    text=paste("v1$",names(v1)[iunit],sep="")
    units <- eval(parse(text=text))
    # hebe added extra units test for unusual strings
    if (units=="") {
      try(tmp <- grep("unit",names(ncatt_get(ncid,param)),value=TRUE),silent = !verbose)
      if ((!inherits(tmp, "try-error")) & (length(tmp)!=0)) {
        units<-gsub(" ","",eval(parse(text=paste('ncatt_get(ncid, param)$',tmp,sep = ""))))
      }
    }
    if (((units=="K") | (units=="degK")) & !grepl("anom",v1$longname)) {
      val <- val - 273 
      units <- "degC"
    }
    if ((length(grep("pa",tolower(units)))>0) &
        (!grepl("vapo",tolower(v1$longname))) |
        (length(grep("N",tolower(units)))>0)) {
      val <- val/100 
      units <- "hPa"
    }
    ## 
    if ( ((units=="Kg/m^2/s" | units=="kg m-2 s-1")) | max(abs(val), na.rm=TRUE) < 0.001 ) {
      # KMP 2024-05-24: Added a check of the values here. If the maximum value is over 1 and the unit is in kg m2 s-2, 
      ## the data has likely already been converted to mm/day but the unit hasn't been changed. 
      ## Otherwise you would get values on the order of 24*60*60 mm/day which makes no sense.
      if(max(abs(val), na.rm=TRUE) < 1) val <- val * (24*60*60) 
      units <- "mm/day"
    } 
    if (verbose) print(paste("Data converted to unit:",units, sep= " "))
  } else if (max(val,na.rm=TRUE)<0.001) {
    if (verbose) print('Variable is likely to be precipitation intensity!')
    val <- val * (24*60*60)
    units <- "mm/day"
  }
  
  ## replace missval by NA
  if (miss2na) {
    imissval <- grep("miss",names(v1))
    if (length(imissval)>0) {
      text=paste("v1$",names(v1)[imissval],sep="")
      missval <- eval(parse(text=text))
      val[val == missval] <- NA
    }
    if (verbose) print(paste(as.character(sum(is.na(val))),"missing values have been replaced by NA" , sep = " "))
  }
  
  ## 
  if (verbose) print("Done!")
  
  ## Copy "filename" attribute to model attributes
  model$filename <- ncid$filename
  
  ## close netcdf file
  nc_close(ncid)
  
  ## Create output and save attributes to the results # 
  ## KMP 2021-08-24: redefining d to solve problems with 
  ## fields with only one lon, lat or time step
  #d <- dim(val)
  #if(is.null(d)) {
  #  d <- c(length(lon$vals),length(lat$vals),time$len)
  #  d <- d[match(seq(length(d)),c(ilon,ilat,itime))]
  #}
  if(is.null(ilon) | is.null(ilat) | is.null(itime)) {
    d <- dim(val)
  } else {
    d <- c(lon$len,lat$len,time$len)
    # REB&AM 2021-12-17: we dont understand why the next line is needed
    #d <- d[match(seq(length(d)),c(ilon,ilat,itime))]
    itime <- 3
  }
  if (verbose) {print("dimensions"); print(d)}
  
  if (!one.cell) {
    if (is.null(ilev)) {
      #HBE added option for 2-D field at one/single time 
      if ((length(d)==2) & (length(time$vdate)==1)) { 
        d<-c(d[ilon],d[ilat],1)
        dim(val) <- c(d[ilon]*d[ilat],1) 
      } else {
        dim(val) <- c(d[ilon]*d[ilat],d[itime])
      }
    } else {
      if (length(lev.w)==1) {
        dim(val) <- c(d[ilon]*d[ilat],d[itime]) ## AM 10.08.2015 Single level selection
        d <- d[-ilev]
      } else {
        dim(val) <- c(d[ilon]*d[ilat]*d[ilev],d[itime])
        print("Warning: 'esd-package' cannot handle more than one level (or height) - Please select one level to retrieve the data (e.g. lev=1000)")
      }   
    }
  }
  
  ## Create a zoo object z
  if (one.cell) {
    z <- zoo(x=val,order.by=time$vdate)  
  } else {
    z <- zoo(x=t(val),order.by=time$vdate)
  }
  ## Add attributes to z
  if(verbose) print("Adding attributes")
  if (!is.null(v1)) {
    attr(z,"variable") <- param
    attr(z,"longname") <- v1$longname
    attr(z,"units") <- units
    attr(z,"dimensions") <- d
  }  
  if (!is.null(ilon)) {
    attr(z,"longitude") <- lon$vals
    attr(z,"greenwich") <- greenwich
  }
  if (!is.null(ilat)) {
    attr(z,"latitude") <- lat$vals
  }
  if (!is.null(ilev)) {
    attr(z,"level") <- lev$vals[lev.w]
    attr(z,"levelUnit") <- lev$units
  }
  if (!is.null(itime)) {
    attr(z,"calendar") <- time$calendar
  }
  ## Add attributes
  if (is.null(model$project_id)) model$project_id <- 'NA'
  attr(z, "file") <- model$filename
  attr(z, "source")         <- model$project_id
  attr(z, 'timeunit')       <- sub('8760hour','annual',model$frequency) ## REB 2024-03-14: sub(...)
  attr(z, 'frequency')      <- 1
  mattr <- names(model)[!names(model) %in% c(names(attributes(z)),"project_id","filename")]
  for(a in mattr) attr(z, a) <- model[[a]]
  attr(z, "model_history") <- model$history
  attr(z, "call")           <- match.call()
  if(is.null(attr(z,"institution"))) attr(z, "institution") <- NA 
  if(is.null(attr(z,"reference"))) attr(z, "reference") <- NA
  attr(z, "history")  <- history.stamp()
  if (one.cell) {
    class(z) <- c("station",model$frequency,"zoo")
    attr(z,'location') <- 'Grid cell'
  } else { 
    class(z) <- c("field",model$frequency,"zoo")
  }
  ## plot the results
  if (plot) map(z,...)
  ## REB 2022-01-19 minor fix
  if (class(z)[2]=="8760hour") class(z)[2] <- 'annual'
  invisible(z)
  #}  
} 

#' Check netcdf file
#'
#' Check content of netcdf file including parameters, units, and time format (frequency, calendar type).
#'
#' @param ncid an object of the class 'ncdf4'
#' @param param meteorological parameter
#' @param verbose if TRUE print progress
#'
#' @export check.ncdf4
check.ncdf4 <- function(ncid, param="auto", verbose=FALSE) {
  if(verbose) print(paste("check.ncdf4",param))
  qf <- NULL
  if (tolower(param) == "auto") {
    if (ncid$nvars > 1) {
      i <- length(names(ncid$var))
      ##i <- grep(param, names(ncid$var))
      ##if (length(i) == 0) i <- as.integer(readline(paste("Choose variable ",paste(namevars,collapse="/") ,"(from 1 - ",length(namevars), "): ",sep = "")))
      ##if (!is.integer(i)) stop("You should introduce an integer value and at least select one variable") 
    } else i <- 1
    param <- names(ncid$var)[i] # ; rm(i)
    v1 <- ncid$var[[i]] 
  } else {
    v1 <- NULL
    i <- grep(param, names(ncid$var))
    v1 <- eval(parse(text=paste("ncid$var[[",i,"]]",sep="")))
    if (is.null(v1))
      stop(paste("Variable ",param," could not be found !",sep=""))
  } 
  ## Checking : Variable dimensions ...
  ndims <- eval(parse(text=paste("ncid$var[[",i,"]]$ndims",sep="")))
  if (verbose) print(paste('check.ncdf4: ndims=',ndims))
  dimnames <- rep(NA,ndims)
  if (ndims>0) {
    for (j in 1:ndims) {
      dimnames[j] <- eval(parse(text=paste("ncid$var[[",i,"]]$dim[[",
                                           j,"]]$name",sep="")))
    }
    if (verbose) print("Checking Dimensions --> [ok]")
    if (verbose) print(paste(as.character(ndims),
                             " dimension(s) has(have) been found :"))
    if (verbose) print(dimnames)
  } else {
    stop("Checking Dimensions --> [fail]")
    if (verbose) print("The variable has no dimensions. The file may be corrupted!")  
  }
  dimnames <- dimnames
  if (verbose) print(paste('check.ncdf4: dimnames',dimnames))
  ## Get all attributes in model, check and update
  model <- ncatt_get(ncid,0)
  
  if (verbose) print(paste('check.ncdf4: ','CMIP checks..'))
  ## Update CMIP3 attributes to match those of CMIP5 
  mnames <- names(model)
  history <- ncatt_get(ncid,0,"history")
  
  if (ncatt_get(ncid,0,"project_id")$hasatt) {
    if (model$project_id=="IPCC Fourth Assessment") {
      model$project_id <- "CMIP3"
      if (verbose) print("project_id IPCC Fourth Assessment set to CMIP3")
    }
  } else if (length(grep("sres",tolower(history$value)))>0) {
    model$project_id <- "CMIP3"
  } else if (length(grep("rcp",tolower(history$value)))>0) {
    model$project_id <- "CMIP5"
  } else {
    if (verbose) print("project_id is missing from file attributes")
    model$project_id <- NULL
  }
  
  if (!ncatt_get(ncid,0,"model_id")$hasatt &
      (ncatt_get(ncid,0,"project_id")$hasatt | !is.null(model$project_id))) {   
    hist2 <- unlist(strsplit(history$value,split=c(" ")))
    ih <- grep("tas",hist2)
    if (model$project_id=="CMIP3") {
      txt <- hist2[ih[1]] 
      model$model_id <- cmip3.model_id(txt)
    } else if (model$project_id=="CMIP5") {
      txt <- hist2[ih[2]]
      model$model_id <- cmip5.model_id(txt)
    }
  }
  ## print(ncatt_get(ncid,0,"project_id")$hasatt)
  if (ncatt_get(ncid,0,"experiment_id")$hasatt &
      (ncatt_get(ncid,0,"project_id")$hasatt | !is.null(model$project_id))) {
    if (tolower(model$project_id)=="cmip3") {
      txt <- unlist(strsplit(tolower(ncatt_get(ncid,0,"history")$value),split="/"))
      model$experiment_id <- paste(txt[grep("20c3m",txt)],txt[grep("sres",txt)],sep="-")
    }
  }
  ## 
  if (ncatt_get(ncid,0,"title")$hasatt) {
    title <- ncatt_get(ncid,0,"title")$value
    modelid <- unlist(strsplit(title,split=c(" ")))
    if (length(grep('-',tolower(modelid))>0) & is.null(model$model_id)) {
      model$model_id <- modelid[grep('-',modelid)]
    }
    if (length(grep('cmip',tolower(modelid))>0) &
        is.null(model$project_id)) {
      model$project_id <- modelid[grep('cmip',tolower(modelid))]
    }
    ## model$model_id <-modelid[1]
    if (length(grep('rcp',tolower(modelid))>0) &
        is.null(model$experiment_id)) {
      model$experiment_id <- modelid[grep('rcp',tolower(modelid))]
    }
    ## model$experiment_id <-modelid[2:3]
    model$type <- modelid[grep('-',modelid)+1] 
  }
  
  ## END CMIP3 MODEL NAMES UPDATE
  ## 
  ## Checking : Time unit and origin
  ## Get system info
  a <- Sys.info()
  ## Get time dimension / val + attributes
  itime <- grep("tim", tolower(dimnames))
  if (length(itime) == 0) {
    itime <- NULL
  } else if (length(itime) > 1) {
    stop("Error in time dim")
  } else if (length(itime)==1) {
    time <- eval(parse(text=paste("v1$dim[[",as.character(itime),"]]",sep="")))
  }
  
  ## Get time unit and origin
  if (verbose) print(paste('check.ncdf4: ','time unit and origin'))
  tatt <- tolower(names(time))
  itunit <- grep(c("unit"),tatt)
  itorigin <- grep(c("orig"),tatt)
  if (length(itunit)>0) {   
    tunit <- eval(parse(text = paste("time$",tatt[itunit],sep="")))
    if (verbose) print(paste("Time unit has been found in time$unit attribute (",tunit,")",sep=""))
  } else {
    tunit <- NULL
  }
  if (!is.null(tunit)) {
    if (verbose) print("Checking Time Unit --> [ok]")
  } else {
    if (verbose) print("Checking Time Unit --> [fail]")
  }
  if (!is.null(tunit)) {
    if (verbose) print("Time unit and origin detected in time$unit attribute")
    if(grepl("since",tunit)) {
      tunit <- time$units
      tsplit <- unlist(strsplit(tunit,split=" "))
      torigin <- time$origin <- paste(tsplit[3:length(tsplit)],collapse=" ")
      tunit <- time$units <- unlist(strsplit(tunit,split=" "))[1]
    } else if(grepl("%Y%m%d",tunit)) {
      tunit < time$units
      tsplit <- unlist(strsplit(tunit,split=" "))
      torigin <- time$origin <- paste(tsplit[grep("%Y%m%d", tsplit)],collapse=" ")
    } else if(grepl("year",tunit)) {
      tunit < time$units
      torigin <- NULL
    }
    if (verbose) print(paste("Updating time$unit (",time$unit,") and creating time$origin (",time$origin,") attribute",sep= ""))
  } else if (length(itorigin)>0) {   
    torigin <- eval(parse(text = paste("time$",tatt[itorigin],sep="")))
    if (verbose) print(paste("Time origin has been found in time origin attribute and set to:",torigin,sep=" "))
  } else {
    torigin <- NULL
  }
  if (is.null(torigin)) {
    if (verbose) print(paste("Time units:", tunit, " l=", 
                             min(time$vals[is.finite(time$vals)]),"-", 
                             max(time$vals[is.finite(time$vals)])))
    if (verbose) warning("Cannot determine the time origin!")
    if (verbose) warning("Example format: '15-Dec-1949'")
    if (verbose) warning("NCEP reanalysis typically: 01-01-01")
    if (verbose) warning("ERA-40 typically: 1900-01-01")
    torigin <- readline("Please enter a valid time origin: ")
  }
  if (!is.null(torigin)) {
    if (torigin=="1-01-01 00:00:00") {
      if (verbose) print("bad time origin")
      torigin <- "0001-01-01 00:00:00"
      if (verbose) print(paste("Re-setting time origin (",torigin,")",sep=""))
    }
  } else {
    torigin <- readline("Give me the time origin (format='YYYY-MM-DD' as '1949-22-01'):")
    if (!is.null(torigin)) {
      if (verbose) print(paste("Time origin set to =", torigin))
    } else {
      stop("Could not determine the time origin. The processing has been stopped !")
    }
  }
  
  if (!is.null(torigin)) {
    if(!grepl("%Y%m%d", as.character(torigin))) {
      yorigin <- format.Date(as.Date(torigin),format="%Y")
      morigin <- format.Date(as.Date(torigin),format="%m")
      dorigin <- format.Date(as.Date(torigin),format="%d")
      if (as.numeric(yorigin) == 0) {
        if (verbose) warning("There is no year zero (Press et al., Numerical recipies)")
        yorigin <- 0
        if (verbose) print(paste("Warning : Year origin has been set to:",as.character(1900),sep="->"))     
      }
      if (is.na(dorigin)) {
        if (verbose) warning("Warning : Day origin is missing !")
        dorigin <- 1
        if (verbose) warning("Warning : Day origin has been set to:",dorigin)
      }
      if (is.na(morigin)) {
        if (verbose) warning("Warning : Month origin is missing !")
        morigin <- 1
        if (verbose) warning("Warning : Month origin has been set to:",morigin)
      }
      torigin1 <- paste(yorigin,morigin,dorigin,sep="-")
      ## REB 2023-02-16
      if (verbose) print(c(torigin1,torigin))
      if (!is.na(unlist(strsplit(torigin,split=" "))[2])) 
        torigin <- paste(torigin1,unlist(strsplit(torigin,split=" "))[2],sep=" ") else
          torigin <- paste(torigin1,'12:00')
    } 
  }
  
  if (!is.null(torigin)) {
    if (verbose) print("Checking Time Origin --> [ok]")
  } else if (verbose) {
    print("Checking Time Origin --> [fail]")
  }
  
  ## Checking : Frequency
  type <- c("year","season","month","day","hour","minute","second")
  type.abb <- substr(tolower(type),1,3)
  ## Initialize
  freq.att <- NULL
  ifreq <- grep("freq",names(model))
  if (length(ifreq)>0) {  
    frq <- tolower(eval(parse(text=paste0("model$",names(model)[ifreq]))))
    if(!frq %in% type & !is.null(model$timeunit)) frq <- paste0(frq, model$timeunit) 
    frq2 <- sub('hr','hou',sub('[[:digit:]]','',frq))
    itype <- grep(frq2, type)
    if (length(itype>0)) {
      freq.att <- type[itype]
      if(grepl('[0-9]',frq) & !grepl('[0-9]',freq.att)) freq.att <- paste0(gsub("[a-z]","",frq), freq.att)
      if(verbose) print(paste0("Frequency has been found in model$frequency attribute (",freq.att,")"))
      if(verbose) print("Checking Frequency from attribute --> [ok]")
    }
  } else {
    if(verbose) print("Checking Frequency from attribute --> [fail]")
    if(verbose) print("Frequency has not been found in the attributes") 
  }
  ## Checking frequency from data
  frequency <- freq.data <- NULL
  if (length(time$vals) > 1) {
    freq.data <- datafrequency(data=as.vector(time$vals),unit=tunit,verbose=FALSE) 
  } else {
    freq.data <- 'none'
  } 
  
  if (!is.null(freq.data)) {
    if (verbose) print("Checking Frequency from the data --> [ok]")
  } else {
    if (verbose) print("Checking Frequency from the data --> [fail]")
  }
  
  ## Checking Calendar attribute if any, otherwise set to "ordinary"  
  # Possible values for CMIP5 files are : "365_day", "standard", "proleptic_gregorian", "360_day" 
  ## REB 2021-05-06 - CMIP6 also uses the Julian Calendar
  ical <- grep(c("calend"),tatt)
  ## 
  if (length(ical)>0) {
    calendar.att <- eval(parse(text = paste("time$",tatt[ical],sep="")))
    if (verbose) print("Checking Calendar from time attribute --> [ok]") 
    if (verbose) print(paste("Calendar attribute has been found in time$calendar (",time$calendar,")",sep =""))
  } else {
    if (verbose) print("Checking Calendar from time attribute --> [fail]")
    calendar.att <- NULL
    print("Warning : Calendar attribute has not been found in the meta data and will be set automatically.")
  }

  ## Identifying starting and ending dates for the data if possible
  if (!is.null(torigin)) {
    if(!grepl("%Y%m%d",as.character(torigin))) {
      yorigin <- as.numeric(format.Date(torigin,format="%Y"))
      morigin <- as.numeric(format.Date(torigin,format="%m"))
      dorigin <- as.numeric(format.Date(torigin,format="%d"))
      horigin <- as.numeric(format.Date(torigin,format="%H"))
    }
  }
  
  ## Get calendar from attribute if any and create vector of dates vdate
  ## 'hou'=strptime(torig,format="%Y-%m-%d %H") + time*3600
  if (verbose) {
    print('<<< Check time metadata >>>')
    print(range(time$vals)); print(torigin); print(tunit)
  }
  if (!is.null(calendar.att)) {
    if (grepl("gregorian|proleptic_gregorian",calendar.att) | grepl("julian",calendar.att) | grepl("standard",calendar.att)) {
      if(grepl("%Y%m%d",tunit)) {
        t.day <- floor(time$vals)
        t.hr <- 24*(time$vals-t.day)
        time$vdate <- strptime(paste(t.day,t.hr), format="%Y%m%d %H")
      } else {
        if (grepl("mon",tunit)) {
          if (sum(round(diff(time$vals)) > 1) < 1) {
            year1 <- time$vals[1]%/%12 + yorigin
            month1 <- morigin
            torigin1 <- paste(as.character(year1),month1,"01",sep="-")
          }
        }
        fmt_torigin <- gsub("M:[0-9]{2}|M[0-9]{2}", "%S",
                          gsub(" [0-9]{2}:[0-9]{2}| [0-9]{2}[0-9]{2}", " %H:%M", 
                            gsub("[0-9]{4}-[0-9]{2}-[0-9]{2}", "%Y-%m-%d", torigin)))
        
        time$vdate <- switch(substr(tunit,1,3),
          'sec'= strptime(torigin,format=fmt_torigin) + time$vals,
          'min'= strptime(torigin,format=fmt_torigin) + time$vals*60,
          'hou'= strptime(torigin,format=fmt_torigin) + time$vals*60*60,
          'day'= strptime(torigin,format=fmt_torigin) + time$vals*60*60*24,
          'mon'=seq(as.Date(torigin1), length.out=length(time$vals), by='month'),
          'yea'= year(as.Date(torigin)) + time$vals)
        #                     'sec'= strptime(torigin,format="%Y-%m-%d %H%M%S") + time$vals,
        #                     'min'= strptime(torigin,format="%Y-%m-%d %H%M%S") + time$vals*60,
        #                     'hou'= strptime(torigin,format="%Y-%m-%d %H:%M:%S") + time$vals*60*60,
			  #  'day'= strptime(torigin,format="%Y-%m-%d %H:%M") + time$vals*60*60*24,
			  #  'mon'=seq(as.Date(torigin1),length.out=length(time$vals),by='month'),
        #                     'yea'= year(as.Date(torigin)) + time$vals)
      }
    } else if (!is.na(strtoi(substr(calendar.att, 1, 3))) | grepl("noleap|365_day|360_day",calendar.att)) {
      if (verbose) print(paste0(substr(calendar.att,1, 3), "-days model year found in calendar attribute"))
      if (grepl("noleap|365_day",calendar.att)) {
        time$daysayear <- 365
      } else {
        time$daysayear <- as.numeric(substr(calendar.att, 1, 3))
      }
      if (!is.null(time$daysayear)) {
        if(verbose) print(paste("Creating time$daysayear attribute and setting attribute to", time$daysayear))
        if (time$daysayear==365) {
          mndays <- c(31,28,31,30,31,30,31,31,30,31,30,31) # Number of days in each month
        } else if (time$daysayear==360) {
          mndays <- rep(30,12) # Number of days in each month
        } else {
          if(verbose) print('Warning! unknown calendar type')
        }
        if (verbose) print(mndays)
        if (!is.null(mndays)) {
          year1 <- time$vals[1]%/%time$daysayear + yorigin
          month1 <- morigin
          if (sum(diff(time$vals)%/%time$daysayear) > 1) {
            if(verbose) print("Warning : Jumps of years has been found in the time series ")
            qf <- c(qf,"jumps of years found in time series")
          }
          if (time$vals[1]%%time$daysayear > 27) {
            year1 <- year1 + 1
            month1 <- month1 + 1
          } 
          if (month1>12) month1 <- month1 - 12 
          # construct vdate
	  ## KMP 2024-06-17: Trying to solve a problem connected to an unusual(?) time format,
          ## a 365 day calendar with values given in hours rather than days
          #years <- time$vals%/%time$daysayear + yorigin
          #dayofyear <- time$vals%%time$daysayear
          tdays <- switch(substr(tunit,1,3),
                          'sec' = time$vals/(60*60*24),
                          'min' = time$vals/(60*24),
                          'hou' = time$vals/(24),
                          'day'= time$vals,
                          time$vals)
          years <- tdays %/% time$daysayear + yorigin
          dayofyear <- tdays %% time$daysayear
	  months <- findInterval(floor(dayofyear), c(0,cumsum(mndays)),
	                         rightmost.closed=FALSE, left.open=FALSE)
	  ## KMP 2023-06-01: the month calculation below gave wrong results for some values, e.g., dayofyear=0 and dayofyear=59
	  #months <- findInterval(ceiling(dayofyear), c(1,cumsum(mndays)),
	  #	     		  rightmost.closed=TRUE, left.open=TRUE)
          days <- dayofyear - (cumsum(mndays)-mndays)[months] + 1
          if (verbose) {print(freq.data); print(median(days,na.rm=TRUE))}
          if(freq.data=='month') {
            ## KMP 2020-05-04: diff stops retrieve from reading 1 timestep data!
            if (length(months)==1) {
              time$vdate <- as.Date(paste(years,months,"01",sep="-"))
            } else if ((sum(diff(months) > 1) > 1) | (sum(diff(years) > 1) > 1) | (sum(round(abs(diff(days)))>2)) > 1) {
              print("Warning : Jumps in data have been found !")
              print("Warning: Trust the first date and force a continuous vector of dates !")
              time$vdate <- seq(as.Date(paste(as.character(year1),month1,"01",sep="-")), by = "month",length.out=time$len)
              qf <- c(qf,"jumps in data found - continuous vector forced")
            } else {
              time$vdate <- as.Date(paste(years,months,"01",sep="-")) #"15",sep="-"
            }
          } else if(freq.data %in% c('season','year')) {
            time$vdate <- as.Date(paste(years,months,"01",sep="-"))
          } else {
            # KMP 2018-10-23: subdaily
            if (verbose) print('Needs to use the PCICt-package...')
            if(!requireNamespace("PCICt",quietly=TRUE)) {
              stop("Package \"PCICt\" needed to retrieve subdaily 360-day calendar data. Please install it.")
            }
            if(length(days)==1) {
              if (verbose) print('Only one day')
              time$vdate <- PCICt::as.PCICt(paste(years,months,floor(days),sep="-"),cal=time$daysayear)
            } else if(median(diff(days),na.rm=TRUE)<1) {
              if (verbose) print('Sub-daily')
              hours <- (days-floor(days))*24
              days <- floor(days)
              time$vdate <- PCICt::as.PCICt(paste(years,months,days,hours,sep=":"),format="%Y:%m:%d:%H",
                                            cal=time$daysayear)
            } else if(median(diff(days),na.rm=TRUE)==1) {
              if (verbose) {print('Daily'); print(table(years)); print(table(months)); print(table(floor(days)))}
              time$vdate <- try(PCICt::as.PCICt(paste(years,months,floor(days),sep="-"),
                                            cal=time$daysayear))
              if (inherits(time$vdate,'try-error')) time$vdate <- seq(1,length(years),by=1)
            } 
          }
        }
      }
    }
    if (verbose) print(paste("Starting date : ",time$vdate[1],
                             "Ending date : ",time$vdate[length(time$vdate)], sep = " "))
  } else {
    if (verbose) print("warnings : Automatic detection of the calendar")
    calendar.detect <- "auto"
    if (grepl("sec",tunit)) time$vdate <- as.POSIXct(torigin,tz='UTC') + time$vals
    if (grepl("min",tunit)) time$vdate <- as.POSIXct(torigin,tz='UTC') + time$vals*60
    if (grepl("hou",tunit)) time$vdate <- as.POSIXct(torigin,tz='UTC') + time$vals*60*60
    if (grepl("day",tunit)) {
      if(length(time$vals)==1) {
        time$vdate <- as.Date((time$vals),origin=as.Date(torigin))
      } else if(median(diff(time$vals))<1) {
        time$vdate <- as.POSIXct(torigin,tz='UTC') + time$vals*60*60*24
      } else if(median(diff(time$vals))>=1) {
        time$vdate <- as.Date((time$vals),origin=as.Date(torigin))
      }
    } 
    if (grepl("mon",tunit)) {
      if (length(time$vals)==1) {
        year1 <- time$vals[1]%/%12 + yorigin
        month1 <- morigin
        time$vdate <- as.Date(paste(as.character(year1),month1,"01",sep="-"))
      } else if (sum(diff(time$vals)>1) < 1) {
        year1 <- time$vals[1]%/%12 + yorigin
        month1 <- morigin
        time$vdate <- seq(as.Date(paste(as.character(year1),month1,"01",sep="-")), by = "month",length.out=length(time$vals))
        #seq(as.Date(paste(as.character(year1),month1,"15",sep="-")), by = "month",length.out=length(time$vals))
      } else print("Warning : Monthly data are mangeled") 
    } 
  }
  if ((length(time$vdate)>1) & (grepl("mon",tunit))) {
    if(sum(diff(as.numeric(format.Date(time$vdate,"%m")))>1)) {
      stop("Vector date is mangeled! Need extra check!")
    }
  }
  ## Checking the data / Extra checks / Automatic calendar detection / etc.
  ## Check 1 # Regular frequency
  if (!is.null(time$vdate)) {
    dt <- as.numeric(rownames(table(diff(time$vdate))))
  } else {
    dt <- NULL
  }
  if (!is.null(time$vdate)) {
    if (verbose) print("Vector of date is in the form :")
    if (verbose) print(str(time$vdate))
    if (verbose & length(time$vdate)>1) print(diff(time$vdate))
  } else {
    if(length(ncid$dim$time$vals)==1) {
      dt <- 1
    } else if (grepl("sec",tunit)) {
      dt <- as.numeric(rownames(table(diff(ncid$dim$time$vals/(24*60*60)))))
    } else if (grepl("min",tunit)) {
      dt <- as.numeric(rownames(table(diff(ncid$dim$time$vals/(24*60)))))
    } else if (grepl("day",tunit)) {
      dt <- as.numeric(rownames(table(diff(ncid$dim$time$vals))))
    } else if (grepl("hou",tunit)) { 
      dt <- as.numeric(rownames(table(diff(ncid$dim$time$vals/24))))
    } else if (grepl("mon",tunit)) {
      dt <- as.numeric(rownames(table(diff(ncid$dim$time$vals))))
    }
    if (length(dt)==1) {
      if (verbose) print("Regular frequency has been detected from the data")
    } else if (verbose) print("Irregular frequency has been detected from the data")
    if ((length(dt)==3) & grepl("day",tunit)) {
      if (verbose) print(paste("Calendar is likely to be a 365-",tunit,"with:",as.character(length(dt)),
                               "irregular frequencies",sep = ""))
      dt <- c(28,30,31)
      if (verbose) print(paste(as.character(dt),tunit,sep="-"))
    }
    if ((length(dt)==4) & grepl("day",tunit)) {
      ## 
      if (verbose) print(paste("Calendar is likely to be a Gregorian 365/366-",tunit,"with:",as.character(length(dt)),
                               "irregular frequencies : ",sep=""))
      dt <- c(28,29,30,31)
      if (verbose) print(paste(as.character(dt),tunit,sep="-")) 
    }
    if (!is.null(time$daysayear)) {
      if ((length(dt)==3) & (time$daysayear != 365) & grepl("day",tunit))
        warning("Calendar does not match with detected frequencies")
    }
    if (length(ical)>0)
      if ((length(dt)!=4) & (grepl("gregorian",calendar.att) | grepl("standard",calendar.att)) & grepl("day",tunit))
        warning("Calendar does not match with detected frequencies")
    if ((length(dt)==2) | (length(dt)>4)) {
      if (verbose) print(paste("Warning : Irregular frequencies have been detected - The data might be corrupted and needs extra Checking !"))   
      if (verbose) print(paste(as.character(dt),tunit,sep=" "))
    }
  }
  
  ## End check 1
  ## Begin check 2 if freq.att matches freq.data
  if (length(time$vals)>1) {
    if (!is.null(freq.att)) {
      model$frequency <- freq.att
      if (!is.null(freq.data)) {
        if (match(freq.att,freq.data,nomatch=FALSE)) {
          if (verbose) print("Frequency found in the attribute matches the frequency detected in data")
          model$frequency <- freq.data <- freq.att
        } else {
          print("Warning : Frequency found in the attribute does not match the frequency detected in data")
          model$frequency <- freq.data
          qf <- c(qf,paste("attribute frequency (",freq.att,") does not match data frequency (",freq.data,")",sep=""))
        }
      } 
    } else if (!is.null(freq.data)) {
      model$frequency <- freq.data
    }
  } else if (sum(is.element(tolower(substr(tunit,1,3)), # REB 2016-03-03
                            c('sec','hou','day','mon','yea')))>0) {
    warning(paste('Need to guess the frequency, based on reckognised time units',tunit))
    model$frequency <- 1
  } else {
    stop("Frequency could not be found, neither detected, the data might be corrupted !")
  }
  if (!is.null(model$frequency)) {
    if (verbose) print(paste("Frequency set to ",model$frequency,sep=""))
    if ( (model$frequency %in% c("month","season","year")) & (!is.null(time$vdate)) ) {
      yr <- year(time$vdate)
      mo <- month(time$vdate)
      dy <- "01"
      time$vdate <- as.Date(paste(yr,mo,dy,sep="-"))
    } #else if (model$frequency=="season") {
    #  yr <-year(time$vdate)
    #  mo <- month(time$vdate)
    #  dy <- "15"
    #  time$vdate <- as.Date(paste(yr,mo,dy,sep="-"))     
    #}
  }
  ## End check 2
  if (verbose) print("Checking --> [Done!]")
  model$qf <- qf
  result <- list(model=model,time=time)
  invisible(result)
}

#' @export retrieve.station
retrieve.station <- function(file,param="auto",path=NULL,is=NULL,stid=NULL,loc=NULL,lon=NULL,lat=NULL,it=NULL,
                             alt=NULL,cntr=NULL,start.year.before=NULL,end.year.after=NULL,
                             nmin=NULL,verbose=FALSE,onebyone=FALSE,...) {
  ncfile <- file
  if (verbose) {
    print(match.call())
    print(paste('retrieve.station',ncfile))
  }
  if (!is.null(path)) ncfile <- file.path(path,ncfile,fsep = .Platform$file.sep)
  
  class.x <- file.class(ncfile)
  if (verbose) {
    print('retrieve.station: Check class:')
    print(class.x$value)
  }
  stopifnot(tolower(class.x$value[1])=='station' | length(is.element(class.x$dimnames,'stid')) > 0)
  ncid <- nc_open(ncfile)
  if (param=='auto') { 
    nvars <- length(names(ncid$var))
    varpick <- 1
    while ( ((ncid$var[[varpick]]$ndims==1) | (sum(is.element(names(ncid$var)[varpick],c('loc','cntr','stationID')))==1)) & 
            (varpick <= nvars) ) varpick <- varpick + 1
    if (verbose) print(paste('retrieve.station:',varpick,names(ncid$var)[varpick]))
    param <- names(ncid$var)[varpick]
  }
  size <- eval(parse(text=paste('ncid$var$',param,'$size',sep='')))
  if (verbose) {
    print(paste('retrieve.station: The variable to read is',param))
    print(paste('retrieve.station: Variable size in netCDF file:',paste(size,collapse=' - ')))
  }
  ## Read the metadata:
  tim <- ncvar_get(ncid,'time'); nt <- length(tim)
  stids <- ncvar_get(ncid,'stationID'); ns <- length(stids)
  if (verbose) {
    print('retrieve.station: Get metadata')
    print(stids)
  }
  tunit <- ncatt_get(ncid,'time','units')
  lons <- ncvar_get(ncid,'lon')
  lats <- ncvar_get(ncid,'lat')
  alts <- ncvar_get(ncid,'alt')
  cntrs <- try(ncvar_get(ncid,'cntr'))
  srcs <- try(ncatt_get(ncid,0,'source')$value)
  if (inherits(cntrs,'try-error')) print('retrieve.stationsummary: Warning - no country information')
  nv <- try(ncvar_get(ncid,'number'))
  if (inherits(nv,'try-error')) print('retrieve.stationsummary: Warning - no valid-data information')
  fyr <- try(ncvar_get(ncid,'first'))
  if (inherits(fyr,'try-error')) print('retrieve.stationsummary: Warning - no start year information')
  lyr <- try(ncvar_get(ncid,'last'))
  if (inherits(lyr,'try-error')) print('retrieve.stationsummary: Warning - no end year information')
  longname <- ncatt_get(ncid,param,'long_name')
  unit <- ncatt_get(ncid,param,'units')
  locs <- try(ncvar_get(ncid,'loc'))
  # if (verbose) testsub <- try(print("retrieve.station: locs <- sub('\xc3','',locs,fixed=TRUE,perl=TRUE)")) else 
  #   testsub <- try(sub('\xc3','',locs,fixed=TRUE))
  # if (!inherits(testsub,'try-error')) locs <- sub('\xc3','',locs,fixed=TRUE) else 
    
  missing <- ncatt_get(ncid,param,'missing_value')
  ## Use the metadata to select the stations to read: there is no need to read
  ## all the stations if only a subset is desired
  if (verbose) print('retrieve.station: Metadata is read - Select selected stations')
  if (is.null(is)) ii <- rep(TRUE,length(stids)) else
    if (!is.logical(is)) {
      ii <- rep(FALSE,length(stids)); ii[is] <- TRUE
    } else ii <- is
  if (!is.null(stid)) ii <- ii & is.element(stids,stid)
  if (!is.null(lon)) ii <- ii & (lons >= min(lon)) & (lons <= max(lon))
  if (!is.null(lat)) ii <- ii & (lats >= min(lat)) & (lats <= max(lat))
  if (!is.null(alt)) { 
    if (length(alt)==2) {
      ii <- ii & (alts >= min(alt)) & (alts >= max(alt)) 
    } else if (alt > 0) {
      ii <- ii & (alts >= alt) 
    } else { 
      ii <- ii & (alts <= abs(alt))
    }
  }
  if (!is.null(loc)) ii <- ii & 
    is.element(tolower(substr(locs,1,nchar(loc))),tolower(loc))
  if (!is.null(cntr)) ii <-ii & 
    is.element(tolower(substr(cntrs,1,nchar(cntr))),tolower(cntr))
  if (verbose) {print('retrieve.station: Read following locations');
    print((1:ns)[ii]); print(locs[ii])}
  if (!is.null(nmin)) ii <- ii & (nv >= nmin)
  if (!is.null(start.year.before)) ii <- ii & (fyr <= start.year.before)
  if (!is.null(end.year.after)) ii <- ii & (lyr >= end.year.after)
  is <- (1:ns)[ii]
  
  if (onebyone) {
    if (verbose) print("retrieve.station: Read stations one by one")
    ## For large files and where the stations are seperated far from each other in the
    ## netCDF file, it may be faster to read the files individually and then combine them 
    ## into one object
    X <- retrieve.station(ncfile,param=param,is=is[1],
                          it=it,start.year.before=start.year.before,
                          end.year.after=end.year.after,
                          nmin=nmin)
    if (length(is)>1) {
      for (ii in is[2:length(is)]) {
        x <- retrieve.station(ncfile,param=param,is=ii,
                              it=it,start.year.before=start.year.before,
                              end.year.after=end.year.after,
                              nmin=nmin)
        X <- combine(X,x)
      }
    }
    return(X)
  }
  
  ## Find the real dates:
  if (sum(tim > 10e7)>0) {
    ## If silly values due to missing data
    print(paste('retrieve.station:',sum(tim > 10e7),'suspect time stamps!'))
    notsuspect <- tim <= 10e7
  } else notsuspect <- rep(TRUE,length(tim))
  if (verbose) {print('retrieve.station: Time information'); print(tunit$value); print(range(tim,na.rm=TRUE))}
  if (length(grep('days since',tunit$value))) 
    t <- as.Date(substr(tunit$value,12,21)) + tim else
      if (length(grep('months since',tunit$value))) 
        t <- seq(as.Date(substr(tunit$value,14,23)),max(tim),'1 month') else
          if (length(grep('years since',tunit$value))) 
            t <- seq(as.Date(substr(tunit$value,14,23)),max(tim),'1 year')
  if (is.null(it)) {
    if (verbose) print('retrieve.station: Read whole record')
    it1 <- 1; it2 <- nt
  } else {
    if (verbose) print('retrieve.station: it is not NULL')
    if ( (is.numeric(it)) & (length(it)==2) ) 
      it <- as.Date(c(paste(it[1],'01-01',sep='-'),(paste(it[2],'12-31',sep='-'))))
    if (is.character(it)) it <- as.Date(it)
    if (verbose) print(paste('retrieve.station: Read selected period',min(it),'-',max(it),
                             'from interval',min(t),max(t)))
    if(it[1]<min(t)) {
      it1 <- 1
    } else {
      it1 <- min(which(t>=it[1] & t<=it[2]))
    }
    if(it[2]>max(t)) {
      it2 <- nt - it1 + 1
    } else {
      it2 <- max(which(t>=it[1] & t<=it[2])) - it1 + 1
    }
    if (verbose) print(paste('retrieve.station:',c(it1,it2)))
    #t <- t[it1:(it1+it2-1)] ## KMP 2022-05-03: this is done on line 1386
  }
  
  if (nt == size[1]) { 
    start <- c(it1, min(is))
    count <- c(it2-it1, max(is) - min(is)+1)
    transpose <- FALSE
  } else { 
    ## The netCDF retrieval is faster when the data is ordered as (space,time)
    start <- c(min(is), it1)
    count <- c(max(is) - min(is)+1, it2)
    transpose <- TRUE
  }
  
  ## Read the actual data:
  if (verbose) {
    print(paste('retrieve.station: reading',param))
    print(paste('retrieve.station: Number of stations in file=',ns,' reading',sum(ii)))
    print('retrieve.station: Confirmation of station IDs:'); print(stids[is])
    print('retrieve.station: start=')
    print(start)
    print('retrieve.station: count=')
    print(count)
    #print(ncid)
  }
  x <- ncvar_get(ncid,param,start=start,count=count)
  ## REB 2022-03-29: needed to add two lines for consistency between x and t.
  it1 <- start[2]; it2 <- start[2]+count[2]-1; it12 <- it1:it2
  if (verbose) {print('retrieve.station: time start & count:'); print(range(it12)); print(length(it12))}
  tim <- tim[it12]
  t <- t[it12]
  if (transpose) x <- t(x)
  nc_close(ncid)
  if (verbose) print('retrieve.station: All data has been extracted from the netCDF file')
  
  if (sum(!notsuspect)>0) {
    if (length(dim(x))==2) x <- x[notsuspect,] else x <- x[notsuspect] 
    t <- t[notsuspect]
    tim <- tim[notsuspect]; nt <- length(tim)
  }
  x[x<=missing$value] <- NA
  
  ## The data matrix is not full and may not necessarily correspond to the selection
  ## Need to remove unwanted stations with station numbers in the range of those selected
  iii <- seq(min((1:ns)[ii]),max((1:ns)[ii]),by=1)
  if (verbose) print(paste(' retrieve.station:  Number of stations read so far',length(iii),
                           ' Total number of stations',length(ii),
                           ' Selected stations=',sum(ii)))
  if (sum(ii)>1) {
    iv <- ii[iii]
    if (length(dim(x))==2) x <- x[,iv] else x <- x[iv]
  } else dim(x) <- NULL
  if (verbose) {
    print(paste('retrieve.station: Dimensions of x is ',paste(dim(x),collapse=' - ')))
    print(summary(c(x))); print(sum(is.finite(x)))
  }
  if (length(dim(x))==2) { 
    nv <- apply(x,1,'nv') 
    jt <- (nv > 0)
  } else if (length(t)>1) jt <- is.finite(x) else jt <- is.finite(t)
  
  if (verbose) print(paste('retrieve.station: Number of valid data points',length(jt), 'remove empty periods'))
  lons <- lons[ii]; lats <- lats[ii]; alts <- alts[ii]; cntrs <- cntrs[ii]
  locs <- locs[ii]; stids <- stids[ii]
  if (verbose) print(paste('retrieve.station: length(t)=',length(t),'length(x)=',length(x),'sum(jt)=',sum(jt)))
  if (length(dim(x))==2) x <- x[jt,] else if (length(t)>1) x <- x[jt]
  tim <- tim[jt]; t <- t[jt]
  
  if (length(t)==1) dim(x) <- c(1,length(x)) else 
    if (is.null(dim(x))) dim(x) <- c(length(x),1)
  if (verbose) print(paste('retrieve.station: Dimensions of x is ',paste(dim(x),collapse=' - '),
                           'and length(t) is',length(t)))
  if (length(t) != dim(x)[1]) {
    print(paste('retrieve.station: Failed sanity check - Dimensions of x is ',paste(dim(x),collapse=' - '),
                'and length(t) is',length(t),' there is a bug in the code'))
    stop('retrieve.station error:')
  }
  
  y <- as.station(zoo(x,order.by=t),loc=locs,lon=lons,lat=lats,alt=alts,
                  cntr = cntrs,stid = stids,longname=longname$value,
                  unit=unit$value,param=param,src=srcs)
  
  if (length(t)>1) {
    if (verbose) print('retrieve.station: Exclude empty time periods')
    if (length(dim(y))==2) iv <- apply(coredata(y),1,FUN='nv') else iv <- nv(y)
    y <- subset(y,it=iv > 0, verbose=verbose)
  }
  if (verbose) print(paste('retrieve.station: as.station',min(index(y)),max(index(y)),'Data size=',length(x),'record length=',length(index(y))))
  if (verbose) print('exit retrieve.station')
  return(y)
}

#' @export retrieve.stationsummary
retrieve.stationsummary <- function(file,path=NULL,stid=NULL,loc=NULL,lon=NULL,lat=NULL,
                                    alt=NULL,cntr=NULL,start.year.before=NULL,end.year.after=NULL,
                                    nmin=NULL,verbose=FALSE,sort=FALSE,...) {
  ncfile <- file
  if (verbose) print(paste('retrieve.stationsummary',ncfile))
  if (!is.null(path)) ncfile <- file.path(path,ncfile,fsep = .Platform$file.sep)
  class.x <- file.class(ncfile)
  if (verbose) {print('Check class'); print(class.x$value)}
  stopifnot(tolower(class.x$value[1])=='station' | length(is.element(class.x$dimnames,'stid')) > 0)
  ncid <- nc_open(ncfile)
  param <- names(ncid$var)[1]
  stats <- names(ncid$var)[grep('summary',names(ncid$var))]
  if (verbose) print(paste('The summary statistics to read is',stats))
  tim <- ncvar_get(ncid,'time'); nt <- length(tim)
  stids <- ncvar_get(ncid,'stationID'); ns <- length(stids)
  if (verbose) {print('Get metadata');print(stids)}
  tunit <- ncatt_get(ncid,'time','units')
  lons <- ncvar_get(ncid,'lon')
  lats <- ncvar_get(ncid,'lat')
  alts <- ncvar_get(ncid,'alt')
  cntrs <- try(ncvar_get(ncid,'cntr'))
  srcs <- try(ncatt_get(ncid,0,'source')$value)
  if (inherits(cntrs,'try-error')) print('retrieve.stationsummary: Warning - no country information')
  nv <- try(ncvar_get(ncid,'number'))
  if (inherits(nv,'try-error')) print('retrieve.stationsummary: Warning - no valid-data information')
  fyr <- try(ncvar_get(ncid,'first'))
  if (inherits(fyr,'try-error')) print('retrieve.stationsummary: Warning - no start year information')
  lyr <- try(ncvar_get(ncid,'last'))
  if (inherits(lyr,'try-error')) print('retrieve.stationsummary: Warning - no end year information')
  daysold <- try(ncvar_get(ncid,'days_old'))
  if (inherits(daysold,'try-error')) 
    {print('retrieve.stationsummary: Warning - no days old information'); daysold <- NULL}
  longname <- ncatt_get(ncid,param,'long_name')
  unit <- ncatt_get(ncid,param,'units')
  locs <- try(ncvar_get(ncid,'loc'))
  ## REB 2023-08-09
  lehr <- try(ncvar_get(ncid,'last_element_highest'))
  if (inherits(lehr,'try-error')) lehr <- NULL
  if (param!='precip') lelr <- try(ncvar_get(ncid,'last_element_lowest')) else
    lelr <- rep(NA,length(lehr))
  if (inherits(lelr,'try-error')) lelr <- NULL
  # if (verbose) testsub <- try(print("retrieve.station: locs <- tolower(sub('\xc3','',locs))")) else 
  #   testsub <- try(sub('\xc3','',locs,fixed=TRUE))
  # if (!inherits(testsub,'try-error')) locs <- tolower(sub('\xc3','',locs,fixed=TRUE))
  
  ## Order alphabetically
  if (verbose) 'Sort alphabetically and in proper order with Scandinavian characters'
  locssrt <- tolower(locs);
  ## KMP 2018-11-02: devtools (run_examples) can only handle ASCII characters so I 
  ## replaced the Scandinavian characters with their escape sequence counterparts
  options(encoding='UTF-8')
  locssrt <- sub("\u00E5",'zzz\u00E5',locssrt)
  locssrt <- sub("\u00E6",'zz\u00E6',locssrt)
  locssrt <- sub("\u00F8",'zzz\u00F8',locssrt)
  locssrt <- sub("\u00E4",'zz\u00E4',locssrt)
  locssrt <- sub("\u00F6",'zzz\u00F6',locssrt)
  srt <- order(locssrt); rm('locssrt')
  locs <- paste(toupper(substr(locs,1,1)),tolower(substr(locs,2,nchar(locs))),sep='')
  
  missing <- ncatt_get(ncid,param,'missing_value')
  y <- data.frame(location=locs,longitude=lons,latitude=lats,altitude=alts,country=cntrs,
                  number.valid=nv,first.year=fyr,last.year=lyr,station.id=stids)
  ## REB 2023-08-09
  if (!is.null(lehr)) y[['lehr']] <- lehr
  if (!is.null(lelr)) y[['lelr']] <- lelr
  if (!is.null(daysold)) y[['daysold']] <- daysold
  attr(y,'variable') <- param
  if (verbose) print('Created data.frame with metadata: y')
  if (length(grep('days since',tunit$value))) 
    t <- as.Date(substr(tunit$value,12,21)) + tim else
      if (length(grep('months since',tunit$value))) 
        t <- seq(as.Date(substr(tunit$value,14,23)),max(tim),'1 month') else
          if (length(grep('years since',tunit$value))) 
            t <- seq(as.Date(substr(tunit$value,14,23)),max(tim),'1 year')
  if (verbose) print(paste('Get the time period',paste(range(t,na.rm=TRUE),collapse=' - ')))
  ok <- ( (t > as.Date('1700-01-01')) & (t < as.Date('2300-01-01')) )
  attr(y,'period') <- range(t[ok])
  attr(y,'longname') <- longname$value
  attr(y,'unit') <- unit$value
  attr(y,'missing_value') <- missing$value
  attr(y,'length') <- length(t)
  attr(y,'source') <- srcs
  if (verbose) print('got metadata attributes: period, unit, missing')
  
  if(verbose) print('Read the summary statistics')
  for (i in 1:length(stats)) {
    if(verbose) print(stats[i])
    sname <- substr(stats[i],9,nchar(stats[i]))
    z <- try(ncvar_get(ncid,stats[i]))
    eval(parse(text=paste('y[["',sname,'"]] <- z',sep='')))
  }
  
  nc_close(ncid)
  if (sort) y <- y[srt,]
  good <- (y$location != "")
  y <- y[good,]
  y$location <- as.character(y$location)
  if (verbose) print('Data is extracted from the netCDF file')
  if(verbose) print(summary(y))
  class(y) <- c("stationsummary","data.frame")
  return(y)
}

# Function that reads data stored on an irregular grid. The data is returned as a 'station' object.
#' @export retrieve.rcm
retrieve.rcm <- function(file,param="auto",...,path=NULL,is=NULL,it=NULL,verbose=FALSE) {
  ncfile <- file
  if(verbose) print("retrieve.rcm")
  if (!is.null(path)) ncfile <- file.path(path,ncfile,fsep = .Platform$file.sep)
  
  if (verbose) print(paste('retrieve ',ncfile))
  ncid <- nc_open(ncfile)
  
  # Extract the time information: unit and time origin
  tatt <- ncatt_get( ncid, varid='time' )
  #if (verbose) print(names(tatt))
  itunit <- (1:length(names(tatt)))[is.element(substr(names(tatt),1,4),'unit')]
  tunit <- tatt[[itunit]]
  tcal <-""
  if (sum(is.element(substr(names(tatt),1,4),'cale'))!=0) {
    if (verbose) print("Calendar found")
    tcal <- tatt$cale
  }
  
  a <- regexpr("since",tunit)
  torg <- substr(tunit,a + attr(a,'match.length')+1,a + attr(a,'match.length')+10)
  torig <- paste(unlist(strsplit(tunit," "))[3:4],collapse=" ")
  tunit <- tolower(substr(tunit,1,a-2))
  if(tolower(param) == "auto") {
    nvars <- length(names(ncid$var))
    if (verbose) print(names(ncid$var))
    varpick <- 1
    while ( ((ncid$var[[varpick]]$ndims==1) | 
             grepl("lon|lat|projection|time|height",names(ncid$var)[varpick])) & 
            (varpick <= nvars) ) varpick <- varpick + 1
    if (verbose) print(paste("Selecting variable",varpick,names(ncid$var)[varpick]))
    param <- names(ncid$var)[varpick]
  }
  if (verbose) print(param)

  # Extract unit etc for the parameter
  vatt <- ncatt_get( ncid, varid=param )
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
  longname <- ncatt_get( ncid, varid=param, attname='long_name')$value
  if (is.null(longname)) {
    longname <- switch(param,'t2m'='temperature','tmax'='maximum temperature','tmin'='minimum temperature',
                       'precip'='precipitation','slp'='mean sea level pressure','pt'='precipitation',
                       'pr'='precipitation')
  }
  # Extract the spatial coordinates:
  vnames <- names(ncid$var)
  ##latid <- vnames[is.element(tolower(substr(vnames,1,3)),'lat')]
  ##lonid <- vnames[is.element(tolower(substr(vnames,1,3)),'lon')]
  #latid <- vnames[grep('lat',tolower(vnames))]
  #lonid <- vnames[grep('lon',tolower(vnames))]
  ## KMP 2016-12-20: grep('lat',...) sometimes finds more than 1 match
  latid <- vnames[tolower(vnames) %in% c("lat","latitude","y")]
  lonid <- vnames[tolower(vnames) %in% c("lon","longitude","x")]
  lat <- ncvar_get(ncid,varid=latid)
  lon <- ncvar_get(ncid,varid=lonid)
  time <- ncvar_get(ncid,varid='time')
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
  if (verbose) print(paste('Find time based on this information:',torg,'and',tunit))
  ## Add a sanity check for poorly designed netCDF files (CARRA)
  if ( (substr(tunit,1,3)=="hou") & (nchar(torg)==10) ) {
    if (verbose) print("Need to add %H to time orgigin")
    torg <- paste(torg,'00')
  }
  if ( (substr(tunit,1,3)=="sec") & (nchar(torg)==10) ) {
    if (verbose) print("Need to add %H:%M:%S to time orgigin")
    torg <- paste(torg,'00:00:00')
  }
  time <- switch(substr(tunit,1,3),
                 'day'=as.Date(time+julian(as.Date(torg))),
                 'mon'=as.Date(julian(as.Date(paste(time%/%12,time%%12+1,'01',sep='-'))) + julian(as.Date(torg))),
                 'hou'=strptime(torg,format="%Y-%m-%d %H") + time*3600,
                 'sec'=strptime(torg,format="%Y-%m-%d %H:%M:%S") + time)
  # next save for later if adding no_leap func
  #if ((tcal %in% c("365_day", "365day", "no_leap", "no leap")) && (any(grepl('hou',tunit))) && ((diff(ttest)>29) && (diff(ttest) <= 31 )) )
  #HBE 2018/1/17 saving POSIX with monthly freq as Date at month start
  #if ( ( diff(time)>=28 && diff(time) <= 31 ) |
  if ( all( diff(time)>=28 & diff(time) <= 31 ) |
       (length(time) == 1) ) {
    time <- as.Date(strftime(time, format="%Y-%m-01"))
    if (verbose) print("monthly frequency, saving as Date Y-m-01")
  #} else if (diff(time)>=360 && diff(time) <= 366) {
  } else if (all(diff(time)>=360 & diff(time) <= 366)) {
    time <- as.Date(strftime(time, format="%Y-%m-01"))
    if (verbose) print("monthly frequency, saving as Date Y-m-01")
  } else 
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
      startt <- min(which(time>=min(it)))#min((1:length(time))[it >= time] )
      stoptt <- max(which(time<=max(it)))#max( (1:length(time))[it <= time] )
      countt <- stoptt - startt + 1#max(it) - startt + 1
    } else if (sum(is.element(it,1600:2200)) > 0) {
      startt <- min(which(year(time)>=min(it)))#min( (1:length(time))[it >= year(time)] )
      stoptt <- max(which(year(time)<=max(it)))#max( (1:length(time))[it <= year(time)] )
      countt <- stoptt - startt + 1#max(it) - startt + 1
    } else if ( (max(it) <= length(time)) & min(it >= 1) ) {
      startt <- min(it)
      countt <- max(it) - startt + 1
    } else {
      print(paste("unkown format of input it:",it))
    }
  } else {
    startt <- 1
    countt <- length(time)
    it <- NA
  }
  
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
  rcm <- ncvar_get(ncid,varid=param,start=start, count=count)
  
  nc_close( ncid )
  if(length(dim(rcm))==3) {
    d <- dim(rcm)
  } else if (length(dim(rcm))!=3 & length(d)==3) {
    ## If there are less than 3 dimensions, add one dummy dimension. To avoid crashes...
    D <- rep(1,3)
    if(any(count>1)) {
      nm <- (1:3)[count>1]
      D[nm] <- d[nm]
    }
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
  RCM <- as.station(zoo(t(rcm), order.by=time), 
                    param=param, unit=vunit, 
                    lon=lon, lat=lat, 
                    longname=longname, src=ncfile,
                    verbose=FALSE)
  #RCM <- zoo(t(rcm),order.by=time)
  #attr(RCM,'longitude') <- c(lon)
  #attr(RCM,'latitude') <- c(lat)
  #attr(RCM,'count') <- count
  #attr(RCM,'start') <- start
  #attr(RCM,'altitude') <- rep(NA,length(lon))
  #attr(RCM,'variable') <- param
  #attr(RCM,'unit') <- vunit
  #attr(RCM,'source') <- ncfile
  #attr(RCM,'location') <- rep(NA,length(lon))
  #attr(RCM,'longname') <- longname
  #attr(RCM,'history') <- history.stamp()
  #class(RCM) <- c('station','day','zoo')
  rm('rcm')
  return(RCM)
}
