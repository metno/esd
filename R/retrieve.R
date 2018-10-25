## Name		: Retrieve
## Description	: S3 method used to retrieve data from netcdf file. The default method was set to the old version of retrieve i.e. retrieve.nc in previous clim.pact package 
## Author 	: A. Mezghani, MET NORWAY
## contact 	: abdelkaderm@met.no
## Created      : 21-03-2013
## Last Update	: 25-04-2013 ; 15-05-2014
## require	: zoo, summary_ncdf4 , check_ncdf4
## input	: a zoo field object / 3 dimensional field with dimensions (time,lon,lat)

## depends on both ncdf and ncdf4 library, one of the two must be installed

## Define retrieve as method
retrieve <- function(ncfile=NULL,...) UseMethod("retrieve")

## Default function
retrieve.default <- function(ncfile,param="auto",type="ncdf4",
                             path=NULL,verbose=FALSE,...) {
  if (verbose) print('retrieve.default')
  ## REB 2018-04-06: Add a check for e.g. station data
  class.x <- file.class(ncfile)
  ##
  X <- NULL
  qf <- NULL
  ## 
  ## Setting the path   (sessionInfo()[[1]]$os=='linux-gnu')?
  #if ( (is.null(path))) { 
  #  path <- dirname(ncfile)
  #  ncfile <- basename(ncfile)
  #}   
  ##if (is.character(ncfile)) {
  ##    fext <- substr(ncfile,nchar(ncfile)-1,nchar(ncfile))
  ##    stopifnot(fext=="nc")
  ##}
  
  ## set path
  ## if (!is.null(path)) {
  ##   path <- gsub("[[:punct:]]$","",path)
  ##   path <- gsub('([[:punct:]])\\1+','\\1',path)
  ## } else if (!is.null(ncfile)){
  ##   i <- max(gregexpr("/",ncfile)[[1]])
  ##   if (i>0) {
  ##     path <- substr(ncfile,1,i-1)
  ##     ncfile <- substr(ncfile,i+1,nchar(ncfile))
  ##   } else {
  ##     path <- getwd()
  ##   }
  ## }
  
  test <- NULL
  
  if ((type=="ncdf") | (class(ncfile)=="ncdf")) { ##(library("ncdf",logical.return=TRUE)) {
    # nc <- open.ncdf(file.path(path,ncfile))
    nc <- open.ncdf(ncfile)
    dimnames <- names(nc$dim)
    ilon <- tolower(dimnames) %in% c("x","i") | grepl("lon",tolower(dimnames))
    ilat <- tolower(dimnames) %in% c("y","j") | grepl("lat",tolower(dimnames))
    lon <- ncvar_get(nc,dimnames[ilon])
    lat <- ncvar_get(nc,dimnames[ilat])
    ## KMP 2017-03-13: grep(lon|x|i) picks out everything containing the letters 
    ## lon, x or i (e.g., 'time') so you can easily end up selecting more than one 
    ## longitude dimension. See solution above.
    #lon <- get.var.ncdf(nc,dimnames[grep("lon|x|i",tolower(dimnames))])
    #lat <- get.var.ncdf(nc,dimnames[grep("lat|y|j",tolower(dimnames))])
    close.ncdf(nc)
    if ( (length(dim(lon))==1) & (length(dim(lat))==1) ) {
      if (verbose) print('Regular grid field found')
      X <- retrieve.ncdf(ncfile,path=path,param=param,verbose=verbose,...)
    } else {
      if (verbose) print('Irregular grid field found')
      class.x <- file.class(ncfile)
      if (tolower(class.x$value[1]=='station') | length(is.element(class.x$dimnames,'stid')) > 0)
        X <- retrieve.station(ncfile,path=path,param=param,verbose=verbose,...) else
          X <- retrieve.rcm(ncfile,path=path,param=param,verbose=verbose,...) 
    }
  } else if ((type=="ncdf4") | (class(ncfile)=="ncdf4")) {##(library("ncdf4",logical.return=TRUE)) {
    #nc <- nc_open(file.path(path,ncfile))
    nc <- nc_open(ncfile)
    dimnames <- names(nc$dim)
    ilon <- tolower(dimnames) %in% c("x","i") | grepl("lon",tolower(dimnames))
    ilat <- tolower(dimnames) %in% c("y","j") | grepl("lat",tolower(dimnames))
    if(any(ilon) & any(ilat)) {
      lon <- ncvar_get(nc,dimnames[ilon])
      lat <- ncvar_get(nc,dimnames[ilat])
      ## KMP 2017-03-13: grep(x|i) is too general - identifies any word with x and i.
      #lon <- ncvar_get(nc,dimnames[grep("lon|x|i",tolower(dimnames))])
      #lat <- ncvar_get(nc,dimnames[grep("lat|y|j",tolower(dimnames))])
      nc_close(nc)
    } else {
      lon <- NULL
      lat <- NULL
    }
    if ( (length(dim(lon))==1) & (length(dim(lat))==1) )  {
      if (verbose) print(paste('Regular grid field found',ncfile))
      X <- retrieve.ncdf4(ncfile,path=path,param=param,verbose=verbose,...)
    }
    else {
      if (verbose) print('Irregular grid field found')
      class.x <- file.class(ncfile)
      if (tolower(class.x$value[1])=='station' | length(is.element(class.x$dimnames,'stid')) > 0)
        X <- retrieve.station(ncfile,path=path,param=param,verbose=verbose,...) else
          X <- retrieve.rcm(ncfile,path=path,param=param,verbose=verbose,...)
    }
  } else {
    print("No suitable ncdf or ncdf4 libraries found to read your file or data")
  }
}

## Set retrieve for ncdf4 object
retrieve.ncdf4 <- function (ncfile = ncfile, path = NULL , param = "auto",
                            lon = NULL, lat = NULL, lev = NULL, it = NULL,
                            miss2na = TRUE, greenwich = FALSE ,
                            ##ncdf4.check = TRUE ,
                            plot = FALSE , verbose = FALSE , ...)  {
  ## Begin of function
  ## Update argument names for internal use only
  require(ncdf4) # REB
  ## REB 2018-04-06: Add a check for e.g .station station data
  class.x <- file.class(ncfile)
  ##
  lon.rng  <- lon
  lat.rng  <- lat
  lev.rng  <- lev
  time.rng <- it
  ## set path
  #    if (!is.null(path)) {
  ## AM this line creates pbms for windows users.
  ##path <- gsub("[[:punct:]]$","",path) 
  ## Update netcdf file name
  #        ncfile <- file.path(path,ncfile,fsep = .Platform$file.sep)
  #    }
  
  ## check if file exists and type of ncfile object
  if (is.character(ncfile)) {
    if (!file.exists(ncfile)) {
      stop(paste("Sorry, the netcdf file '", ncfile,
                 "' does not exist or the path has not been set correctly !",sep =""))}
    ncid <- nc_open(ncfile)     
  } else if (class(ncfile) == "ncdf4")
    ncid <- ncfile
  else stop("ncfile format should be a valid netcdf filename or a netcdf id of class 'ncdf4'")  
  ## Read and put attributes in model
  model <- ncatt_get(ncid,0)
  ## Get variable attributes in v1
  namevars <- names(ncid$var)
  if (tolower(param) == "auto") {
    if (ncid$nvars > 1) {
      i <- length(namevars)
      ## print(i)
      ##i <- grep(param, names(ncid$var))
      ##if (length(i) == 0) i <- as.integer(readline(paste("Choose variable ",paste(namevars,collapse="/") ,"(from 1 - ",length(namevars), "): ",sep = "")))
      ##if (!is.integer(i)) stop("You should introduce an integer value and at least select one variable") 
    } else i <- 1
    param <- names(ncid$var)[i] # ; rm(i)
    v1 <- ncid$var[[i]] 
  } else {
    v1 <- NULL
    i <- grep(param,namevars)
    v1 <- eval(parse(text=paste("ncid$var[[",i,"]]",sep="")))
    if (is.null(v1))
      stop(paste("Variable ",param," could not be found !",sep=""))
  }
  ## Get dimensions
  ## Get dimension names
  dimnames <- rep(NA,v1$ndims)
  for (i in 1:v1$ndim)
    dimnames[i] <- tolower(v1$dim[[i]]$name)
  ## Get lon, lat, lev, time attr and values and update values if necessary
  ## Longitudes
  ilon <- which(tolower(dimnames) %in% c("x","i") | grepl("lon|ncells",tolower(dimnames)))
  ## KMP 2017-03-13: grep(x|i) is too general - identifies any word with x or i. 
  #ilon <- grep("lon|x|ncells|i", dimnames)
  if (length(ilon) ==0)
    ilon <- NULL
  else if (length(ilon)>1)
    stop("Error in dim lon")
  if (!is.null(ilon))
    lon <- eval(parse(text=paste("v1$dim[[",as.character(ilon),"]]",sep="")))
  else
    lon <- NULL
  if (!is.null(ilon)) {
    ilonunit <- grep("unit",names(lon))
    if (length(ilonunit>1)) {
      if (verbose) print(paste("Longitude unit is :",lon$unit,sep=" "))
      lonunit <- eval(parse(text = paste("lon$",names(lon)[ilonunit],sep="")))
      if (length(grep("degree.*.east",lonunit))<1)
        stop("'retrieve.ncdf4' is not suited to extract longitude units different from 'degrees_east'")
    }
  }
  ## 
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
  }##else if (!(sum(id) > 0)) lon$vals <- lon$vals + 180
  
  ## Latitudes
  ilat <- which(tolower(dimnames) %in% c("y","j") | grepl("lat",tolower(dimnames)))
  ## KMP 2017-03-13: grep(y|j) is too general - identifies any word with y or j. 
  #ilat <- grep("lat|y|i", dimnames)
  if (length(ilat) ==0)
    ilat <- NULL
  else if (length(ilat) > 1)
    stop("Error in dim lat")
  if (!is.null(ilat))
    lat <- eval(parse(text=paste("v1$dim[[",as.character(ilat),"]]",sep="")))
  else
    lat <- NULL
  ## 
  ## Pressure Level if pressure variable / not used for surface variables
  ilev <- grep("lev|hei", dimnames)
  if (length(ilev) ==0)
    ilev <- NULL
  else if (length(ilev)>1)
    stop("Error in dim lev")
  if (!is.null(ilev))
    lev <- eval(parse(text=paste("v1$dim[[",as.character(ilev),"]]",sep="")))
  else
    lev <- NULL
  ## 
  ## Time
  itime <- grep("tim", dimnames)
  if (length(itime) ==0) itime <- NULL
  else if (length(itime)>1)
    stop("Error in dim time")
  if (!is.null(itime))
    time <- eval(parse(text=paste("v1$dim[[",as.character(itime),"]]",sep="")))
  else
    time <- NULL
  ## Check & update meta data from the data itself
  ncid2 <- check.ncdf4(ncid,param=param,verbose=verbose) 
  if (length(grep("model",ls())) > 0) model <- ncid2$model 
  if (!is.null(itime)) time <- ncid2$time
  rm(ncid2)
  
  if (verbose) print(model$frequency)
  ## Subselect a spatial and a temporal domain
  ##
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
      if (length(lon.rng) > 2)
        stop("lon.rng should be in the form of c(x1,x2)")
      else if (length(lon.rng) == 1) {
        lon.w <- which((lon$vals-lon.rng) == min(abs(lon$vals-lon.rng)))
        if (verbose)
          print(paste("Single point extraction / Selected nearest grid cell lon :",
                      as.character(lon$vals[lon.w]),lon$unit,sep=" "))
      }
      else if (length(lon.rng) == 2)  {
        lon.w <- which((lon$vals >= lon.rng[1]) &
                         (lon$vals <= lon.rng[length(lon.rng)]))
        if (verbose)
          print(paste("Selected longitudes:",paste(as.character(sort(lon$vals[lon.w])),
                                                   collapse="/"),lon$units,sep=" "))
      }
    } else lon.w <- seq(1,length(lon$vals),1)
    ## lon$vals <- as.vector(lon$vals[lon.w])
    lon$len <- length(lon.w)
  }
  
  ## latitude extract range
  if (!is.null(ilat)) {
    if (!is.null(lat.rng)) {
      if (length(lat.rng) > 2) stop("lat.rng should be in the form of c(y1,y2)")
      if (length(lat.rng) == 1) {
        lat.w <- which((lat$vals-lat.rng) == min(abs(lat$vals-lat.rng)))
        if (verbose)
          print(paste("Single point extraction / Selected nearest grid cell lat :",
                      as.character(lat$vals[lat.w]),lat$unit,sep=" "))
      }
      if (length(lat.rng) == 2) { 
        lat.w <- which((lat$vals >= lat.rng[1]) &
                         (lat$vals <= lat.rng[length(lat.rng)]))
        if (verbose)
          print(paste("Selected Latitudes:",paste(as.character(lat$vals[lat.w]),
                                                  collapse="/"),lat$units,sep=" "))
      }
    } else lat.w <- seq(1,length(lat$vals),1)
    ## lat$vals <- as.vector(lat$vals[lat.w])
    lat$len <- length(lat.w)
  }
  
  ## time extract range
  if (!is.null(itime)) {
    if (!is.null(time.rng)) {
      if (length(time.rng) > 2) stop("time.rng should be in the form of c(year1,year2)")
      if (length(time.rng) == 1) {
        time.w <- which((time$vals-time.rng) == min(abs(time$vals-time.rng)))
        if (verbose)
          print(paste("Single time extraction:",as.character(time$vals[time.w]),
                      time$unit,sep=" "))
      }
      if (length(time.rng) == 2) {
        if (sum(is.element(time.rng,format.Date(time$vdate,"%Y"))) < 1)
          stop("Selected time interval is outside the range of the data") 
        time.w <- which((format.Date(time$vdate,"%Y") >= time.rng[1]) &
                          (format.Date(time$vdate,"%Y") <= time.rng[length(time.rng)]))
        if (verbose) {
          if (model$frequency == "mon")
            print(paste("Selected time values:",
                        paste(as.character(format.Date(time$vdate[time.w],"%Y-%m")),
                              collapse="/"),model$frequency,sep=" "))
          else
            print(paste("Selected time values:",
                        paste(as.character(time$vdate[time.w]),collapse="/"),
                        model$frequency,sep=" "))
        }
        if ((length(grep("time.w",ls())) < 1) | (length(time.w)<1))
          stop("No time overlapping with selected time interval")
      }
    } else time.w <- seq(1,length(time$vals),1)
    ## Updating time$vals and time$vdate
    time$vdate <- time$vdate[time.w]
    ## time$vals <- as.vector(time$vals[time.w])
    time$len <- length(time.w)
  } 
  
  ## level extract range
  if (!is.null(ilev)) {
    if (is.null(lev.rng)) {
      lev.rng <- as.integer(readline(paste("Warning: 'esd-package' cannot handle more than one pressure level, specify one level from the list and type 'Enter'",
                                           paste(param,"(",
                                                 paste(lev$val,collapse="/"),
                                                 lev$levelUnit,")",sep=""))))
      if (length(lev.rng)>1)
        lev.rng <- as.integer(readline(paste("Warning: 'esd-package' cannot handle more than one pressure level, enter a single level value from the list and type 'Enter'",
                                             paste(param,"(",
                                                   paste(lev$val,collapse="/"),
                                                   lev$levelUnit,")",sep=""))))
    }
    if (length(lev.rng) > 2)
      stop("lev.rng should be in the form of c(z1,z2)")
    else if (length(lev.rng) == 1) {
      lev.w <- which((lev$vals-lev.rng) == min(abs(lev$vals-lev.rng)))
      if (verbose)
        print(paste("Single level extraction:",
                    as.character(lev$vals[lev.w]),
                    lev$unit,sep=" "))
    } else if (length(lev.rng) == 2) { 
      lev.w <- which((lev$vals >= lev.rng[1]) &
                       (lev$vals <= lev.rng[length(lev.rng)]))
      if (verbose)
        print(paste("Selected Levels:",
                    paste(as.character(lev$vals[lev.w]),
                          collapse="/"),lev$units,sep=" "))
    } else
      lev.w <- seq(1,length(lev$vals),1)
    ## lev$vals <- as.vector(lev$vals[lev.w])
    lev$len <- length(lev.w)
  }
  ##
  ## Extract values and add Scale Factor and offset if any
  if (verbose) print(paste("Reading data for ",v1$longname,sep=""))
  if ((one.cell) & (!is.null(itime))) {
    if (!is.null(ilev)) {
      start <- c(lon.w,lat.w,lev.w[1],time.w[1])
      count <- c(1,1,length(lev.w),length(time.w))
      val <- ncvar_get(ncid,param,start,count)
    } else {
      start <- c(lon.w,lat.w,time.w[1])
      count <- c(1,1,length(time.w))
      val <- ncvar_get(ncid,param,start,count)
    }
    lon$vals <- lon$vals[lon.w]
    lat$vals <- lat$vals[lat.w]
  } else if ((!is.null(ilon)) & (!is.null(itime))) {  
    diff.lon.w <- diff(rank(lon$vals[lon.w]))
    id2 <- which(diff.lon.w!=1)
    if (!is.null(ilev)) {
      if ((sum(id) > 0) & (sum(id2)!=0)) { ## & !greenwich    
        count <- c(length(lon.w),length(lat.w),length(lev.w),length(time.w))
        lon.w1 <-lon.w[1:id2]
        lon.w2 <- lon.w[(id2+1):length(lon.w)]
        start1 <- c(lon.w1[1],lat.w[1],lev.w[1],time.w[1])
        count1 <- c(length(lon.w1),length(lat.w),length(lev.w),length(time.w))
        val1 <- ncvar_get(ncid,param,start1,count1,collapse_degen=FALSE)
        d1 <- dim(val1)
        dim(val1) <- c(d1[1],prod(d1[2:length(d1)]))
        start2 <- c(lon.w2[1],lat.w[1],lev.w[1],time.w[1])
        count2 <- c(length(lon.w2),length(lat.w),length(lev.w),length(time.w))
        val2 <- ncvar_get(ncid,param,start2,count2,collapse_degen=FALSE)
        d2 <- dim(val2)
        dim(val2) <- c(d2[1],prod(d2[2:length(d2)]))
        val <- rbind(val1,val2)
      } else {
        start <- c(lon.w[1],lat.w[1],lev.w[1],time.w[1])
        count <- c(length(lon.w),length(lat.w),length(lev.w),length(time.w))
        val <- ncvar_get(ncid,param,start,count,collapse_degen=FALSE)
      }
      dim(val) <- count
      ## sort longitudes ...
      lon$vals <- lon$vals[lon.w]
      lon.srt <- order(lon$vals)
      if (sum(diff(lon.srt)!=1)) {
        if (verbose) print("Sort Longitudes") 
        ## lon.srt <- order(lon$vals)
        lon$vals <- lon$vals[lon.srt]
      } ## else lon.srt <- seq(1,length(lon$vals),1)
      lat$vals <- lat$vals[lat.w]
      lat.srt <- order(lat$vals)
      if (sum(diff(lat.srt)!=1)) {
        if (verbose) print("Sort Latitudes") 
        ## lat.srt <- order(lat$vals)
        lat$vals <- lat$vals[lat.srt]
      } ## else lat.srt = seq(1,length(lat$vals),1)
      val <- val[lon.srt,lat.srt,,]
      dim(val) <- count
      ##if (length(lev.w)==1) dim(val) <- count[-3]
    }
    else {
      if ((sum(id) > 0) & (sum(id2)!=0)) { ## & !greenwich
        count <- c(length(lon.w),length(lat.w),length(time.w))
        lon.w1 <-lon.w[1:id2]
        lon.w2 <- lon.w[(id2+1):lon$len]
        start1 <- c(lon.w1[1],lat.w[1],time.w[1])
        count1 <- c(length(lon.w1),length(lat.w),length(time.w))
        val1 <- ncvar_get(ncid,param,start1,count1,collapse_degen=FALSE)
        d1 <- dim(val1)
        dim(val1) <- c(d1[1],prod(d1[2:length(d1)]))
        start2 <- c(lon.w2[1],lat.w[1],time.w[1])
        count2 <- c(length(lon.w2),length(lat.w),length(time.w))
        val2 <- ncvar_get(ncid,param,start2,count2,collapse_degen=FALSE)
        d2 <- dim(val2)
        dim(val2) <- c(d2[1],prod(d2[2:length(d2)]))
        val <- rbind(val1,val2)
        stopifnot((d1[2]==d2[2]) | (d1[3]==d2[3]))
        ##dim(val) <- c((d1[1]+d2[1]),d1[2],d2[3])
      } else {
        start <- c(lon.w[1],lat.w[1],time.w[1])
        count <- c(length(lon.w),length(lat.w),length(time.w))
        val <- ncvar_get(ncid,param,start,count)
      }
      dim(val) <- count   
      lon$vals <- lon$vals[lon.w]
      lon.srt <- order(lon$vals)
      if (sum(diff(lon.srt)!=1)!=0) {
        if (verbose) print("Sort Longitudes") 
        ##lon.srt <- order(lon$vals)
        lon$vals <- lon$vals[lon.srt]
      } ## else lon.srt <- seq(1,length(lon$vals),1)
      lat$vals <- lat$vals[lat.w]
      lat.srt <- order(lat$vals)
      if (sum(diff(lat.srt)!=1)!=0) {
        if (verbose) print("Sort Latitudes") 
        ## lat.srt <- order(lat$vals)
        lat$vals <- lat$vals[lat.srt]
      } ## else lat.srt = seq(1,length(lat$vals),1)
      val <- val[lon.srt,lat.srt,]
      ##dim(val) <- count
    }
  }
  ##
  ## Convert units
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
    if ((length(grep("pa",tolower(units)))>0) | (length(grep("N",tolower(units)))>0)) {
      val <- val/100 
      units <- "hPa"
    }
    ## 
    if ((units=="Kg/m^2/s") | (units=="kg m-2 s-1") | (max(abs(val),na.rm=TRUE)<0.001)) {
      val <- val * (24*60*60) 
      units <- "mm/day"
    }
    if (verbose) print(paste("Data converted to unit:",units, sep= " "))
  } else if (max(val,na.rm=TRUE)<0.001) {
    if (verbose) print('Variable is likely to be precipitation intensity !')
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
  if (verbose) print("Done !")
  
  ## Copy "filename" attribute to model attributes
  model$filename <- ncid$filename
  
  ## close netcdf file
  nc_close(ncid)
  
  ## Create output and save attributes to the results # 
  ## 
  ## 
  d <- dim(val)
  if(is.null(d)) {
    d <- c(length(lon$vals),length(lat$vals),length(time$vals))
    d <- d[match(seq(length(d)),c(ilon,ilat,itime))]
  }
  if (verbose) {print("dimensions"); print(d)}
  ##    
  if (!one.cell) {
    if (is.null(ilev)) {
      #HBE added option for 2-D field at one/single time 
      if ((length(d)==2) & (length(time$vdate)==1)) { 
        d<-c(d[ilon],d[ilat],1)
        dim(val) <- c(d[ilon]*d[ilat],1) 
      } else  dim(val) <- c(d[ilon]*d[ilat],d[itime]) 
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
  ## d <- dim(val)
  ##create a zoo object z
  
  if (one.cell) {
    z <- zoo(x=val,order.by=time$vdate)  
  } else {
    ##create a zoo object z
    z <- zoo(x=t(val),order.by=time$vdate)
  }
  ## z <- zoo(x=t(val),order.by=time$vdate)
  ## Add attributes to z
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
  ##attr(z, "project_id")     <- ifelse(!is.null(model$project_id), model$project_id, NA)
  attr(z, "file") <- model$filename
  attr(z, "source")         <- model$project_id
  attr(z, 'timeunit')       <- model$frequency
  attr(z, 'frequency')      <- 1
  ## KMP 2017-03-22: There is a lot of information specific to files 
  ##    from, e.g., CORDEX and CMIP5 that is not passed on to the object here.
  mattr <- names(model)[!names(model) %in% c(names(attributes(z)),"project_id","filename")]
  for(a in mattr) attr(z, a) <- model[[a]]
  ## not needed with the mattr loop:
  #attr(z, "title") <- model$title
  #attr(z, "model_id")       <- model$model_id
  #attr(z, "experiment_id")  <- model$experiment_id
  #attr(z, "realization")    <- model$realization
  #attr(z, "initialization_method") <- model$initialization_method
  #attr(z, "physics_version") <- model$physics_version
  #attr(z, "parent_experiment_rip") <- model$parent_experiment_rip
  #attr(z, 'type')           <- model$type
  ## attr(z, "timestamp")      <- date()
  ## attr(z, "anomaly")        <- FALSE
  ## attr(z, "time:method")    <- NA
  ## attr(z, "spatial:method") <- NA
  attr(z, "URL")            <- "http://climexp.knmi.nl/"
  attr(z, "call")           <- match.call()
  ## attr(z, "history")        <- NA
  if(is.null(attr(z,"institution"))) attr(z, "institution") <- NA 
  if(is.null(attr(z,"reference"))) attr(z, "reference") <- NA
  attr(z, "history")  <- history.stamp()
  if (one.cell) {
    class(z) <- c("station",model$frequency,"zoo")
    attr(z,'location') <- 'Grid cell'
  } else 
    class(z) <- c("field",model$frequency,"zoo")
  
  ## plot the results
  if (plot) map(z,...)
  
  invisible(z)
  
} # End of the function


## Set retrieve for ncdf3 object
retrieve.ncdf <- function (ncfile = ncfile, path = NULL , param = "auto",
                           lon = NULL, lat = NULL, lev = NULL, it = NULL,
                           miss2na = TRUE, greenwich = FALSE , ##ncdf.check = TRUE ,
                           plot = FALSE , verbose = FALSE , ...) {
  
  ## Update argument names for internal function use only
  require(ncdf)
  if (verbose) print('retrieve.ncdf')
  ## REB 2018-04-06: Add a check for e.g. station data
  class.x <- file.class(ncfile,type='ncdf3')
  lon.rng  <- lon
  lat.rng  <- lat
  lev.rng  <- lev
  time.rng <- it
  ##
  
  ## set path
  if (!is.null(path)) {
    ## AM this line creates pbms for windows users.
    ## path <- gsub("[[:punct:]]$","",path)
    ncfile <- file.path(path,ncfile,fsep = .Platform$file.sep)
  }
  
  ## check if file exists and type of ncfile object
  if (is.character(ncfile)) {
    if (!file.exists(ncfile)) {
      stop(paste("Sorry, the netcdf file '", ncfile,
                 "' does not exist or the path has not been set correctly !",sep =""))
    }
    ncid <- open.ncdf(ncfile)     
  } else if (class(ncfile) == "ncdf") {
    ncid <- ncfile
  } else {
    stop("ncfile format should be a valid netcdf filename or a netcdf id of class 'ncdf'")  
  }
  
  ## Read and put attributes in model
  ## model <- att.get.ncdf(ncid,0,"global")
  globalatt <- list('CDI','history','institution','Conventions','references',
                    'institute_id','experiment_id','model_id','forcing',
                    'parent_experiment_id','parent_experiment_rip','branch_time',
                    'contact','initialization_method','physics_version',
                    'tracking_id','version_number','product', 'experiment',
                    'frequency','creation_date','project_id','table_id','title',
                    'parent_experiment','modeling_realm','realization',
                    'cmor_version','CDO','NCO')
  get.att.val <- function(x,nc) return(att.get.ncdf(nc,0,x)$value)
  model <- lapply(as.matrix(globalatt),get.att.val,nc=ncid)
  names(model) <- globalatt
  ##
  
  ## Get variable attributes in v1
  namevars <- names(ncid$var)
  if (tolower(param) == "auto") {
    if (ncid$nvars > 1) {
      i <- length(namevars)
      ## print(i)
      ##i <- grep(param, names(ncid$var))
      ##if (length(i) == 0) i <- as.integer(readline(paste("Choose variable ",paste(namevars,collapse="/") ,"(from 1 - ",length(namevars), "): ",sep = "")))
      ##if (!is.integer(i)) stop("You should introduce an integer value and at least select one variable") i <- grep(param, names(ncid$var))
    } else i <- 1
    param <- names(ncid$var)[i] # ; rm(i)
    v1 <- ncid$var[[i]] 
  } else {
    v1 <- NULL
    i <- grep(param, names(ncid$var))
    v1 <- eval(parse(text=paste("ncid$var[[",i,"]]",sep="")))
    if (is.null(v1)) stop(paste("Variable ",param," could not be found !",sep=""))
  }
  
  ## Get dimensions
  ## Get dimension names
  dimnames <- rep(NA,v1$ndims)
  for (i in 1:v1$ndim) {
    dimnames[i] <- tolower(v1$dim[[i]]$name)
  }
  ## Get lon, lat, lev, time attr and values and update values if necessary
  ## Longitudes
  ilon <- which(tolower(dimnames) %in% c("x","i") | grepl("lon|ncells",tolower(dimnames)))
  ## KMP 2017-03-13: grep(x|i) is too general - identifies any word with x or i. 
  #ilon <- grep("lon|x|i", dimnames)
  if (length(ilon) ==0) {
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
      if (length(grep("degree.*.east",lonunit)) < 1) {
        stop("'retrieve.ncdf4' is not suited to extract longitude units different from 'degrees_east'")
      }
    }
  }
  ## 
  ## Update longitude values if greenwich is FALSE
  if (!greenwich) {
    id <- lon$vals > 180
    if (sum(id) > 0) {
      if (verbose) print("Convert to non-Greenwich as left boundary")
      lon$vals[id] <- lon$vals[id] - 360
    }
  } else {
    id <- lon$vals < 0
    if (sum(id) > 0) {
      if (verbose) print("Convert to Greenwich as left boundary")
      lon$vals[id] <- lon$vals[id] + 360
    }
  }##else if (!(sum(id) > 0)) lon$vals <- lon$vals + 180
  
  ## Latitudes
  ilat <- which(tolower(dimnames) %in% c("y","j") | grepl("lat",tolower(dimnames)))
  ## KMP 2017-03-13: grep(y|j) is too general - will identify any word with an y or j. 
  #ilat <- grep("lat|y|j", dimnames)
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
  ##
  ## Pressure Level if pressure variable / not used for surface variables
  ilev <- grep("lev", dimnames)
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
  ## 
  ## Time 
  itime <- grep("tim", dimnames)
  if (length(itime) ==0) {
    itime <- NULL
  } else if (length(itime)>1) {
    stop("Error in dim time")
  }
  if (!is.null(itime)) {
    time <- eval(parse(text=paste("v1$dim[[",as.character(itime),"]]",sep="")))
  } else time <- NULL
  ## Check and update info 
  ##  
  ##if (ncdf.check) { 
  ncid2 <- F(ncid,param=param,verbose=verbose) 
  if (length(grep("model",ls())) > 0) model <- ncid2$model 
  if (!is.null(itime)) time <- ncid2$time
  rm(ncid2)
  ##}
  ## 
  if (verbose) print(model$frequency)
  ## Subselect a spatial and a temporal domain
  ## longitude extract range
  if (!is.null(ilon)) {
    if (!is.null(lon.rng)) {
      if (length(lon.rng) > 2) stop("lon.rng should be in the form of c(x1,x2)")
      else if (length(lon.rng) == 1) {
        lon.w <- which((lon$vals-lon.rng) == min(abs(lon$vals-lon.rng))) 
        if (verbose) print(paste("Single lon extraction / Selected nearest grid cell to lon :",
                                 as.character(lon$vals[lon.w]),lon$unit,sep=" "))
      } else if (length(lon.rng) == 2) {
        lon.w <- which((lon$vals >= lon.rng[1]) &
                         (lon$vals <= lon.rng[length(lon.rng)]))
        if (verbose) print(paste("Selected longitudes:",
                                 paste(as.character(sort(lon$vals[lon.w])),
                                       collapse="/"),lon$units,sep=" "))
      }
    } else lon.w <- seq(1,length(lon$vals),1)
    ## lon$vals <- as.vector(lon$vals[lon.w])
    lon$len <- length(lon.w)
  }
  
  ## latitude extract range
  if (!is.null(ilat)) {
    if (!is.null(lat.rng)) {
      if (length(lat.rng) > 2)
        stop("lat.rng should be in the form of c(y1,y2)")
      if (length(lat.rng) == 1) {
        lat.w <- which((lat$vals-lat.rng) == min(abs(lat$vals-lat.rng)))
        if (verbose)
          print(paste("Single point extraction / Selected nearest grid cell lat :",
                      as.character(lat$vals[lat.w]),lat$unit,sep=" "))
      }
      if (length(lat.rng) == 2) { 
        lat.w <- which((lat$vals >= lat.rng[1]) & (lat$vals <= lat.rng[length(lat.rng)]))
        if (verbose)
          print(paste("Selected Latitudes:",
                      paste(as.character(lat$vals[lat.w]),
                            collapse="/"),lat$units,sep=" "))
      }
    } else lat.w <- seq(1,length(lat$vals),1)
    ## lat$vals <- as.vector(lat$vals[lat.w])
    lat$len <- length(lat.w)
  }
  
  ##
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
  ## time extract range
  if (!is.null(itime)) {
    if (!is.null(time.rng)) {
      if (length(time.rng) > 2)
        stop("time.rng should be in the form of c(year1,year2)")
      if (length(time.rng) == 1) {
        time.w <- which((time$vals-time.rng) == min(abs(time$vals-time.rng)))
        if (verbose)
          print(paste("Single time extraction:",as.character(time$vals[time.w]),
                      time$unit,sep=" "))
      }
      if (length(time.rng) == 2) {
        if (sum(is.element(time.rng,format.Date(time$vdate,"%Y"))) < 1)
          stop("Selected time interval is outside the range of the data") 
        time.w <- which((format.Date(time$vdate,"%Y") >= time.rng[1]) &
                          (format.Date(time$vdate,"%Y") <= time.rng[length(time.rng)]))
        if (verbose) {
          if (model$frequency == "mon")
            print(paste("Selected time values:",
                        paste(as.character(format.Date(time$vdate[time.w],"%Y-%m")),
                              collapse="/"),model$frequency,sep=" "))
          else
            print(paste("Selected time values:",
                        paste(as.character(time$vdate[time.w]),collapse="/"),
                        model$frequency,sep=" "))
        }
        if ((length(grep("time.w",ls())) < 1) | (length(time.w)<1))
          stop("No time overlapping with selected time interval")
      }
    } else time.w <- seq(1,length(time$vals),1)
    ## Updating time$vals and time$vdate
    time$vdate <- time$vdate[time.w]
    ## time$vals <- as.vector(time$vals[time.w])
    time$len <- length(time.w)
  } 
  ## 
  ## level extract range
  if (!is.null(ilev)) {
    if (is.null(lev.rng))
      lev.rng <- as.integer(readline(paste("Warning: 'esd-package' cannot handle more than one pressure level, enter one level from the list and type 'Enter'",
                                           paste(param,"(",
                                                 paste(lev$val,collapse="/"),
                                                 lev$levelUnit,")",sep=""))))
    else {
      if (length(lev.rng)>1)
        lev.rng <- as.integer(readline(paste("Warning: 'esd-package' cannot handle more than one pressure level, enter one level from the list and type 'Enter'", paste(param,"(",paste(lev$val,collapse="/"),lev$levelUnit,")",sep=""))))
    }
    if (!is.null(lev.rng)) {
      if (length(lev.rng) > 2)
        stop("lev.rng should be in the form of c(z1,z2)")
      if (length(lev.rng) == 1) {
        lev.w <- which((lev$vals-lev.rng) == min(abs(lev$vals-lev.rng)))
        if (verbose) print(paste("Single level extraction:",
                                 as.character(lev$vals[lev.w]),
                                 lev$unit,sep=" "))
      }
      if (length(lev.rng) == 2) { 
        lev.w <- which((lev$vals >= lev.rng[1]) &
                         (lev$vals <= lev.rng[length(lev.rng)]))
        if (verbose)
          print(paste("Selected Levels:",
                      paste(as.character(lev$vals[lev.w]),
                            collapse="/"),lev$units,sep=" "))
      }
    } else lev.w <- seq(1,length(lev$vals),by=1)
    ## lev$vals <- as.vector(lev$vals[lev.w])
    lev$len <- length(lev.w)
  }
  
  ##
  ## Extract values and add Scale Factor and offset if any
  if (verbose) print(paste("Reading data for ",v1$longname,sep=""))
  if ((one.cell) & (!is.null(itime))) {
    if (!is.null(ilev)) {
      start1 <- c(lon.w,lat.w,lev.w[1],time.w[1])
      count1 <- c(1,1,length(lev.w),length(time.w))
      val <- get.var.ncdf(ncid,param,start1,count1)
    } else {
      start1 <- c(lon.w,lat.w,time.w[1])
      count1 <- c(1,1,length(time.w))
      val <- get.var.ncdf(ncid,param,start1,count1)
    }
    lon$vals <- lon$vals[lon.w]
    lat$vals <- lat$vals[lat.w]
  } else if ((!is.null(ilon)) & (!is.null(itime))) {
    diff.lon.w <- diff(rank(lon$vals[lon.w]))
    id2 <- which(diff.lon.w!=1)
    if (!is.null(ilev)) {
      if ((sum(id) > 0) & (sum(id2)!=0)) { ## & !greenwich
        count <- c(length(lon.w),length(lat.w),length(lev.w),
                   length(time.w))
        lon.w1 <- lon.w[1:id2]
        lon.w2 <- lon.w[(id2+1):length(lon.w)]
        start1 <- c(lon.w1[1],lat.w[1],lev.w[1],time.w[1])
        count1 <- c(length(lon.w1),length(lat.w),length(lev.w),
                    length(time.w))
        val1 <- get.var.ncdf(ncid,param,start1,count1)
        ## ,collapse_degen=FALSE)
        d1 <- dim(val1)
        dim(val1) <- c(d1[1],prod(d1[2:length(d1)]))
        start2 <- c(lon.w2[1],lat.w[1],lev.w[1],time.w[1])
        count2 <- c(length(lon.w2),length(lat.w),length(lev.w),
                    length(time.w))
        val2 <- get.var.ncdf(ncid,param,start2,count2)
        ##,collapse_degen=FALSE)
        d2 <- dim(val2)
        dim(val2) <- c(d2[1],prod(d2[2:length(d2)]))
        val <- rbind(val1,val2)
      } else {
        start <- c(lon.w[1],lat.w[1],lev.w[1],time.w[1])
        count <- c(length(lon.w),length(lat.w),length(lev.w),
                   length(time.w))
        val <- get.var.ncdf(ncid,param,start,count)
        ##,collapse_degen=FALSE)
      }
      dim(val) <- count
      ## sort longitudes ...
      lon$vals <- lon$vals[lon.w]
      lon.srt <- order(lon$vals)
      if (sum(diff(lon.srt)!=1)) {
        if (verbose) print("Sort Longitudes") 
        ## lon.srt <- order(lon$vals)
        lon$vals <- lon$vals[lon.srt]
      } ## else lon.srt <- seq(1,length(lon$vals),1)
      lat$vals <- lat$vals[lat.w]
      lat.srt <- order(lat$vals)
      if (sum(diff(lat.srt)!=1)) {
        if (verbose) print("Sort Latitudes") 
        ## lat.srt <- order(lat$vals)
        lat$vals <- lat$vals[lat.srt]
      } ## else lat.srt = seq(1,length(lat$vals),1)
      val <- val[lon.srt,lat.srt,,]
      dim(val) <- count
    }
    else {
      if ((sum(id) > 0) & (sum(id2)!=0)) { ## & !greenwich
        count <- c(length(lon.w),length(lat.w),length(time.w))
        lon.w1 <-lon.w[1:id2]
        lon.w2 <- lon.w[(id2+1):lon$len]
        start1 <- c(lon.w1[1],lat.w[1],time.w[1])
        count1 <- c(length(lon.w1),length(lat.w),length(time.w))
        val1 <- get.var.ncdf(ncid,param,start1,count1)## ,collapse_degen=FALSE)
        dim(val1) <- count1
        d1 <- dim(val1)
        dim(val1) <- c(d1[1],prod(d1[2:length(d1)]))
        start2 <- c(lon.w2[1],lat.w[1],time.w[1])
        count2 <- c(length(lon.w2),length(lat.w),length(time.w))
        val2 <- get.var.ncdf(ncid,param,start2,count2)##,collapse_degen=FALSE)
        dim(val2) <- count2
        d2 <- dim(val2)
        dim(val2) <- c(d2[1],prod(d2[2:length(d2)]))
        val <- rbind(val1,val2)
        stopifnot((d1[2]==d2[2]) | (d1[3]==d2[3]))
        ##dim(val) <- c((d1[1]+d2[1]),d1[2],d2[3])
      } else {
        start <- c(lon.w[1],lat.w[1],time.w[1])
        count <- c(length(lon.w),length(lat.w),length(time.w))
        val <- get.var.ncdf(ncid,param,start,count)
      }
      dim(val) <- count   
      lon$vals <- lon$vals[lon.w]
      lon.srt <- order(lon$vals)
      if (sum(diff(lon.srt)!=1)!=0) {
        if (verbose) print("Sort Longitudes") 
        ##lon.srt <- order(lon$vals)
        lon$vals <- lon$vals[lon.srt]
      } ## else lon.srt <- seq(1,length(lon$vals),1)
      lat$vals <- lat$vals[lat.w]
      lat.srt <- order(lat$vals)
      if (sum(diff(lat.srt)!=1)!=0) {
        if (verbose) print("Sort Latitudes") 
        ## lat.srt <- order(lat$vals)
        lat$vals <- lat$vals[lat.srt]
      } ## else lat.srt = seq(1,length(lat$vals),1)
      val <- val[lon.srt,lat.srt,]
      dim(val) <- count
    }
  }
  
  ## 
  ## Convert units
  iunit <- grep("unit",names(v1))
  if (length(iunit)>0) {
    text=paste("v1$",names(v1)[iunit],sep="")
    units <- eval(parse(text=text))
    if ((units=="K") | (units=="degK") & !grepl("anom",v1$longname)) {
      val <- val - 273 
      units <- "degC"
    }
    if ((length(grep("pa",tolower(units)))>0) | (length(grep("N",tolower(units)))>0)) {
      val <- val/100 
      units <- "hPa"
    }
    ## 
    if ((units=="Kg/m^2/s") | (units=="kg m-2 s-1") | (max(abs(val),na.rm=TRUE)<0.1)) {
      val <- val * (24*60*60) 
      units <- "mm"
    }
    if (verbose) print(paste("Data converted to unit:",units, sep= " "))
  } else if (max(val,na.rm=TRUE)<0.001) {
    if (verbose) print('Variable is likely to be precipitation intensity !')
    val <- val * (24*60*60) 
    units <- "mm"
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
  if (verbose) print("Done !")
  
  ## Copy "filename" attribute to model attributes
  model$filename <- ncid$filename
  
  ## close netcdf4 file
  close.ncdf(ncid)
  
  ## Create output and save attributes to the results # 
  d <- dim(val)
  if (verbose) {print("dimensions")
    print(d)
  }
  ## Convert into 1D or 2D object
  if (!one.cell) {
    if (is.null(ilev)) {
      #HBE added option for 2-D field at one/single time 
      if ((length(d)==2) & (length(time$vdate)==1)) { 
        d<-c(d[ilon],d[ilat],1)
        dim(val) <- c(d[ilon]*d[ilat],1) 
      } else  dim(val) <- c(d[ilon]*d[ilat],d[itime]) 
    } else {
      if (length(lev.w)==1) {
        dim(val) <- c(d[ilon]*d[ilat],d[itime]) ## AM 10.08.2015 Single level selection
        d <- d[-ilev]
      } else {
        dim(val) <- c(d[ilon]*d[ilat]*d[ilev],d[itime])
        print("Warning: 'esd-package' cannot handle more than one level (or heigth) - Please select one level to retrieve the data (e.g. lev=1000)")
      }   
    }
  }
  
  ##create a zoo object z
  if (one.cell) {
    z <- zoo(x=val,order.by=time$vdate)  
  } else {
    z <- zoo(x=t(val),order.by=time$vdate)
  }
  ## d <- dim(val)
  
  ## Add attributes to z
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
    attr(z,"level") <- lev$vals
    attr(z,"levelUnit") <- lev$units
  }
  if (!is.null(itime)) {
    attr(z,"calendar") <- model$calendar
  }
  ## Add attributes
  ##attr(z, "project_id")     <- ifelse(!is.null(model$project_id), model$project_id, NA)
  attr(z, "file") <- model$filename
  attr(z,'source') <- model$project_id
  attr(z,'timeunit') <- model$frequency
  attr(z,'frequency') <- 1
  ## KMP 2017-03-22: There is a lot of information specific to files 
  ##    from, e.g., CORDEX and CMIP5 that is not passed on to the object here.
  mattr <- names(model)[!names(model) %in% c(names(attributes(z)),"filename","project_id")]
  for(a in mattr) attr(z, a) <- model[[a]]
  ## not needed with the mattr loop:
  #attr(z, "title") <- model$title
  #attr(z,'model_id') <- model$model_id
  #attr(z,'experiment_id') <- model$experiment_id
  #attr(z,'realization') <- model$realization
  #attr(z, "initialization_method") <- model$initialization_method
  #attr(z, "physics_version") <- model$physics_version
  #attr(z, "parent_experiment_rip") <- model$parent_experiment_rip
  #attr(z,'type') <- model$type
  ## attr(z, "timestamp")      <- date()
  ## attr(z, "anomaly")        <- FALSE
  ## attr(z, "time:method")    <- NA
  ## attr(z, "spatial:method") <- NA
  ##attr(z, "title")          <- model$title
  attr(z, "URL")            <- "http://climexp.knmi.nl/"
  attr(z, "call")           <- match.call()
  ## attr(z, "history")        <- NA
  if(is.null(attr(z,"institution"))) attr(z, "institution") <- NA 
  if(is.null(attr(z,"reference"))) attr(z, "reference") <- NA
  attr(z, "history")        <- history.stamp()
  
  if (one.cell) {
    class(z) <- c("station",model$frequency,"zoo")
    attr(z,'location') <- 'Grid cell'
  } else 
    class(z) <- c("field",model$frequency,"zoo")
  
  ## plot the results
  if (plot) map(z,...)
  
  invisible(z)
  
} # End of the function

## ------ check.ncdf --------
## Name		: check.ncdf4
## Description	: check the netcdf file attributes and data and update attributes if necessary for further use within ESD package.
## Author 	: A. Mezghani, METNO
## contact 	: abdelkaderm@met.no
## Last Update	: 11-04-2013 ; 17-03-2014
## require	: ncdf4

## Define check as a method
## check <- function(ncid,...) UseMethod("check")

## source("esd/R/frequency.R")

check.ncdf4 <- function(ncid, param="auto",verbose = FALSE) { ## use.cdfcont = FALSE - AM 22-10-2013 not used any more ! 
  
  ## Checking : Number of variables and select only one from the netcdf file, get variable attributes in v1. The user should decide between the list of variables
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
  dimnames <- rep(NA,ndims)
  if (ndims>0) {
    for (j in 1:ndims)
      dimnames[j] <- eval(parse(text=paste("ncid$var[[",i,"]]$dim[[",
                                           j,"]]$name",sep="")))
    if (verbose)
      print("Checking Dimensions --> [ok]")
    if (verbose)
      print(paste(as.character(ndims),
                  " dimension(s) has(have) been found :"))
    if (verbose) print(dimnames)
  } else {
    stop("Checking Dimensions --> [fail]")
    if (verbose) print("The variable has no dimensions. The file may be corrupted!")  
  }
  dimnames <- dimnames
  ## Get all attributes in model, check and update
  model <- ncatt_get(ncid,0)
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
  itime <- grep("tim", dimnames)
  if (length(itime) == 0) {
    itime <- NULL
  } else if (length(itime) > 1) {
    stop("Error in time dim")
  } else if (length(itime)==1) {
    time <- eval(parse(text=paste("v1$dim[[",as.character(itime),"]]",sep="")))
  }
  
  ## Get time unit and origin
  tatt <- tolower(names(time))
  itunit <- grep(c("unit"),tatt)
  itorigin <- grep(c("orig"),tatt)
  if (length(itunit)>0) {   
    tunit <- eval(parse(text = paste("time$",tatt[itunit],sep="")))
    if (verbose) print(paste("Time unit has been found in time$unit attribute (",tunit,")",sep=""))
  } else tunit <- NULL
  if (!is.null(tunit)) {
    if (verbose) print("Checking Time Unit --> [ok]")}
  else if (verbose) print("Checking Time Unit --> [fail]")
  if (!is.null(tunit) & (!is.null(grep("since",tunit)))) {
    if (verbose) print("Time unit and origin detected in time$unit attribute")
    tunit <- time$units
    tsplit <- unlist(strsplit(tunit,split=" "))
    torigin <- time$origin <- paste(tsplit[3:length(tsplit)],collapse=" ")
    tunit <- time$units <- unlist(strsplit(tunit,split=" "))[1]
    if (verbose) print(paste("Updating time$unit (",time$unit,") and creating time$origin (",time$origin,") attribute",sep= ""))
  } else if (length(itorigin)>0) {   
    torigin <- eval(parse(text = paste("time$",tatt[itorigin],sep="")))
    if (verbose) print(paste("Time origin has been found in time origin attribute and set to:",torigin,sep=" "))
  } else torigin <- NULL
  ##if (is.null(torigin) & is.null(tunit)) {
  ##   if (verbose) print("Attributes time unit and origin have not been found -> Sys.time() will be used !")
  ##   if ((tolower(a[1]) == "linux") & (use.cdfcont)) {
  ##      if (verbose) print("Linux & use.cdfcont : call cdfcont()")
  ##      time$origin <- torigin <- cdfcont(ncid$filename)$time.origin
  ##      time$units <- tunit <- cdfcont(ncid$filename)$time.unit
  ##   }
  #}
  if (is.null(torigin)) {
    if (verbose) print(paste("Time units:", tunit, " l=", min(tim[is.finite(time$vals)]),"-", max(tim[is.finite(time$vals)])))
    if (verbose) warning("Cannot determine the time origin!")
    if (verbose) warning("Example format: '15-Dec-1949'")
    if (verbose) warning("NCEP reanalysis typically: 01-01-01")
    if (verbose) warning("ERA-40 typically: 1900-01-01")
    torigin <- readline("Please enter a valid time origin: ")
  }
  
  if (!is.null(torigin)) {
    if (torigin == "1-01-01 00:00:00") {
      if (verbose) print("bad time origin")
      torigin <- "0001-01-01 00:00:00"
      if (verbose) print(paste("Re-setting time origin (",torigin,")",sep=""))
    }
  } else {
    torigin <- readline("Give me the time origin (format='YYYY-MM-DD' as '1949-22-01'):")
    if (!is.null(torigin)) {
      if (verbose) print(paste("Time origin set to =", torigin))
      else stop("Could not determine the time origin. The processing has been stopped !")
    }
  }
  
  if (!is.null(torigin)) {
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
    torigin <- paste(torigin1,unlist(strsplit(torigin,split=" "))[2],sep=" ") 
  }
  
  if (!is.null(torigin)) {
    if (verbose) print("Checking Time Origin --> [ok]")
  } else if (verbose)
    print("Checking Time Origin --> [fail]")
  
  ## Checking : Frequency
  type <- c("year","season","months","Days","hours","minutes","seconds")
  type.abb <- substr(tolower(type),1,3)
  ## Initialize
  freq.att <- NULL
  ifreq <- grep("freq",names(model))
  if (length(ifreq)>0) {  
    itype <- grep(tolower(eval(parse(text=paste("model$",names(model)[ifreq],sep="")))),tolower(type))
    if (length(itype>0)) {
      if (verbose) print(paste("Frequency has been found in model$frequency attribute (",type[itype],")",sep="")) 
      freq.att <- frequency.name[grep(model$frequency,frequency.abb)]
    }
    if (verbose) print("Checking Frequency from attribute --> [ok]")
  } else {
    if (verbose) print("Checking Frequency from attribute --> [fail]")
    if (verbose) print("Frequency has not been found in the attributes") 
  }
  ##  
  ## Checking frequency from data
  frequency <- freq.data <- NULL
  if (length(time$vals) > 1)
    freq.data <- frequency.data(data=as.vector(time$vals),unit=tunit,verbose=FALSE) else
      freq.data <- 'none'
  if (!is.null(freq.data)) {
    if (verbose)
      print("Checking Frequency from the data --> [ok]")
  } else if (verbose) print("Checking Frequency from the data --> [fail]")
  ## Checking Calendar attribute if any, otherwise set to "ordinary"  # Possible values for CMIP5 files are : "365_day" , "standard" , "proleptic_gregorian" , "360_day"
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
    yorigin <- as.numeric(format.Date(torigin,format="%Y"))
    morigin <- as.numeric(format.Date(torigin,format="%m"))
    dorigin <- as.numeric(format.Date(torigin,format="%d"))
    horigin <- as.numeric(format.Date(torigin,format="%H"))
  }
  
  ## Get calendar from attribute if any and create vector of dates vdate
  ## 'hou'=strptime(torig,format="%Y-%m-%d %H") + time*3600
  ##
  if (!is.null(calendar.att)) {
    if (grepl("gregorian",calendar.att) | grepl("standard",calendar.att)) {
      ## if (grepl("sec",tunit))
      ##     time$vdate <- as.Date((time$vals/(24*60*60)),
      ##                           origin=as.Date(torigin))
      ## if (grepl("min",tunit))
      ##     time$vdate <- as.Date((time$vals/(24*60)),
      ##                           origin=as.Date(torigin))
      ## if (grepl("hou",tunit))
      ##     time$vdate <- as.Date((time$vals/24),origin=as.Date(torigin))
      ## if (grepl("day",tunit)) {
      ##     time$vdate <- as.Date((time$vals),origin=as.Date(torigin))   
      ## }
      if (grepl("mon",tunit)) {
        if (sum(round(diff(time$vals)) > 1) < 1) {
          year1 <- time$vals[1]%/%12 + yorigin
          month1 <- morigin
          torigin1 <- paste(as.character(year1),month1,"01",sep="-")
          ##         time$vdate <- seq(as.Date(torigin1), by = "month",length.out=length(time$vals))
          
          ##     } else print("Warning : Monthly data are Mangeled")
        }
      } 
      
      time$vdate <- switch(tunit,'seconds'= strptime(torigin,format="%Y-%m-%d %H%M%S") + time$vals,
                           'minutes'= strptime(torigin,format="%Y-%m-%d %H%M%S") + time$vals*60,
                           'hours'= strptime(torigin,format="%Y-%m-%d %H:%M:%S") + time$vals*3600,
                           'days'= as.Date(torigin) + time$vals,
                           'months'= seq(as.Date(torigin1),length.out=length(time$vals),by='month'),
                           'years'= year(as.Date(torigin)) + time$vals)
      
    } else if (!is.na(strtoi(substr(calendar.att, 1, 3))) | grepl("noleap",calendar.att)) {
      if (verbose) print(paste(substr(calendar.att,1, 3), "-days' model year found in calendar attribute"))
      if (grepl("noleap",calendar.att))
        time$daysayear <- 365
      else
        time$daysayear <- as.numeric(substr(calendar.att, 1, 3))
      if (!is.null(time$daysayear)) if (verbose) print(paste("Creating time$daysayear attribute and setting attribute to ", time$daysayear, sep=" "))
    }
    if (!is.null(time$daysayear) & tunit!='hours') {
      if (time$daysayear==365) 
        mndays <- c(31,28,31,30,31,30,31,31,30,31,30,31) # Number of days in each month
      else if (time$daysayear==360)
        mndays <- rep(30,12) # Number of days in each month
      ##else if
      ##mndays <- c(29.5,29.5,30.5,30.5,30.5,30.5,31.0,30.5,30.5,30.5,30.5,31.0)
      
      if (!is.null(time$daysayear) & !is.null(mndays)) {
        year1 <- time$vals[1]%/%time$daysayear + yorigin
        month1 <- morigin
        
        if (sum(diff(time$vals)%/%time$daysayear) > 1 & (verbose))
          print("Warning : Jumps of years has been found in the time series ")
        qf <- c(qf,"jumps of years found in time series")
        if (time$vals[1]%%time$daysayear > 27) {
          year1 <- year1 + 1
          month1 <- month1 + 1
        } 
        if (month1>12) month1 <- month1 - 12 
        # construct vdate
        months <- ((time$vals%%time$daysayear)%/%round(mean(mndays))) + 1
        years <- time$vals%/%time$daysayear + yorigin
        #shifting mndays by month1 to start with different initial months than january (1)
        ## KMP 2016-11-08 this doesn't work. why add a 0?
        #mndays <- c(0,mndays[month1:length(mndays)-1],mndays[1:month1-1])
        if(month1>1) mndays <- c(mndays[month1:length(mndays)],mndays[1:(month1-1)])
        days <- time$vals%%time$daysayear - rep(cumsum(mndays),time$len/12)
        #HBE added seperate test for seasonal data May 2nd 2018
        if (freq.data!='season') {
          if ((sum(diff(months) > 1) > 1) | (sum(diff(years) > 1) > 1) | (sum(round(abs(diff(days)))>2)) > 1) {
            print("Warning : Jumps in data have been found !")
            print("Warning: Trust the first date and force a continuous vector of dates !")
            time$vdate <- seq(as.Date(paste(as.character(year1),month1,"01",sep="-")), by = "month",length.out=time$len)
            qf <- c(qf,"jumps in data found - continuous vector forced")
          } else time$vdate <- as.Date(paste(years,months,"01",sep="-")) #round (days)
        } else time$vdate <- as.Date(paste(years,months,"15",sep="-")) 
      }
      # HBE added fix for no_leap or 365_day hourly calander
    } else if ( !is.null(time$daysayear) & (tunit =='hours') ) if (time$daysayear==365) {
      rankp<- seq(as.POSIXct(torigin), as.POSIXct("2200-01-01 00:00:00"), by="hour")
      rankp <- rankp[(month(rankp)!=2 | day(rankp)!=29)]
      time$vdate <- rankp[floor(time$vals)+1]
    } else if (verbose) {
      print(time$vdate[1])
      print(paste("Starting date : ",time$vdate[1],"Ending date : ",time$vdate[length(time$vdate)], sep = " "))
    }
  } else {
    if (verbose) print("warnings : Automatic detection of the calendar")
    calendar.detect <- "auto"
    ##                                     # NOT COMPLETE ...
    ##======= KMP 2016-05-04: as.Date does not work for sec, min, hour ========== 
    #if (grepl("sec",tunit)) time$vdate <- as.Date((time$vals/(24*60*60)),origin=as.Date(torigin))
    #if (grepl("min",tunit)) time$vdate <- as.Date((time$vals/(24*60)),origin=as.Date(torigin))
    #if (grepl("hou",tunit)) time$vdate <- as.Date((time$vals/24),origin=as.Date(torigin))
    if (grepl("sec",tunit)) time$vdate <- as.POSIXct(torigin,tz='UTC') + time$vals
    if (grepl("min",tunit)) time$vdate <- as.POSIXct(torigin,tz='UTC') + time$vals*60
    if (grepl("hou",tunit)) time$vdate <- as.POSIXct(torigin,tz='UTC') + time$vals*60*60
    ##===========================================================================
    if (grepl("day",tunit)) time$vdate <- as.Date((time$vals),origin=as.Date(torigin))   
    if (grepl("mon",tunit)) {
      if (sum(diff(time$vals>1)) < 1) {
        year1 <- time$vals[1]%/%12 + yorigin
        month1 <- morigin
        time$vdate <- seq(as.Date(paste(as.character(year1),month1,"15",sep="-")), by = "month",length.out=length(time$vals))
      } else print("Warning : Monthly data are mangeled") 
    } 
  }
  #HBE 11/04/18 added check for monthly time-unit in test below
  if ((length(time$vdate)>0) & (grepl("mon",tunit)) & (sum(diff(as.numeric(format.Date(time$vdate,"%m")))>1)) & (verbose)) stop("Vector date is mangeled ! Need extra check !")
  ## Checking the data / Extra checks / Automatic calendar detection / etc.
  ## Check 1 # Regular frequency
  ## 
  if (!is.null(time$vdate))
    dt <- as.numeric(rownames(table(diff(time$vdate))))
  else
    dt <- NULL
  if (!is.null(time$vdate)) {
    if (verbose) print("Vector of date is in the form :")
    if (verbose) print(str(time$vdate))
    if (verbose) print(diff(time$vdate))
  } else {
    if (grepl("sec",tunit))
      dt <- as.numeric(rownames(table(diff(ncid$dim$time$vals/(24*60*60)))))
    if (grepl("min",tunit))
      dt <- as.numeric(rownames(table(diff(ncid$dim$time$vals/(24*60)))))
    if (grepl("day",tunit))
      dt <- as.numeric(rownames(table(diff(ncid$dim$time$vals))))
    if (grepl("hou",tunit))
      dt <- as.numeric(rownames(table(diff(ncid$dim$time$vals/24))))
    if (grepl("mon",tunit))
      dt <- as.numeric(rownames(table(diff(ncid$dim$time$vals))))
    if (length(dt)==1) {
      if (verbose) print("Regular frequency has been detected from the data")
    } else if (verbose) print("Irregular frequency has been detected from the data")
    if ((length(dt)==3) & grepl("day",tunit)) {
      if (verbose) print(paste("Calendar is likely to be a 365-",tunit," with: ",as.character(length(dt))," irregular frequencies",sep = ""))
      dt <- c(28,30,31)
      if (verbose) print(paste(as.character(dt),tunit,sep="-"))
    }
    if ((length(dt)==4) & grepl("day",tunit)) {
      ## 
      if (verbose) print(paste("Calendar is likely to be a Gregorian 365/366-",tunit," with: ",as.character(length(dt))," irregular frequencies : ",sep=""))
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
    } else if (!is.null(freq.data)) model$frequency <- freq.data
    else if (sum(is.element(tolower(substr(tunit,1,3)), # REB 2016-03-03
                            c('sec','hou','day','mon','yea')))>0) {
      ## If the frequency is not provided, try to derive it from the time
      ## information - Time difference
      warning(paste('Need to guess the frequency, based on reckognised time units',tunit))
      model$frequency <- 1
    } else {
      stop("Frequency could not be found, neither detected, the data might be corrupted !")
    }
    if (!is.null(model$frequency)) {
      if (verbose) print(paste("Frequency set to ",model$frequency,sep=""))
      if (model$frequency=="month") {
        yr <-year(time$vdate)
        mo <- month(time$vdate)
        dy <- "01"
        time$vdate <- as.Date(paste(yr,mo,dy,sep="-"))     
      }
      #HBE 11/04/18 added regular season time-stamp
      if (model$frequency=="season") {
        yr <-year(time$vdate)
        mo <- month(time$vdate)
        dy <- "15"
        time$vdate <- as.Date(paste(yr,mo,dy,sep="-"))     
      }
    }
  }
  
  
  
  ## End check 2
  if (verbose) print("Checking --> [Done!]")
  ## use zoo library to format the data
  model$qf <- qf
  
  ## Extra Checking
  ##y.test <- data.e <- ncvar_get(ncid, v1$name, start = c(1,1,1), count = c(1,1,-1))
  ##ac.gcm <- data.frame(y = y.test, x1 = as.vector(cos(2 * pi * tim/daysayear)), x2 = as.vector(sin(2 * pi * tim/daysayear)))
  result <- list(model=model,time=time)
  invisible(result)
  
}


check.ncdf <- function(ncid, param="auto",verbose = FALSE) { ## use.cdfcont = FALSE - AM 22-10-2013 not used any more ! 
  ##
  ## Checking : Number of variables and select only one from the netcdf file, get variable attributes in v1. The user should decide between the list of variables
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
    if (is.null(v1)) stop(paste("Variable ",param," could not be found !",sep=""))
  }
  ## Checking : Variable dimensions ...
  ndims <- eval(parse(text=paste("ncid$var[[",i,"]]$ndims",sep="")))
  dimnames <- rep(NA,ndims)
  if (ndims>0) {
    for (j in 1:ndims) dimnames[j] <- eval(parse(text=paste("ncid$var[[",i,"]]$dim[[",j,"]]$name",sep="")))
    if (verbose) print("Checking Dimensions --> [ok]")
    if (verbose) print(paste(as.character(ndims), " dimension(s) has(have) been found :"))
    if (verbose) print(dimnames)
  } else {
    stop("Checking Dimensions --> [fail]")
    if (verbose) print("The variable has no dimensions. The file may be corrupted!")  
  }
  dimnames <- tolower(dimnames)
  ## Get global attributes in model, check and update
  model <- att.get.ncdf(ncid,0,"model")
  globalatt <- list('CDI','history','institution','Conventions','references',
                    'institute_id','experiment_id','model_id','forcing',
                    'parent_experiment_id','parent_experiment_rip','branch_time',
                    'contact','initialization_method','physics_version',
                    'tracking_id','version_number','product', 'experiment',
                    'frequency','creation_date','project_id','table_id','title',
                    'parent_experiment','modeling_realm','realization',
                    'cmor_version','CDO','NCO')
  get.att.val <- function(x,nc) return(att.get.ncdf(nc,0,x)$value)
  model <- lapply(as.matrix(globalatt),get.att.val,nc=ncid)
  names(model) <- globalatt
  ## if (verbose) print(model)
  ## Update CMIP3 attributes to match those of CMIP5 
  mnames <- names(model)
  history <- model$history ##att.get.ncdf(ncid,0,"history")
  ##
  ## Get project id 
  project.id <- model$project_id ##att.get.ncdf(ncid,0,"project_id")
  if ((is.null(project.id) | project.id==0) & (grepl(project.id,"CMIP3"))) {
    ##if (!is.null(model)) model$project_id <- project.id$value
    if (project.id=="IPCC Fourth Assessment") {
      model$project_id <- "CMIP3"
      if (verbose)
        print("project_id : IPCC Fourth Assessment set to CMIP3")
    }
  } else if (length(grep("sres",tolower(history)))>0)
    model$project_id <- "CMIP3"
  else if (length(grep("rcp",tolower(history)))>0)
    model$project_id <- "CMIP5"
  
  if (verbose)
    print(paste("project_id : ",model$project_id))
  
  ## Get model ID
  model.id <- model$model_id ##att.get.ncdf(ncid,0,"model_id")
  if ((is.null(model.id) | model.id==0) & !is.null(project.id)) {   
    ##hist2 <- unlist(strsplit(history$value,split=c(" ")))
    hist2 <- unlist(strsplit(history,split=c(" ")))
    ih <- grep("tas",hist2)
    if (model$project_id=="CMIP3") {
      txt <- hist2[ih[1]] 
      model$model_id <- cmip3.model_id(txt)
    }
    else if (model$project_id=="CMIP5") {
      txt <- hist2[ih[2]]
      model$model_id <- cmip5.model_id(txt)
    }
  }
  if (verbose)
    print(paste("Model_id : ",model$model_id))
  ## print(ncatt_get(ncid,0,"project_id")$hasatt)
  
  ## Get Experiment ID
  ##experiment.id <- att.get.ncdf(ncid,0,"experiment_id")
  experiment.id <- model$experiment_id
  if ((is.null(experiment.id) | experiment.id==0) & !is.null(model$project_id)) {
    model$experiment_id <- experiment.id
    if (tolower(model$project_id)=="cmip3") {
      txt <- unlist(strsplit(tolower(history),split="/"))
      model$experiment_id <- paste(txt[grep("20c3m",txt)],txt[grep("sres",txt)],sep="-")
    }
  }
  if (verbose)
    print(paste("experiment_id : ",model$experiment_id))
  ## 
  
  ## Get title 
  ##experiment.title <- model$title ## att.get.ncdf(ncid,0,"title")$hasatt
  ## 
  if (att.get.ncdf(ncid,0,"title")$hasatt) {
    exp.title <- att.get.ncdf(ncid,0,"title")$value
    modelid <- unlist(strsplit(exp.title,split=c(" ")))
    if (length(grep('-',tolower(modelid))>0) & 
        (is.null(model$model_id) | model$model_id==0))
      model$model_id <- modelid[grep('-',modelid)]
    if (length(grep('cmip',tolower(modelid))>0) &
        (is.null(model$project_id) | model$project_id==0))
      model$project_id <- modelid[grep('cmip',tolower(modelid))]
    if ( length(grep('rcp',tolower(modelid))>0) &
         (is.null(model$experiment_id)| model$experiment_id==0))
      model$experiment_id <- modelid[grep('rcp',tolower(modelid))]
    model$type <- modelid[grep('-',modelid)+1]
  }
  if (verbose)
    print(paste("experiment_title : ",exp.title))
  ## END CMIP3 MODEL NAMES UPDATE
  ##
  ## Checking : Time unit and origin
  ## Get system info
  a <- Sys.info()
  ## Get time dimension / val + attributes 
  itime <- grep("tim", dimnames)
  if (length(itime) == 0)
    itime <- NULL
  else if (length(itime) > 1) stop("Error in time dim")
  else if (length(itime)==1) time <- eval(parse(text=paste("v1$dim[[",as.character(itime),"]]",sep="")))
  
  ## Get time unit and origin
  tatt <- tolower(names(time))
  itunit <- grep(c("unit"),tatt)
  itorigin <- grep(c("orig"),tatt)
  if (length(itunit)>0) {   
    tunit <- eval(parse(text = paste("time$",tatt[itunit],sep="")))
    if (verbose)
      print(paste("Time unit has been found in time$unit attribute (",tunit,")",sep=""))
  } else tunit <- NULL
  if (!is.null(tunit)) {
    if (verbose) print("Checking Time Unit --> [ok]")}
  else if (verbose) print("Checking Time Unit --> [fail]")
  if (!is.null(tunit) & (!is.null(grep("since",tunit)))) {
    if (verbose) print("Time unit and origin detected in time$unit attribute")
    tunit <- time$units
    tsplit <- unlist(strsplit(tunit,split=" "))
    torigin <- time$origin <- paste(tsplit[3:length(tsplit)],collapse=" ")
    tunit <- time$units <- unlist(strsplit(tunit,split=" "))[1]
    if (verbose) print(paste("Updating time$unit (",time$unit,") and creating time$origin (",time$origin,") attribute",sep= ""))
  } else if (length(itorigin)>0) {   
    torigin <- eval(parse(text = paste("time$",tatt[itorigin],sep="")))
    if (verbose) print(paste("Time origin has been found in time origin attribute and set to:",torigin,sep=" "))
  } else torigin <- NULL
  ##if (is.null(torigin) & is.null(tunit)) {
  ##   if (verbose) print("Attributes time unit and origin have not been found -> Sys.time() will be used !")
  ##   if ((tolower(a[1]) == "linux") & (use.cdfcont)) {
  ##      if (verbose) print("Linux & use.cdfcont : call cdfcont()")
  ##      time$origin <- torigin <- cdfcont(ncid$filename)$time.origin
  ##      time$units <- tunit <- cdfcont(ncid$filename)$time.unit
  ##   }
  #}
  if (is.null(torigin)) {
    if (verbose) print(paste("Time units:", tunit, " l=", min(tim[is.finite(time$vals)]),"-", max(tim[is.finite(time$vals)])))
    if (verbose) warning("Cannot determine the time origin!")
    if (verbose) warning("Example format: '15-Dec-1949'")
    if (verbose) warning("NCEP reanalysis typically: 01-01-01")
    if (verbose) warning("ERA-40 typically: 1900-01-01")
    torigin <- readline("Please enter a valid time origin: ")
  }
  
  if (!is.null(torigin)) {
    if (torigin == "1-01-01 00:00:00") {
      if (verbose) print("bad time origin")
      torigin <- "0001-01-01 00:00:00"
      if (verbose) print(paste("Re-setting time origin (",torigin,")",sep=""))
    }
  } else {
    torigin <- readline("Give me the time origin (format='YYYY-MM-DD' as '1949-22-01'):")
    if (!is.null(torigin)) {
      if (verbose) print(paste("Time origin set to =", torigin))
      else stop("Could not determine the time origin. The processing has been stopped !")
    }
  }
  
  if (!is.null(torigin)) {
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
    if (is.na(dorigin)) {
      if (verbose) warning("Warning : Month origin is missing !")
      morigin <- 1
      if (verbose) warning("Warning : Month origin has been set to:",morigin)
    }
    torigin1 <- paste(yorigin,morigin,dorigin,sep="-")
    torigin <- paste(torigin1,unlist(strsplit(torigin,split=" "))[2],sep=" ") 
  }
  
  if (!is.null(torigin)) {
    if (verbose)
      print("Checking Time Origin --> [ok]")}
  else if (verbose)
    print("Checking Time Origin --> [fail]")
  ##
  ## Checking : Frequency
  type <- c("annual","season","month","day","hour","minute","second")
  type.abb <- substr(tolower(type),1,3)
  ## Initialize
  freq.att <- NULL
  ifreq <- grep("freq",names(model))
  if (length(ifreq)>0) {  
    itype <- grep(tolower(eval(parse(text=paste("model$",names(model)[ifreq],sep="")))),tolower(type))
    if (length(itype>0)) {
      if (verbose) print(paste("Frequency has been found in model$frequency attribute (",type[itype],")",sep="")) 
      freq.att <- frequency.name[grep(model$frequency,frequency.abb)]
    }
    if (verbose) print("Checking Frequency from attribute --> [ok]")
  } else {
    if (verbose) print("Checking Frequency from attribute --> [fail]")
    if (verbose) print("Frequency has not been found in the attributes") 
  }
  ##
  ## Checking frequency from data
  frequency <- freq.data <- NULL
  freq.data <- frequency.data(data=as.vector(time$vals),unit=tunit,verbose=FALSE)
  if (!is.null(freq.data)) {if (verbose) print("Checking Frequency from the data --> [ok]")} else if (verbose) print("Checking Frequency from the data --> [fail]")
  ## Checking Calendar attribute if any, otherwise set to "ordinary"  # Possible values for CMIP5 files are : "365_day" , "standard" , "proleptic_gregorian" , "360_day"
  ical <- grep(c("calend"),tatt)
  ## 
  if (length(ical)>0) {   
    calendar.att <- eval(parse(text = paste("time$",tatt[ical],sep="")))
    if (verbose) print("Checking Calendar from time attribute --> [ok]") 
    if (verbose) print(paste("Calendar attribute has been found in time$calendar (",time$calendar,")",sep =""))
  } else {
    if (verbose) print("Checking Calendar from time attribute --> [fail]")
    calendar.att <- NULL
    print("Warning : Calendar attribute has not been found in the meta data")
  }
  ## Identifying starting and ending dates for the data if possible
  if (!is.null(torigin)) {
    yorigin <- as.numeric(format.Date(torigin,format="%Y"))
    morigin <- as.numeric(format.Date(torigin,format="%m"))
    dorigin <- as.numeric(format.Date(torigin,format="%d"))
  }
  ## Get calendar from attribute if any and create vector of dates vdate
  
  if (!is.null(calendar.att)) {
    if (grepl("gregorian",calendar.att) | grepl("standard",calendar.att)) {
      if (grepl("sec",tunit)) time$vdate <- as.Date((time$vals/(24*60*60)),origin=as.Date(torigin))
      if (grepl("min",tunit)) time$vdate <- as.Date((time$vals/(24*60)),origin=as.Date(torigin))
      if (grepl("hou",tunit)) time$vdate <- as.Date((time$vals/24),origin=as.Date(torigin))
      if (grepl("day",tunit)) {
        time$vdate <- as.Date((time$vals),origin=as.Date(torigin))   
      }
      if (grepl("mon",tunit)) {
        if (sum(round(diff(time$vals)) > 1) < 1) {
          year1 <- time$vals[1]%/%12 + yorigin
          month1 <- morigin
          time$vdate <- seq(as.Date(paste(as.character(year1),month1,"01",sep="-")), by = "month",length.out=length(time$vals))
        } else print("Warning : Monthly data are Mangeled") 
      } 
    } else if (!is.na(strtoi(substr(calendar.att, 1, 3))) | grepl("noleap",calendar.att)) {
      if (verbose) print(paste(substr(calendar.att,1, 3), "-days' model year found in calendar attribute"))
      if (grepl("noleap",calendar.att))
        time$daysayear <- 365
      else
        time$daysayear <- as.numeric(substr(calendar.att, 1, 3))
      if (!is.null(time$daysayear)) if (verbose) print(paste("Creating time$daysayear attribute and setting attribute to ", time$daysayear, sep=" "))
    }
    if (!is.null(time$daysayear)) {
      if (time$daysayear==365) 
        mndays <- c(31,28,31,30,31,30,31,31,30,31,30,31) # Number of days in each month
      else if (time$daysayear==360)
        mndays <- rep(30,12) # Number of days in each month
      ##else if
      ##mndays <- c(29.5,29.5,30.5,30.5,30.5,30.5,31.0,30.5,30.5,30.5,30.5,31.0)
      if (!is.null(time$daysayear) & !is.null(mndays)) {
        year1 <- time$vals[1]%/%time$daysayear + yorigin
        month1 <- morigin
        
        if (sum(diff(time$vals%/%time$daysayear) > 1) & (verbose)) print("Warning : Jumps of years has been found in the time series ")
        if (time$vals[1]%%time$daysayear > 27) {
          year1 <- year1 + 1
          month1 <- month1 + 1
        } 
        if (month1>12) month1 <- month1 - 12 
        # construct vdate
        months <- ((time$vals%%time$daysayear)%/%round(mean(mndays))) + 1
        years <- time$vals%/%time$daysayear + yorigin
        #shifting mndays by month1 to start with different initial months than january (1)
        ## KMP 2016-11-08 doesn't work
        #mndays <- c(0,mndays[month1:length(mndays)-1],
        #            mndays[1:month1-1])
        if(month1>1) mndays <- c(mndays[month1:length(mndays)],mndays[1:(month1-1)])
        days <- time$vals%%time$daysayear - rep(cumsum(mndays),
                                                time$len/12)
        if ((sum(diff(months) > 1) > 1) | (sum(diff(years) > 1) > 1) | (sum(round(abs(diff(days)))>2)) > 1) {
          print("Warning: Jumps in data have been found !")
          print("Warning: Trust the first date and force a continuous vector of dates !")
          time$vdate <- seq(as.Date(paste(as.character(year1),
                                          month1,"01",sep="-")),
                            by = "month",length.out=time$len)
          qf <- c(qf,"jumps in data found - continuous vector forced")
        } else
          time$vdate <- as.Date(paste(years,months,"01",sep="-")) #round (days)                  
      }  
    } else  
      
      
      if (verbose) {
        print(time$vdate[1])
        print(paste("Starting date : ",time$vdate[1],"Ending date : ",
                    time$vdate[length(time$vdate)], sep = " "))
      }
  } else {
    if (verbose) print("warnings : Automatic detection of the calendar")
    calendar.detect <- "auto"
    ##                                     # NOT COMPLETE ...
    ##======= KMP 2016-05-04: as.Date does not work for sec, min, hour ========== 
    #if (grepl("sec",tunit)) time$vdate <- as.Date((time$vals/(24*60*60)),origin=as.Date(torigin))
    #if (grepl("min",tunit)) time$vdate <- as.Date((time$vals/(24*60)),origin=as.Date(torigin))
    #if (grepl("hou",tunit)) time$vdate <- as.Date((time$vals/24),origin=as.Date(torigin))
    if (grepl("sec",tunit)) time$vdate <- as.POSIXct(torigin,tz='UTC') + time$vals
    if (grepl("min",tunit)) time$vdate <- as.POSIXct(torigin,tz='UTC') + time$vals*60
    if (grepl("hou",tunit)) time$vdate <- as.POSIXct(torigin,tz='UTC') + time$vals*60*60
    ##===========================================================================
    if (grepl("day",tunit)) time$vdate <- as.Date((time$vals),origin=as.Date(torigin))   
    if (grepl("mon",tunit)) {
      if (sum(diff(time$vals>1)) < 1) {
        year1 <- time$vals[1]%/%12 + yorigin
        month1 <- morigin
        time$vdate <- seq(as.Date(paste(as.character(year1),
                                        month1,"15",sep="-")),
                          by = "month",length.out=length(time$vals))
      } else print("Warning : Monthly data are Mangeled") 
    } 
  }
  
  ## Checking the data / Extra checks / Automatic calendar detection / etc.
  ## Check 1 # Frequency
  ## 
  if (!is.null(time$vdate))
    dt <- as.numeric(rownames(table(diff(time$vdate))))
  else
    dt <- NULL
  
  if (length(dt)==1) {
    if (verbose)
      print("Regular frequency has been detected from the data")
  } else
    if (verbose)
      print("Irregular frequency has been detected from the data")
  
  if (!is.null(time$vdate)) {
    if (verbose) print("Vector of date is in the form :")
    if (verbose) str(time$vdate)
    if (verbose)
      print(paste("Time difference dt : ",paste(dt,collapse="/"))) ## diff(time$vdate))
  }
  ## 
  if (!is.null(time$vdate)) {
    if (grepl("sec",tunit))
      dt <- as.numeric(rownames(
        table(diff(ncid$dim$time$vals/(24*60*60)))))
    if (grepl("min",tunit))
      dt <- as.numeric(rownames(table(diff(ncid$dim$time$vals/(24*60)))))
    if (grepl("day",tunit))
      dt <- as.numeric(rownames(table(diff(ncid$dim$time$vals))))
    if (grepl("hou",tunit))
      dt <- as.numeric(rownames(table(diff(ncid$dim$time$vals/24))))
    if (grepl("mon",tunit))
      dt <- as.numeric(rownames(table(diff(ncid$dim$time$vals))))
    
    if ((length(dt)==3)) {## & grepl("day",tunit)) {
      if (verbose)
        print(paste("Calendar is likely to be a 365-",
                    tunit," with: ",as.character(length(dt)),
                    " irregular frequencies",sep = ""))
      calendar.att <- "365-days"
      dt <- c(28,30,31)
      if (verbose)
        print(paste(as.character(dt),tunit,sep="-"))
    }
    if (length(dt)==4) {## & grepl("day",tunit)) {
      if (verbose)
        print(paste("Calendar is likely to be a Gregorian 365/366-",
                    tunit," with: ",as.character(length(dt)),
                    " irregular frequencies : ",sep=""))
      calendar.att <- "Gregorian"
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
      if (verbose)
        print(paste("Warning : Irregular frequencies have been detected - The data might be corrupted and needs extra Checking !"))   
      if (verbose)
        print(paste(as.character(dt),tunit,sep=" "))
    }
    if (sum((dt >= 28) & (dt<=31))>0) freq.data <- "month"
    if (sum((dt >= 350) & (dt<=31))>0) freq.data <- "year"
  }
  
  ## set calendar automatically
  model$calendar <- calendar.att
  if (verbose) print(paste("Calendar set to : ",model$calendar))
  ##
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
        }
      } 
    } else if (!is.null(freq.data))
      model$frequency <- freq.data
    else
      stop("Frequency could not be found, neither detected, the data might be corrupted !")
    
    if (!is.null(model$frequency)) {
      if (verbose)
        print(paste("Frequency set to ",model$frequency,sep=""))
      if (model$frequency=="month") {
        yr <-year(time$vdate)
        mo <- month(time$vdate)
        dy <- "01"
        time$vdate <- as.Date(paste(yr,mo,dy,sep="-"))     
      }
    }
  }
  
  ## Mangeled months
  dm <- as.numeric(diff(as.numeric(format.Date(time$vdate,"%m"))))
  sdm <- sum(dm!=1 & dm!=-11)
  if ((length(time$vdate)>0) & (sum(dm!=1 & dm!=-11)>0)) {
    print(paste("converted months are mangeled (",sdm," jumps are found in the data)!",sep=""))
    year1 <- year(time$vdate)[1]
    month1 <- month(time$vdate)[1]
    time$vdate <- seq(as.Date(as.character(paste(year1,month1,"01",sep="-"))), by = freq.data,length.out=length(time$vals))
    print("Trusting first date and frequency to generate a new sequence of dates")
    qf <- c(qf,"jumps in data found - continuous vector forced")
  }
  if (median(as.numeric(row.names(table(diff(ncid$dim$time$vals))))) > 100
      & grepl('day',tunit)) {
    if (verbose) print("Looks like the data contains annual values")
    freq.data <- 'annual'
  }
  
  ## End check 2
  if (verbose) print("Checking --> [Done!]")
  ## use zoo library to format the data
  
  ## Extra Checking
  ##y.test <- data.e <- ncvar_get(ncid, v1$name, start = c(1,1,1), count = c(1,1,-1))
  ##ac.gcm <- data.frame(y = y.test, x1 = as.vector(cos(2 * pi * tim/daysayear)), x2 = as.vector(sin(2 * pi * tim/daysayear)))
  result <- list(model=model,time=time)
  invisible(result)
}

retrieve.station <- function(ncfile,param="auto",type="ncdf4",
                             path=NULL,is=NULL,stid=NULL,loc=NULL,lon=NULL,lat=NULL,it=NULL,
                             alt=NULL,cntr=NULL,start.year.before=NULL,end.year.after=NULL,
                             nmin=NULL,verbose=FALSE,onebyone=FALSE,...) {
  if (verbose) print(paste('retrieve.station',ncfile))
  
  ## REB 2018-04-06: Add a check for e.g. station data
  class.x <- file.class(ncfile)
  if (verbose) {print('Check class'); print(class.x$value)}
  stopifnot(tolower(class.x$value[1])=='station' | length(is.element(class.x$dimnames,'stid')) > 0)
  ncid <- nc_open(ncfile)
  if (param=='auto') param <- names(ncid$var)[1]
  size <- eval(parse(text=paste('ncid$var$',param,'$size',sep='')))
  if (verbose) {
    print(paste('The variable to read is',param))
    print(paste('Variable size in netCDF file:',paste(size,collapse=' - ')))
  }
  ## Read the metadata:
  tim <- ncvar_get(ncid,'time'); nt <- length(tim)
  stids <- ncvar_get(ncid,'stationID'); ns <- length(stids)
  if (verbose) {print('Get metadata');print(stids)}
  tunit <- ncatt_get(ncid,'time','units')
  lons <- ncvar_get(ncid,'lon')
  lats <- ncvar_get(ncid,'lat')
  alts <- ncvar_get(ncid,'alt')
  cntrs <- try(ncvar_get(ncid,'cntr'))
  nv <- try(ncvar_get(ncid,'number'))
  fyr <- try(ncvar_get(ncid,'first'))
  lyr <- try(ncvar_get(ncid,'last'))
  longname <- ncatt_get(ncid,param,'long_name')
  unit <- ncatt_get(ncid,param,'units')
  locs <- try(ncvar_get(ncid,'loc'))
  locs <- sub('\xc3','',locs)
  missing <- ncatt_get(ncid,param,'missing_value')
  ## Use the metadata to select the stations to read: there is no need to read
  ## all the stations if only a subset is desired
  if (verbose) print('Metadata is read - Select selected stations')
  if (is.null(is)) ii <- rep(TRUE,length(stids)) else
    if (!is.logical(is)) {
      ii <- rep(FALSE,length(stids)); ii[is] <- TRUE
    } else ii <- is
  if (!is.null(stid)) ii <- ii & is.element(stids,stid)
  if (!is.null(lon)) ii <- ii & (lons >= min(lon)) & (lons <= max(lon))
  if (!is.null(lat)) ii <- ii & (lats >= min(lat)) & (lats <= max(lat))
  if (!is.null(alt)) { 
    if (length(alt)==2) ii <- ii & (alts >= min(alt)) & (alts >= max(alt)) else
      if (alt > 0) ii <- ii & (alts >= alt) else 
        ii <- ii & (alts <= abs(alt))
  }
  #browser()
  if (!is.null(loc)) ii <- ii & 
    is.element(tolower(substr(locs,1,nchar(loc))),tolower(loc))
  if (!is.null(cntr)) ii <-ii & 
    is.element(tolower(substr(cntrs,1,nchar(cntr))),tolower(cntr))
  if (verbose) {print('Read following locations');
    print((1:ns)[ii]); print(locs[ii])}
  if (!is.null(nmin)) ii <- ii & (nv >= nmin)
  if (!is.null(start.year.before)) ii <- ii & (fy <= start.year.before)
  if (!is.null(end.year.after)) ii <- ii & (ly >= end.year.after)
  is <- (1:ns)[ii]
  
  if (onebyone) {
    if (verbose) print("Read stations one by one")
    ## For large files and where the stations are seperated far from each other in the
    ## netCDF file, it may be faster to read the files individually and then combine them 
    ## into one object
    X <- retrieve.station(ncfile,param=param,type=type,path=path,is=is[1],
                          it=it,start.year.before=start.year.before,
                          end.year.after=end.year.after,
                          nmin=nmin)
    if (length(is)>1) {
      for (ii in is[2:length(is)]) {
        x <- retrieve.station(ncfile,param=param,type=type,path=path,is=ii,
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
    print(paste(sum(tim > 10e7),'suspect time stamps!'))
    notsuspect <- tim <= 10e7
  } else notsuspect <- rep(TRUE,length(tim))
  if (verbose) {print('Time information'); print(tunit$value); print(range(tim))}
  if (length(grep('days since',tunit$value))) 
    t <- as.Date(substr(tunit$value,12,21)) + tim else
      if (length(grep('months since',tunit$value))) 
        t <- seq(as.Date(substr(tunit$value,14,23)),max(tim),'1 month') else
          if (length(grep('years since',tunit$value))) 
            t <- seq(as.Date(substr(tunit$value,14,23)),max(tim),'1 year')
  
  if (is.null(it)) {
    if (verbose) print('Read whole record')
    it1 <- 1; it2 <- nt
  } else {
    if (verbose) print('it is not NULL')
    if (is.character(it)) it <- as.Date(it)
    if (verbose) print(paste('Read selected period',min(it),'-',max(it),
                             'from interval',min(t),max(t)))
    it1 <- (1:length(t))[is.element(t,it)][1]
    it2 <- length(it)
    if (verbose) print(c(it1,it2))
    t <- t[it1:(it1+it2-1)]
  }
  
  start <- c(it1,min(is))
  count <- c(it2,max(is) - min(is)+1)
  
  ## Read the actual data:
  if (verbose) {
    print(paste('reading',param))
    print(paste('Number of stations in file=',ns,' reading',sum(ii)))
    print('Confirmation of station IDs:'); print(stids[is])
    print('start=')
    print(start)
    print('count=')
    print(count)
  }
  
  x <- ncvar_get(ncid,param,start=start,count=count)
  nc_close(ncid)
  if (verbose) print('All data has been extracted from the netCDF file')
  
  if (sum(!notsuspect)>0) {
    if (length(dim(x))==2) x <- x[notsuspect,] else x <- x[notsuspect] 
    t <- t[notsuspect]
    tim <- tim[notsuspect]; nt <- length(tim)
  }
  x[x<=missing$value] <- NA
  
  ## The data matrix is not full and may not necessarily correspond to the selection
  ## Need to remove unwanted stations with station numbers in the range of those selected
  iii <- seq(min((1:ns)[ii]),max((1:ns)[ii]),by=1)
  if (verbose) print(paste(' Number of stations read so far',length(iii),
                           ' Total number of stations',length(ii),
                           ' Selected stations=',sum(ii)))
  if (sum(ii)>1) {
    iv <- ii[iii]
    if (length(dim(x))==2) x <- x[,iv] else x <- x[iv]
  } else dim(x) <- NULL
  if (verbose) {
    print(paste('Dimensions of x is ',paste(dim(x),collapse=' - ')))
    print(summary(c(x))); print(sum(is.finite(x)))
  }
  if (length(dim(x))==2) { 
    nv <- apply(x,1,'nv') 
    jt <- (nv > 0)
  } else if (length(t)>1) jt <- is.finite(x) else jt <- is.finite(t)
  
  if (verbose) print(paste('Number of valid data points',length(jt)))
  lons <- lons[ii]; lats <- lats[ii]; alts <- alts[ii]; cntrs <- cntrs[ii]
  locs <- locs[ii]; stids <- stids[ii]
  if (verbose) print(paste('length(t)=',length(t),'length(x)=',length(x),'sum(jt)=',sum(jt)))
  if (length(dim(x))==2) x <- x[jt,] else if (length(t)>1) x <- x[jt]
  tim <- tim[jt]; t <- t[jt]
  
  if (length(t)==1) dim(x) <- c(1,length(x)) else 
    if (is.null(dim(x))) dim(x) <- c(length(x),1)
  if (verbose) print(paste('Dimensions of x is ',paste(dim(x),collapse=' - '),
                           'and length(t) is',length(t)))
  if (length(t) != dim(x)[1]) {
    print(paste('Dimensions of x is ',paste(dim(x),collapse=' - '),
                'and length(t) is',length(t)))
    stop('retrieve.station error:')
  }
  y <- as.station(zoo(x,order.by=t),loc=locs,lon=lons,lat=lats,alt=alts,
                  cntr = cntrs,stid = stids,longname=longname,
                  unit=unit$value,param=param)
  
  ## Weed out empty spaces
  if (length(t)>1) {
    if (verbose) print('Exclude empty time periods')
    if (length(dim(y))==2) iv <- apply(coredata(y),1,FUN='nv') else iv <- nv(y)
    y <- subset(y,it=iv > 0)
  }
  if (verbose) print(paste('as.station',min(index(y)),max(index(y)),'Data size=',length(x),'record length=',length(index(y))))
  if (verbose) print('exit retrieve.station')
  return(y)
}

retrieve.stationsummary <- function(ncfile,type="ncdf4",
                                    path=NULL,stid=NULL,loc=NULL,lon=NULL,lat=NULL,
                                    alt=NULL,cntr=NULL,start.year.before=NULL,end.year.after=NULL,
                                    nmin=NULL,verbose=FALSE,...) {
  if (verbose) print(paste('retrieve.stationsummary',ncfile))
  ## REB 2018-04-06: Add a check for e.g. station data
  class.x <- file.class(ncfile)
  if (verbose) {print('Check class'); print(class.x$value)}
  stopifnot(tolower(class.x$value[1])=='station' | length(is.element(class.x$dimnames,'stid')) > 0)
  ncid <- nc_open(ncfile)
  param <- names(ncid$var)[1]
  stats <- names(ncid$var)[grep('summary',names(ncid$var))]
  if (verbose) print(paste('The summary statistics to read is',stats))
  ## Read the metadata:
  tim <- ncvar_get(ncid,'time'); nt <- length(tim)
  stids <- ncvar_get(ncid,'stationID'); ns <- length(stids)
  if (verbose) {print('Get metadata');print(stids)}
  tunit <- ncatt_get(ncid,'time','units')
  lons <- ncvar_get(ncid,'lon')
  lats <- ncvar_get(ncid,'lat')
  alts <- ncvar_get(ncid,'alt')
  cntrs <- try(ncvar_get(ncid,'cntr'))
  nv <- try(ncvar_get(ncid,'number'))
  fyr <- try(ncvar_get(ncid,'first'))
  lyr <- try(ncvar_get(ncid,'last'))
  longname <- ncatt_get(ncid,param,'long_name')
  unit <- ncatt_get(ncid,param,'units')
  locs <- try(ncvar_get(ncid,'loc'))
  ## Fix names: remove unacceptable coding
  locs <- tolower(sub('\xc3','',locs))
  
  ## Order alphabetically
  if (verbose) 'Sort alphabetically and in proper order with Scandinavian characters'
  locssrt <- tolower(locs); locssrt <- sub('å','zzzå',locssrt)
  locssrt <- sub('æ','zzæ',locssrt); locssrt <- sub('ø','zzzø',locssrt)
  locssrt <- sub('ä','zzä',locssrt); locssrt <- sub('ö','zzzö',locssrt)
  srt <- order(locssrt); rm('locssrt')
  locs <- paste(toupper(substr(locs,1,1)),tolower(substr(locs,2,nchar(locs))),sep='')
  
  missing <- ncatt_get(ncid,param,'missing_value')
  y <- data.frame(location=locs,longitude=lons,latitude=lats,altitude=alts,country=cntrs,
                  number.valid=nv,first.year=fyr,last.year=lyr,station.id=stids)
  attr(y,'variable') <- param
  if (verbose) print('Created data.frame with metadata: y')
  if (length(grep('days since',tunit$value))) 
    t <- as.Date(substr(tunit$value,12,21)) + tim else
      if (length(grep('months since',tunit$value))) 
        t <- seq(as.Date(substr(tunit$value,14,23)),max(tim),'1 month') else
          if (length(grep('years since',tunit$value))) 
            t <- seq(as.Date(substr(tunit$value,14,23)),max(tim),'1 year')
  if (verbose) print(paste('Get the time period',paste(range(t),collapse=' - ')))
  ok <- ( (t > as.Date('1700-01-01')) & (t < as.Date('2300-01-01')) )
  attr(y,'period') <- range(t[ok])
  attr(y,'unit') <- unit
  attr(y,'missing_value') <- missing
  attr(y,'length') <- length(t)
  if (verbose) print('got metadata attributes: period, unit, missing')
  
  if(verbose) print('Read the summary statistics')
  for (i in 1:length(stats)) {
    if(verbose) print(stats[i])
    sname <- substr(stats[i],9,nchar(stats[i]))
    z <- try(ncvar_get(ncid,stats[i]))
    eval(parse(text=paste('y[["',sname,'"]] <- z',sep='')))
  }
  
  nc_close(ncid)
  y <- y[srt,]
  good <- (y$location != "")
  y <- y[good,]
  y$location <- as.character(y$location)
  if (verbose) print('Data is extracted from the netCDF file')
  if(verbose) print(summary(y))
  return(y)
}

## Used to check the contents in netCDF file - to use in retrieve to call retrieve.dsenemble,
## retrieve.eof or retrieve.station rather than the standard form to read field objects.
## Assumes that empty class attribute means a field object
file.class <- function(ncfile,type="ncdf4") {
  if (type=='ncdf4') {
    nc <- nc_open(ncfile)
    dimnames <- names(nc$dim)
    class.x <- ncatt_get(nc,0,'class')
    nc_close(nc)
  } else {
    nc <- open.ncdf(ncfile)
    dimnames <- names(nc$dim)
    class.x <- get.att.ncdf(nc,0,'class')
    close.ncdf(nc)
  }
  #attr(class.x,'dimnames') <- dimnames
  class.x$dimnames <- dimnames
  return(class.x)
}

