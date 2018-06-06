## Author=? Date 

## https://www.unidata.ucar.edu/software/netcdf/docs/netcdf/CDF-Data-Types.html:
## short: 16-bit signed integers. The short type holds values between -32768 and 32767.

## Help functions 
firstyear <- function(x) {
  yrs <- year(x)
  if (is.null(dim(x))) y <- min(yrs[is.finite(x)]) else
                       y <- apply(x,2,function(x,yrs=yrs) min(yrs[is.finite(x)]),yrs)
  return(y)
}

lastyear <- function(x) {
  yrs <- year(x)
  if (is.null(dim(x))) y <- max(yrs[is.finite(x)]) else
                       y <- apply(x,2,function(x,yrs=yrs) max(yrs[is.finite(x)]),yrs)
  return(y)
}

write2ncdf4 <- function(x,...) UseMethod("write2ncdf4")

write2ncdf4.default <- function(x,...) {
}

write2ncdf4.list <- function(x,fname='field.nc',prec='short',scale=0.1,offset=NULL,
                             torg="1970-01-01",missval=-999,verbose=FALSE) {
  if (verbose) print('write2ncdf4.list')
  stopifnot(inherits(x[[1]],'field'))
  ## Write.list is meant to add several fields to one netCDF file
  if (verbose) print(names(x))
  n <- length(x)
  ## Accomodate for the possibility with different precisions, scaling factor, etc
  ## If one is given, use it for all variables
  if (is.null(scale)) scale <- 1
  if (is.null(offset)) offset <- 0
  if (length(prec)==1) prec <- rep(prec,n)
  if (length(scale)==1) scale <- rep(scale,n)
  if (length(offset)==1) offset <- rep(offset,n)
  
  if (verbose) print(attr(x[[1]],'dimensions'))
  
  dimlon <- ncdim_def( "longitude", "degree_east", lon(x[[1]]) )
  dimlat <- ncdim_def( "latitude", "degree_north", lat(x[[1]]) )
  if (inherits(index(x),c('numeric','integer')))
    index(x[[1]]) <- as.Date(paste(index(x[[1]]),'-01-01',sep=''))
  
  dimtim <- ncdim_def( "time", paste("days since",torg),
                       as.numeric(as.Date(index(x[[1]]),origin=torg)) )
  varids <- unlist(lapply(x,function(x) varid(x)[1]))
  if (length(varids) != length(names(x))) varids <- names(x)
  units <- unlist(lapply(x,function(x) unit(x)[1]))
  if (verbose) {print(varids); print(units); print(n)}
  x4nc <- list()
  for (i in 1:n) {
    if (verbose) print(paste(i,'ncvar_def',varids[i]))
    x4nc[[varids[i]]] <- ncvar_def(varids[i], units[i], list(dimlon,dimlat,dimtim), -1, 
                                   longname=attr(x[[i]],'longname'), prec=prec[i])
  }
  
  # Create a netCDF file with this variable
  ncnew <- nc_create( fname, x4nc )
  
  # Write some values to this variable on disk.
  for (i in 1:n) {
    if (verbose) print(names(x)[i])
    y <- coredata(x[[i]])
    if (is.null(offset[i])) offset[i] <- mean(y,na.rm=TRUE)
    if (is.null(scale[i])) scale[i] <- 1
    y <- t(y)
    y[!is.finite(y)] <- missval
    y <- round((y-offset[i])/scale[i])
    if (verbose) {
      print(dim(y)); print(attr(y,'dimensions'))
      print(offset[i]); print(scale[i])
      }
    dim(y) <- attr(x,'dimensions')
    if (verbose) print(summary(round(y)))
    ncvar <- x4nc[[varids[i]]]
    ncvar_put( ncnew, ncvar, round(y) )
    ncatt_put( ncnew, ncvar, "add_offset", offset[i], prec="float" )
    ncatt_put( ncnew, ncvar, "scale_factor", scale[i], prec="float" ) 
    ncatt_put( ncnew, ncvar, "_FillValue", missval, prec="float" ) 
    ncatt_put( ncnew, ncvar, "missing_value", missval, prec="float" ) 
    history <- toString(attr(x[[i]],'history')$call)
    ncatt_put( ncnew, ncvar, "history", history, prec="text" ) 
  }
  ncatt_put( ncnew, 0, 'class', class(x))
  ncatt_put( ncnew, 0, "description", 
             paste("Saved from esd using write2ncdf4",date()))
  ncatt_put( ncnew, 0, "esd-version", attr(x[[1]],'history')$session$esd.version)
  
  nc_close(ncnew)
  if (verbose) print('netCDF file saved')
}

write2ncdf4.field <- function(x,fname='field.nc',prec='short',scale=NULL,offset=NULL,
                              torg="1970-01-01",missval=-999,ncclose=TRUE,verbose=FALSE) {
  if (verbose) {print('write2ncdf4.field'); print(names(attributes(x)))}

  y <- coredata(x)
  if (is.null(offset)) offset <- mean(y,na.rm=TRUE)
  if (is.null(scale)) scale <- (max(abs(c(y)),na.rm=TRUE) - offset)/10000
  y <- t(y)
  y <- round((y-offset)/scale)
  y[!is.finite(y)] <- missval
  if (verbose) {
    print(attr(y,'dimensions')); print(c(scale,offset))
    print(range(c(y))); print(range(c(x),na.rm=TRUE))
  }
  dim(y) <- attr(x,'dimensions')

  dimlon <- ncdim_def( "longitude", "degree_east", lon(x) )
  dimlat <- ncdim_def( "latitude", "degree_north", lat(x) )
  if (inherits(index(x),c('numeric','integer')))
      index(x) <- as.Date(paste(index(x),'-01-01',sep=''))
  
  dimtim <- ncdim_def( "time", paste("days since",torg),
                      as.numeric(as.Date(index(x),origin=torg)) )
  x4nc <- ncvar_def(varid(x)[1], unit(x)[1], list(dimlon,dimlat,dimtim), -1, 
                    longname=attr(x,'longname'), prec=prec)
     
     # Create a netCDF file with this variable
  ncnew <- nc_create( fname, x4nc )

  # Write some values to this variable on disk.
  ncvar_put( ncnew, x4nc, round(y) )
  ncatt_put( ncnew, x4nc, "add_offset", offset, prec="float" )
  ncatt_put( ncnew, x4nc, "scale_factor", scale, prec="float" ) 
  ncatt_put( ncnew, x4nc, "_FillValue", missval, prec="float" ) 
  ncatt_put( ncnew, x4nc, "missing_value", missval, prec="float" ) 
  history <- toString(attr(x,'history')$call)
  ncatt_put( ncnew, x4nc, "history", history, prec="text" ) 
  ncatt_put( ncnew, 0, 'class', class(x))
  ncatt_put( ncnew, 0, "description", 
             paste("Saved from esd using write2ncdf4",date()))
  if (verbose) print(attr(x,'history'))
  ncatt_put( ncnew, 0, "esd-version", attr(x,'history')$session$esd.version)
  nc_close(ncnew)
}



# https://www.unidata.ucar.edu/software/netcdf/docs/netcdf/CDL-Data-Types.html:
# short: 16-bit signed integers. The short type holds values between -32768 and 32767. 

write2ncdf4.station <- function(x,fname,prec='short',offset=0, missval=-999,tim=NULL,stano=NULL,append=FALSE,
                                scale=0.1,torg='1899-12-31',verbose=FALSE) {
  #require(ncdf4)

  if (!inherits(x,"station")) stop('x argument must be a station object') 
  
  if (verbose) print('write2ncdf4.station')
  
  ## Write a station object as a netCDF file using the short-type combined with add_offsetet and scale_factor
  ## to reduce the size.   
  ## Examine the station object: dimensions and attributes  

  ## Get time 
  if (is.null(tim)) nt <- dim(x)[1] else nt <- length(tim)
  if (!is.null(stano))  {
    if (!is.null(dim(x))) nstations <- dim(x)[2] else if (!is.null(x)) nstations <- 1 else nstations <- 0
    if (append & (length(stano) != nstations)) {
      print(dim(x)); print(length(stano))
      stop('write2ncdf4.station: stano argument does not match x')
    }
    ns <- length(stano) 
  } else {
    ns <- dim(x)[2] 
    stano <- 1:ns
  }
  if (verbose) print(c(nt,ns))
  
  ## if (is.null(d)) d <- c(length(x),1)
  if (verbose) print(paste('Number of stations: ',paste(ns)))
  
  atts <- names(attributes(x))
 
  if (verbose) print(atts)
  attr2attr <- is.element(atts,c('station_id','variable','unit','longname','location','country'))
  ##atts <- atts[iattr2ncdf]
  attr2var <- is.element(atts,c('longitude','latitude','altitude'))
  na <- length(atts); la <- rep(0,na)
  attrprec <- rep('character',na)
  attr(x,'quality') <- as.character(attr(x,'quality'))
  if (verbose) print(paste('attributes:', paste(atts, collapse=', '),
                           '; types:',paste(attrprec, collapse=', ')))
  
  #fyr <- firstyear(x)
  #lyr <- lastyear(x)
  #nv <- apply(coredata(x),2,'nv')
  
  y <- coredata(x)
  y[!is.finite(y)] <- missval
  y <- round((y - offset)/scale)
  
  if (is.null(attr(x,'calendar')))
      attr(x,'calendar') <- 'standard'

  if (!is.na(attr(x,'calendar')))
      calendar <- attr(x,'calendar')
  else calendar <- 'standard'

  if (class(index(x))=='Date')
      if (is.null(tim)) time <- julian(index(x)) - julian(as.Date(torg)) else {
        if (verbose) print('Use prescribed time coordinates')
        y <- zoo(y,order.by=index(x))
        x2 <- merge(zoo(rep(NA,nt),order.by=tim),zoo(y),all=FALSE)
        x2 <- window(x2[,-1],start=tim[1],end=tim[length(tim)])
        x2 <- attrcp(x,x2); class(x2) <- class(x); y <- x2; rm('x2')
        time <- julian(index(y)) - julian(as.Date(torg))
        if (verbose) {print(range(index(x))); print(dim(x))}
      }
  else if (inherits(x,'annual'))
      time <- julian(as.Date(paste(year(x),'01-01',sep='-')))-julian(as.Date(torg))
       
# Attributes with same number of elements as stations are saved as variables
  
  if (is.null(tim)) tim <- index(y)
  start <- c( (1:length(tim))[is.element(tim,index(y)[1])],stano[1] )
  if (!is.null(dim(y))) count <- dim(y) else count <- c(length(y),1)
  if (verbose) {
    print("start & count"); print(start); print(count); 
    print("dim(y)"); print(dim(y))
    print("netCDF dimensions"); print(c(nt,ns)); print(start+count-c(1,1))
    }
  
# Define the dimensions
  if (!append) {
    if (verbose) print('Define dimensions')
    if (verbose) print(stid(x))
    dimS <- ncdim_def( name="stid", units="number",vals=1:ns,unlim=TRUE)
    dimT <- ncdim_def( name="time", units=paste("days since",torg), vals=1:nt, calendar=calendar,unlim=TRUE)
    dimnchar   <- ncdim_def("nchar",   "", 1:12, create_dimvar=FALSE )
    #dimstation <- ncdim_def("station", "", 1:ns, create_dimvar=FALSE )
  
    if (verbose) {
      print('Define variable')
      print(paste('create netCDF-file',fname))
      print(summary(c(y)))
    }
  
    if (verbose) print('Define the netCDF structure')
    latid <- ncvar_def(name="lat",dim=list(dimS), units="degrees_north", missval=missval,longname="latitude", 
                       prec="float",verbose=verbose)

    lonid <- ncvar_def(name="lon",dim=list(dimS), units="degrees_east", missval=missval,longname="longitude", 
                       prec="float",verbose=verbose)
    altid <- ncvar_def(name="alt",dim=list(dimS), units="meters", missval=missval,longname="altitude", 
                       prec=prec,verbose=verbose)

    #locid <- ncvar_def(name="loc",dim=list(dimnchar,dimstation),units="NA",prec="char",longname="location",
    #                   verbose=verbose)
    #stid <- ncvar_def(name="stationID",dim=list(dimnchar,dimstation),units="NA",prec="char",longname="station_id",
    #                  verbose=verbose)
    #cntrid <- ncvar_def(name="cntr",dim=list(dimnchar,dimstation),units="NA",prec="char",longname="country",
    #                    verbose=verbose)
    locid <- ncvar_def(name="loc",dim=list(dimnchar,dimS),units="name",prec="char",longname="location",
                       verbose=verbose)
    stid <- ncvar_def(name="stationID",dim=list(dimnchar,dimS),units="number",prec="char",longname="station_id",
                      verbose=verbose)
    cntrid <- ncvar_def(name="cntr",dim=list(dimnchar,dimS),units="name",prec="char",longname="country",
                        verbose=verbose)
  
    fyrid <- ncvar_def(name="first",dim=list(dimS), units="year", missval=missval,longname="first_year", 
                       prec="short",verbose=verbose)
    lyrid <- ncvar_def(name="last",dim=list(dimS), units="year", missval=missval,longname="last_year", 
                       prec="short",verbose=verbose)
    nvid <- ncvar_def(name="number",dim=list(dimS), units="count", missval=missval,longname="number_valid_data", 
                       prec="float",verbose=verbose)
    
    ncvar <- ncvar_def(name=varid(x)[1],dim=list(dimT,dimS), units=ifelse(unit(x)[1]=="°C", "degC",unit(x)[1]),
                         longname=attr(x,'longname')[1], prec=prec,compression=9,verbose=verbose)
    if (verbose) print('The variables have been defined')
  } 
  
  if (append & file.exists(fname)) {
    if (verbose) print(paste('Appending',fname))
    ncid <- nc_open(fname, write=TRUE)
    ncvar <- ncid$var[[1]]
    lonid <- ncid$var[["lon"]]
    latid <- ncid$var[["lat"]]
    altid <- ncid$var[["alt"]]
    locid <- ncid$var[["loc"]]
    cntrid <- ncid$var[["cntr"]]
    fyrid <- ncid$var[["first"]]
    lyrid <- ncid$var[["last"]]
    nvid <- ncid$var[["number"]]
    stid <- ncid$var[["stationID"]]
    dimS <- ncid$dim[["stid"]]; ns <- dimS$len
    dimT <- ncid$dim[["time"]]; nt <- dimT$len
  } else {
    if (verbose) print(paste('Creating file',fname))
    ncid <- nc_create(fname,vars=list(ncvar,lonid,latid,altid,locid,stid,cntrid, 
                                      fyrid,lyrid,nvid)) ## vars)
  }
 
  if (verbose) print('Saving the variables:')
  ncvar_put( ncid, ncvar, coredata(y),start=start,count=count)
  ncatt_put( ncid, ncvar, 'add_offset',offset,prec='float')
  ncatt_put( ncid, ncvar, 'scale_factor',scale,prec='float')
  ncatt_put( ncid, ncvar, 'missing_value',missval,prec='float')
  ncvar_put( ncid, lonid, lon(y),start=start[2],count=count[2])
  ncvar_put( ncid, latid, lat(y),start=start[2],count=count[2])
  ncvar_put( ncid, altid, alt(y),start=start[2],count=count[2])

  ncvar_put( ncid, fyrid, firstyear(x),start=start[2],count=count[2])
  ncvar_put( ncid, lyrid, lastyear(x),start=start[2],count=count[2])
  if (is.null(dim(x))) number <- sum(is.finite(coredata(x))) else
  if (length(dim(x))==2) number <- apply(coredata(x),2,FUN='nv') else number <- -1
  ncvar_put( ncid, nvid, number,start=start[2],count=count[2])
  
  ## There are some times problems saving text data, and there seems to be some 
  ## inconsistency in the ncdf4 package. To by-pass this problem, we had to make
  ## the following code more complicated. There seems to be a mix-up between the 
  ## dimensions sometimes.
  if (verbose) print('Saving textual information')
  test <- try(ncvar_put( ncid, locid, loc(y),start=c(1,start[2]),count=c(12,count[2])))
  if (inherits(test,'try-error'))
    try(ncvar_put( ncid, locid, loc(y),start=c(start[2],1),count=c(count[2],12)))
  test <- try(ncvar_put( ncid, stid, as.character(stid(y)),c(1,start[2]),count=c(12,count[2])))
  if (inherits(test,'try-error'))
    try(ncvar_put( ncid, stid, as.character(stid(y)),c(start[2],1),count=c(count[2],12)))
  test <- try(ncvar_put( ncid, cntrid, cntr(y),start=c(1,start[2]),count=c(12,count[2])))
  if (inherits(test,'try-error'))
    try(ncvar_put( ncid, cntrid, cntr(y),start=c(start[2],1),count=c(count[2],12)))
  if (verbose) print('textual data saved')
  
  if (!append) {
  ## global attributes
    if (verbose) print('global attributes')
    ncatt_put( ncid, 0, 'class', class(x))
    ncatt_put( ncid, 0, 'title', paste(levels(factor(attr(x,"info"))),collapse="/"))
    ncatt_put( ncid, 0, 'source', paste(levels(factor(attr(x,"source"))),collapse="/"))
    ncatt_put( ncid, 0, 'history', paste(unlist(attr(x,"history")),collapse="/"))
    ncatt_put( ncid, 0, 'references', paste(levels(factor(attr(x,"reference"))),collapse="/"))
    ncatt_put( ncid, 0, "esd-version", attr(x,'history')$session$esd.version)
  }
  nc_close(ncid)
  if (verbose) print('close')
}


## These small functions are common code that simplify saving data as netCDF 
write2ncdf4.pca <- function(x,fname='esd.pca.nc',prec='short',verbose=FALSE,scale=0.01,offset=0,missval=-99) {
  if (verbose) print('write2ncdf4.pca')
  pcaatts <- names(attributes(x))
  pattern <- attr(x,'pattern')
  pattern[!is.finite(pattern)] <- missval
  dpat <- dim(pattern); attributes(pattern) <- NULL; dpat -> dim(pattern)
  if (verbose) print(pcaatts)
  if (class(index(x))=='Date') index(x) <- year(x) + (month(x)-1)/12
  ## set up the dimensions of the PCA
  dimpca <- ncdim_def( "i_pca", "index", 1:dim(x)[2] )
  dimtim <- ncdim_def( "time_pca", "year", index(x) )
  dimxy <- ncdim_def( "space_pca", "index", 1:dim(pattern)[1] )
  dimsea <- ncdim_def( "season", "index", 1:4 )
  pca <- ncvar_def("pca", "weights", list(dimtim,dimpca), missval, 
                    longname='principal components', prec=prec)
  pat <- ncvar_def("pattern_pca", "weights", list(dimxy,dimpca), missval, 
                    longname='principal component analysis patterns', prec=prec)
  lon <- ncvar_def("longitude", "degree_east", dimxy,missval, 
                    longname='longitude', prec='float')
  lat <- ncvar_def("latitude", "degree_north", dimxy,missval, 
                    longname='latitude', prec='float')
  alt <- ncvar_def("altitude", "m", dimxy,missval, 
                    longname='altitude', prec='float')
  stid <- ncvar_def("station_id", "number", dimxy,missval, 
                    longname='station ID', prec="integer")
#  loc <- ncvar_def("location", "name", dimxy,"NA", 
#                    longname='location name', prec="char")
#  cntr <- ncvar_def("country", "name", dimxy,"NA", 
#                    longname='country name', prec="char")
#  src <- ncvar_def("src", "name", dimxy,"NA", 
#                    longname='source', prec="char")
  lambda <- ncvar_def("lambda", "number", dimpca,missval, 
                    longname='eigenvalues', prec="float")
  if (is.numeric(index(x))) index(x) <- year(x)
  dpca <- dim(pca); attributes(pca) <- NULL; dpca -> dim(pca)
  ncvar_put( nc, pca, round((pca - offset)/scale) )
  ncatt_put( nc, pca, "add_offset", offset, prec="float" )
  ncatt_put( nc, pca, "scale_factor", scale, prec="float" ) 
  ncatt_put( nc, pca, "_FillValue", missval, prec="float" ) 
  ncatt_put( nc, pca, "missing_value", missval, prec="float" ) 
  history <- toString(attr(x,'history')$call)
  ncatt_put( nc, pca, "history", history, prec="text" ) 
  ncvar_put( nc, pat, round((pattern - offset)/scale) )
  ncatt_put( nc, pat, "add_offset", offset, prec="float" )
  ncatt_put( nc, pat, "scale_factor", scale, prec="float" ) 
  ncatt_put( nc, pat, "_FillValue", missval, prec="float" ) 
  ncatt_put( nc, pat, "missing_value", missval, prec="float" ) 
  ncatt_put( nc, pat, "dimensions_pca", paste(attr(pattern,'dimensions'),collapse=', '), prec="char" ) 
  ncatt_put( nc, pat, "locations", paste(loc(x),collapse=','), prec="char" ) 
  ncatt_put( nc, pat, "country", paste(cntr(x),collapse=','), prec="char" ) 
  ncatt_put( nc, pat, "source", paste(src(x),collapse=','), prec="char" ) 
  ncvar_put( nc, tim, index(x) )
  ncvar_put( nc, lon, lon(x) )
  ncvar_put( nc, lat, lat(x) )
  ncvar_put( nc, alt, alt(x) )
#  ncvar_put( nc, loc, loc(x) )
  ncvar_put( nc, stid, stid(x) )
#  ncvar_put( nc, cntr, cntr(x) )
#  ncvar_put( nc, src, src(x) )
  ncvar_put( nc, lambda, attr(x,'eigenvalues') )
  ncatt_put( nc, pca, "history", paste(attr(x,'history'),collapse=';'), prec="char" )
  ncatt_put( nc, 0, 'class', class(x))
  ncatt_put( nc, 0, "esd-version", attr(x,'history')$session$esd.version)
}

write2ncdf4.eof <- function(x,fname='eof.nc',prec='short',scale=10,offset=NULL,torg="1970-01-01",missval=-999) {
}

  
write2ncdf4.dsensemble <- function(x,fname='esd.dsensemble.nc',prec='short',offset=0,scale=0.1,
                              torg="1970-01-01",missval=-99,verbose=TRUE) {
  ## prec - see http://james.hiebert.name/blog/work/2015/04/18/NetCDF-Scale-Factors/
  if (verbose) print('write2ncdf4.field')
  class.x <- class(x)
  ngcms <- length(x) - 2
  ## Get the two first elements of the list which contain information common to the rest
  info <- x$info
  if (!is.null(x$pca)) pca <- x$pca
  if (!is.null(x$eof)) pca <- x$eof
  lons <- lon(pca)
  lats <- lat(pca)
  
  ## Clear these so that the rest are zoo objects describing the model runs
  x$info <- NULL
  x$pca <- NULL
  x$eof <- NULL
  names.x <- names(x)
  if (verbose) print(names.x)
  
  ## Define the variables and dimensions
  if (verbose) {print('Check index type - set to year'); print(index(x[[1]]))}
  if (is.list(x)) {
    if (class(index(x))=='Date') tim <- year(x[[1]]) + (month(x[[1]])-1)/12 else
                                 tim <- year(x[[1]])
    dimgcm <- dim(x[[1]])
    nloc <- length(x)
  } else {
    if (class(index(x))=='Date') tim <- year(x) + (month(x)-1)/12 else
                                 tim <- year(x)
    dimgcm <- dim(x)
    nloc <- 1
  }
  if (!is.null(pca)) {
    pcaatts <- names(attributes(pca))
    pattern <- attr(pca,'pattern')
    dpat <- dim(pattern)
  } else dpat <- 1
  if (verbose) {print(pcaatts); print(dim(pattern))}
  
  ## Set dimensions  
  dimtim <- ncdim_def( "time", 'year', as.integer(tim) )
  dimens <- ncdim_def( "ensemble_member", 'index', 1:nloc )
  if (!is.null(pca)) {
    if (class(index(pca))=='Date') index(pca) <- year(pca) + (month(pca)-1)/12 else
                                   index(pca) <- year(pca)
    dimpca <- ncdim_def( "i_pca", "index", 1:dimgcm[2] )
    dimtimpca <- ncdim_def( "time_pca", "year", index(pca) )
    if (length(dpat)==2) {
      dimxy <- ncdim_def( "space_pca", "index", 1:dim(pattern)[1] )
    } else {
      dimx <- ncdim_def( "longitude", "degrees_east", lons )
      dimy <- ncdim_def( "latitude", "degrees_north", lats )
    }
  }
  dimsea <- ncdim_def( "season", "index", 1:4 )
  if (verbose) print(tim)
  varlist <- list() # For single-station dsensemble objects
  varlist$gcm <- ncvar_def("gcm", "weights", list(dimtim,dimpca,dimens), missval, 
                           longname='principal components', prec=prec)

  ## set up the dimensions of the PCA
  if (!is.null(pca)) {
    varlist$pca <- ncvar_def("PC", "weights", list(dimtimpca,dimpca), missval, 
                             longname='principal components', prec=prec)
    ## EOFs and PCA have different number of dimensions
    if (length(dpat)==2) {
      if (verbose) print('--- PCA ---')
      varlist$pat <- ncvar_def("pattern", "weights", list(dimxy,dimpca), missval, 
                               longname='principal component analysis patterns', prec=prec)
      varlist$lon <- ncvar_def("longitude", "degree_east", dimxy,missval, 
                               longname='longitude', prec='float')
      varlist$lat <- ncvar_def("latitude", "degree_north", dimxy,missval, 
                               longname='latitude', prec='float')
      varlist$alt <- ncvar_def("altitude", "m", dimxy,missval, 
                               longname='altitude', prec='float')
      varlist$stid <- ncvar_def("station_id", "number", dimxy,"NA", 
                                longname='station ID', prec="char")
      varlist$loc <- ncvar_def("location", "name", dimxy,"NA", 
                               longname='location name', prec="char")
      varlist$cntr <- ncvar_def("country", "name", dimxy,"NA", 
                                longname='country name', prec="char")
      varlist$src <- ncvar_def("src", "name", dimxy,"NA", 
                               longname='source', prec="char")
    } else
    if (length(dpat)==3) {
      if (verbose) print('--- EOF ---')
      varlist$pat <- ncvar_def("pattern", "weights", list(dimx,dimy,dimpca), missval, 
                               longname='principal component analysis patterns', prec=prec)
          }
    varlist$lambda <- ncvar_def("lambda", "number", dimpca,missval, 
                                longname='eigenvalues', prec="float")
  }

  if (verbose) print(names(varlist))
     
  ## Create a netCDF file with this variable
  if (verbose) print('Create netCDF-file')

  ncnew <- nc_create( fname, varlist,verbose=verbose)
  
  if (verbose) print('write pca/eof data')                   
  ## Add the information stored in the list elements as 2D zoo objects
  X <- unlist(lapply(x,function(x) x[1:239,]))
  X <- round((X - offset)/scale)
  if (verbose) print(c(length(X),dim(x[[1]]),length(x)))
  if (verbose) print(table(unlist(lapply(x,dim))))
  dim(X) <- c(dim(x[[1]]),length(x))
  if (verbose) print('write zoo data')     
  ## Write some values to this variable on disk: GCM results
  ncvar_put( ncnew, varlist$gcm, X )
  ncatt_put( ncnew, varlist$gcm, "add_offset", offset, prec="float" )
  ncatt_put( ncnew, varlist$gcm, "scale_factor", scale, prec="float" ) 
  ncatt_put( ncnew, varlist$gcm, "_FillValue", missval, prec="float" ) 
  ncatt_put( ncnew, varlist$gcm, "missing_value", missval, prec="float" )
  ## GCM names
  if (verbose) print('write GCM names')  
  ncatt_put( ncnew, varlist$gcm, "GCM runs", paste(names.x,collapse=','),prec="char" )
  ## PCA/EOF variables:
  if (verbose) print('EOF/PCA variables')
  ncvar_put( ncnew, varlist$pca, round((pca - offset)/scale) )
  ncatt_put( ncnew, varlist$pca, "add_offset", offset, prec="float" )
  ncatt_put( ncnew, varlist$pca, "scale_factor", scale, prec="float" ) 
  ncatt_put( ncnew, varlist$pca, "_FillValue", missval, prec="float" ) 
  ncatt_put( ncnew, varlist$pca, "missing_value", missval, prec="float" ) 
  ncvar_put( ncnew, varlist$pat, round((pattern - offset)/scale) )
  ncatt_put( ncnew, varlist$pat, "add_offset", offset, prec="float" )
  ncatt_put( ncnew, varlist$pat, "scale_factor", scale, prec="float" ) 
  ncatt_put( ncnew, varlist$pat, "_FillValue", missval, prec="float" ) 
  ncatt_put( ncnew, varlist$pat, "missing_value", missval, prec="float" ) 
  ncatt_put( ncnew, varlist$pat, "dimensions_pca", paste(attr(pattern,'dimensions'),collapse=', '), prec="char" )
  ## If the object contains PCAs
  if (length(dpat)==2) {
    if (verbose) print('PCA only variables')
    ncvar_put( ncnew, varlist$lon, lon(x) )
    ncvar_put( ncnew, varlist$lat, lat(x) )
    ncvar_put( ncnew, varlist$alt, alt(x) )
    ncvar_put( ncnew, varlist$loc, loc(x) )
    ncvar_put( ncnew, varlist$stid, as.character(stid(x)) )
    ncvar_put( ncnew, varlist$cntr, cntr(x) )
    ncvar_put( ncnew, varlist$src, src(x) )
    ncvar_put( ncnew, varlist$lambda, attr(pca,'eigenvalues') )
  } else if (length(dpat)==3)
    ncvar_put( ncnew, varlist$lambda, attr(pca,'eigenvalues') )
  if (verbose) print('history')
  ncatt_put( ncnew, varlist$pca, "history", paste(attr(x,'history'),collapse=';'), prec="char" )  
  ## Global attributes:
  ncatt_put( ncnew, 0, 'class', class(x))
  ncatt_put( ncnew, 0, "description", "Saved from esd using write2ncdf4.dsensemble")
  #ncatt_put( ncnew, 0, "class", paste(class.x,collapse='-'))
  ncatt_put( ncnew, 0, "esd-version", attr(x,'history')$session$esd.version)
  nc_close(ncnew)
  if (verbose) print(paste('Finished sucessfully - file', fname))  
}
