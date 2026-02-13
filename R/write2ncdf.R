#' Saves climate data as netCDF.
#' 
#' Method to save data as netCDF, making sure to include the data
#' structure and meta-data (attributes). The code tries to follow the netCDF
#' 'CF' convention (\url{https://cfconventions.org/}). The method is built on the \code{ncdf4} package.
#' 
#' @seealso write2ncdf4.station nc4combine.stations 
#' write2ncdf4.field write2ncdf4.list
#' write2ncdf4.station write2ncdf4.dsensemble
#' write2ncdf4.eof write2ncdf4.pca
#' 
#' @param x data object
#' @param \dots additional arguments
#' 
#' @return None
#' 
#' @keywords netcdf ncdf4 save
#' 
#' @examples
#' 
#' \dontrun{
#' nacd <- station(src='nacd')
#' X <- annual(nacd)
#' write2ncdf4(X,file='test.nc')
#' }
#' 
#' @export write2ncdf4
write2ncdf4 <- function(x,...) UseMethod("write2ncdf4")

#' @exportS3Method
#' @export
write2ncdf4.default <- function(x,...) {
  args <- list(...)
  if (!is.null(args$verbose)) verbose <- args$verbose else verbose <- FALSE
  if (!is.null(args$missval)) missval <- args$missval else missval <- -32767
  if (is.null(args$file)) stop("Need the argument 'file'") else file <- args$file 
  if (verbose) cat('write2ncdf4.default ',names(args),' \n')
  if (verbose) cat(dim(x),'\n')
  
  dimlon <- ncdim_def( "longitude", "degree_east", lon(x) )
  dimlat <- ncdim_def( "latitude", "degree_north", lat(x) )
  varid <- varid(x)
  unit <- unit(x)
  if (verbose)  cat(varid,unit,'\n')
  x4nc <- ncvar_def(varid, unit, list(dimlon,dimlat), -1, 
                    longname=attr(x,'longname'), prec='float')
  
  # Create a netCDF file with this variable
  ncnew <- nc_create( file, x4nc )
  x[!is.finite(x)] <- missval
  ncvar_put( ncnew, varid, x )
  ncatt_put( ncnew, varid, "_FillValue", missval, prec="float" ) 
  ncatt_put( ncnew, varid, "missing_value", missval, prec="float" ) 
  history <- toString(attr(x,'history')$call)
  ncatt_put( ncnew, varid, "history", history, prec="text" ) 
  ncatt_put( ncnew, 0, 'class', class(x))
  ncatt_put( ncnew, 0, "description", 
             paste("Saved from esd using write2ncdf4",date()))
  ncatt_put( ncnew, 0, "esd-version", attr(x,'history')$session$esd.version)
  add_ADC_meta(ncnew,x,conventions=NA,title=NA,summary=NA,project=NA,license=NA,verbose=verbose)
  nc_close(ncnew)
  if (verbose) print('netCDF file saved')
}

#' Saves climate data as netCDF.
#' 
#' Method to save data as netCDF, making sure to include the data
#' structure and meta-data (attributes). The code tries to follow the netCDf
#' 'CF' convention. The method is built on the \code{ncdf4} package.
#' 
#' @aliases write2ncdf4.field
#' @seealso write2ncdf4
#' 
#' @param x data object
#' @param file file name
#' @param prec Precision: see \code{\link[ncdf4]{ncvar_def}}
#' @param scale Sets the atttribute 'scale_factor' which is used to scale
#' (multiply) the values stored (to save space may be represented as 'short').
#' @param offset Sets the attribute 'add_offset' which is added to the values
#' stored (to save space may be represented as 'short').
#' @param torg Time origin
#' @param missval Missing value: see \code{\link[ncdf4]{ncvar_def}}
#' @param verbose if TRUE print progress
#' @param start defines the start of the year (default is 'January')
#' @param \dots additional arguments
#' 
#' @return None
#' 
#' @keywords netcdf ncdf4 save
#' 
#' @exportS3Method
#' @export write2ncdf4.list
write2ncdf4.list <- function(x,...,file='field.nc',prec='short',scale=0.1,offset=NULL,
                             torg="1970-01-01",missval=-32767,verbose=FALSE) {
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
  units <- unlist(lapply(x,function(x) attr(x,"unit")[1]))
  if (verbose) {print(varids); print(units); print(n)}
  x4nc <- list()
  for (i in 1:n) {
    if (verbose) print(paste('write2ncdf4.list: ',i,'ncvar_def',varids[i]))
    if (verbose) print(paste('write2ncdf4.list: longname=',attr(x[[i]],'longname')))
    x4nc[[varids[i]]] <- ncvar_def(varids[i], units[i], list(dimlon,dimlat,dimtim), -1, 
                                   longname=attr(x[[i]],'longname'), prec=prec[i])
  }
  
  # Create a netCDF file with this variable
  ncnew <- nc_create( file, x4nc )
  
  # Write some values to this variable on disk.
  for (i in 1:n) {
    if (verbose) print(names(x)[i])
    y <- coredata(x[[i]])
    if (is.null(offset[i])) offset[i] <- mean(y,na.rm=TRUE)
    if (is.null(scale[i])) scale[i] <- 1
    y <- t(y)
    y <- round((y-offset[i])/scale[i])
    y[!is.finite(y)] <- missval
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
  add_ADC_meta(ncnew,x,conventions=NA,title=NA,summary=NA,project=NA,license=NA,verbose=verbose)
  nc_close(ncnew)
  if (verbose) print('netCDF file saved')
}

#' @exportS3Method
#' @export write2ncdf4.field
write2ncdf4.field <- function(x,...,file='field.nc',prec='short',scale=NULL,offset=NULL,
                              torg="1970-01-01",missval=-32767,verbose=FALSE) {
  if (verbose) {print('write2ncdf4.field'); print(names(attributes(x)))}
  
  y <- coredata(x)
  if (is.null(offset)) offset <- mean(y,na.rm=TRUE)
  if (is.null(scale)) scale <- (max(abs(c(y)),na.rm=TRUE) - offset)/10000
  y <- t(y)
  y <- round((y-offset)/scale)
  y[!is.finite(y)] <- missval
  if (verbose) {
    print(attr(y,'dimensions')); print(c(scale,offset))
    #print(range(c(y))); print(range(c(x),na.rm=TRUE))
  }
  dim(y) <- attr(x,'dimensions')
  
  dimlon <- ncdim_def( "longitude", "degree_east", lon(x) )
  dimlat <- ncdim_def( "latitude", "degree_north", lat(x) )
  if (inherits(index(x),c('numeric','integer')))
    index(x) <- as.Date(paste(index(x),'-01-01',sep=''))
  
  dimtim <- ncdim_def( "time", paste("days since",torg),
                       as.numeric(as.Date(index(x),origin=torg)) )
  if (verbose) print(paste('write2ncdf4.field: longname=',attr(x,'longname')[1]))
  x4nc <- ncvar_def(varid(x)[1], attr(x,"unit")[1], list(dimlon,dimlat,dimtim), -1, 
                    longname=attr(x,'longname')[1], prec=prec)
  
  # Create a netCDF file with this variable
  ncnew <- nc_create( file, x4nc )
  
  # Write some values to this variable on disk.
  ncvar_put( ncnew, x4nc, round(y) )
  ncatt_put( ncnew, x4nc, "add_offset", offset, prec="float" )
  ncatt_put( ncnew, x4nc, "scale_factor", scale, prec="float" ) 
  ncatt_put( ncnew, x4nc, "_FillValue", missval, prec="float" ) 
  ncatt_put( ncnew, x4nc, "missing_value", missval, prec="float" ) 
  history <- toString(attr(x,'history')$call)
  ncatt_put( ncnew, x4nc, "history", history, prec="text" ) 
  #ncatt_put( ncnew, 0, 'class', class(x))
  ncatt_put( ncnew, 0, "description", 
             paste("Saved from esd using write2ncdf4",date()))
  if (verbose) print(attr(x,'history'))
  ncatt_put( ncnew, 0, "esd-version", attr(x,'history')$session$esd.version)
  
  # Add the attributes of object x
  attnames <- names(attributes(x))
  if (verbose) print(attnames)
  attnames <- attnames[-grep('history',attnames)]
  attnames <- attnames[-grep('unit',attnames)]
  attnames <- attnames[-grep('variable',attnames)]
  attnames <- attnames[-grep('dim',attnames)]
  attnames <- attnames[-grep('index',attnames)]
  attnames <- attnames[-grep('longitude',attnames)]
  attnames <- attnames[-grep('latitude',attnames)]
  if (length(grep('greenwich',attnames))>0) 
    attnames <- attnames[-grep('greenwich',attnames)]
  if (length(grep('call',attnames))>0) 
    attnames <- attnames[-grep('call',attnames)]
  for (ia in 1:length(attnames)) {
    if (verbose) print(paste0("write2ncdf4.field: ",attnames[ia],": ", attr(x,attnames[ia])))
    ## REB 2025-06-16
    if ( (!is.null(attr(x,attnames[ia]))) & (length(attr(x,attnames[ia])>0)) ) 
      ncatt_put( ncnew, 0, attnames[ia], as.character(attr(x,attnames[ia])), 
                 prec="text")
  }
  add_ADC_meta(ncnew,x,conventions=NA,title=NA,summary=NA,project=NA,license=NA,verbose=verbose)
  nc_close(ncnew)
}

#' Saves climate data as netCDF.
#' 
#' Method to save station data as netCDF, making sure to include the data
#' structure and meta-data (attributes). The code tries to follow the netCDf
#' 'CF' convention. The method is built on the \code{ncdf4} package.
#'
#' To save space, the values are saved as short (16-bit signed integer that
#' can hold values between -32768 and 32767).
#' (see NC_SHORT in \url{https://www.unidata.ucar.edu/software/netcdf/docs/data_type.html}).
#'
#' @seealso write2ncdf4
#' 
#' @param x data object
#' @param file file name
#' @param prec Precision: see \code{\link[ncdf4]{ncvar_def}}
#' @param scale Sets the attribute 'scale_factor' which is used to scale
#' (multiply) the values stored (to save space may be represented as 'short').
#' @param it a time index, see \code{\link{subset}}
#' @param stid station id
#' @param append a Boolean; if TRUE append output to existing file - not used currently
#' @param stid_unlim a Boolean; if TRUE the stid dimension is unlimited
#' @param namelength a numeric specifying the number of characters in dimension and variable names
#' @param nmin Only calculate summary statistics for stations with nmin years of data (e.g. 30 years).
#' @param offset Sets the attribute 'add_offset' which is added to the values
#' stored (to save space may be represented as 'short').
#' @param torg Time origin
#' @param missval Missing value: see \code{\link[ncdf4]{ncvar_def}}
#' @param verbose TRUE - clutter the screen.
#' @param doi - Data ID. All the following arguments are meant to accommodate for the convention described 
#' at \url{https://adc.met.no/node/4}
#' @param namingauthority The organization that provides the initial id (see above) for the dataset (see \url{https://adc.met.no/node/4}).
#' @param processinglevel A textual description of the processing (or quality control) level of the data. (see \url{https://adc.met.no/node/4}).
#' @param creatortype Specifies type of creator with one of the following: 'person', 'group', 'institution', or 'position'. If this attribute is not specified, the creator is assumed to be a person. (see \url{https://adc.met.no/node/4}).
#' @param creatoremail The email address of the person (or other creator type specified by the creator_type attribute) principally responsible for creating this data. (see \url{https://adc.met.no/node/4}).
#' @param institution The name of the institution principally responsible for originating this data. (see \url{https://adc.met.no/node/4}).
#' @param publishername The name of the person (or other entity specified by the publisher_type attribute) responsible for publishing the data file or product to users, with its current metadata and format. (see \url{https://adc.met.no/node/4}).
#' @param publisheremail The email address of the person (or other entity specified by the publisher_type attribute) responsible for publishing the data file or product to users, with its current metadata and format. (see \url{https://adc.met.no/node/4}).
#' @param publisherurl The URL of the person (or other entity specified by the publisher_type attribute) responsible for publishing the data file or product to users, with its current metadata and format. (see \url{https://adc.met.no/node/4}).
#' @param project The name of the project(s) principally responsible for originating this data. (see \url{https://adc.met.no/node/4}).
#' @param nspc The size of chunks of data processes sequentially in order to limit the need of computer memory. 
#' Smaller number requires less memory.
#' @param start specifies the start month of the year (default is January)
#' @param \dots additional arguments
#' 
#' @return None
#' 
#' @keywords netcdf ncdf4 save
#' 
#' @exportS3Method
#' @export write2ncdf4.station
write2ncdf4.station <- function(x,...,file='station.nc',prec='short',offset=0, 
                                missval=-32767,it=NULL,stid=NULL,append=FALSE,
                                scale=0.1,torg='1899-12-31',stid_unlim=FALSE,namelength=24,nmin=30,verbose=FALSE,
                                doi='NA',namingauthority='NA',processinglevel='NA',creatortype='NA',
                                creatoremail='NA',institution='NA',publishername='NA',publisheremail='NA',
                                publisherurl='NA',project='NA',nspc=100,start='Jan') {
  
  if (!inherits(x,"station")) stop('x argument must be a station object') 
  
  if (verbose) {
    print('write2ncdf4.station'); print(range(index(x)))
    print(class(x)); print(dim(x)); print(range(c(coredata(x))))
    #print(range(c(coredata(x)),na.rm=TRUE)); print('---')
  }
  ## Quality check - remove obviously unrealistic values
  if (is.precip(x)) {
    #  cx <- coredata(x); cx[cx < 0] <- NA; cx[cx > 1500] <- NA; cx -> coredata(x); rm('cx') 
    x[coredata(x) <0] <- NA
  }
  if (is.T(x)) {
    #  cx <- coredata(x); cx[cx < -100] <- NA; cx[cx > 100] <- NA; cx -> coredata(x); rm('cx')
    x[abs(coredata(x)) > 100] <- NA
  }
  
  #cx <- coredata(x); cx[cx < -100] <- NA;
  x[coredata(x) < -180] <- NA
  if (prec=='short') x[coredata(x) > 3200] <- NA;
  #cx -> coredata(x); rm('cx') 
  
  ## Don't save empty space:
  x0 <- x
  if (length(dim(x))==2) {
    x[coredata(x) <= missval] <- NA
  }
  if (sum(is.finite(x))==0) {
    print('write2ncdf4.station: No valid data after weeding missing values')
    print(summary(x0)); print(summary(x))
    return()
  }
  rm('x0'); gc(reset=TRUE)
  ## Weed out stations with short time series:
  if (length(dim(x))==2) {
    if (verbose) print('Excude short series')
    good <- apply(coredata(x),2,FUN='nv')
    #x <- subset(x,is=good > 365) # doesn't work for annual/seasonal data
    if(inherits(x,c("annual","season"))) {
      x <- subset(x,is=good > 5)
    } else if (inherits(x,"month")) {
      x <- subset(x,is=good > 24)
    } else {
      x <- subset(x,is=good > 365)
    }
  }
  
  if (length(dim(x))==2) {
    good <- apply(coredata(x),1,FUN='nv') 
  } else {
    good <- nv(x)
  }
  
  if (is.null(it)) x <- subset(x,it=good>0)
  if (verbose) {print('time period after missing data have been removed'); print(range(index(x)))}
  
  ## Write a station object as a netCDF file using the short-type combined with add_offsetet and scale_factor
  ## to reduce the size.   
  ## Examine the station object: dimensions and attributes  
  
  ## Get time interval: 
  if (is.null(it)) nt <- dim(x)[1] else nt <- max(it) - min(it) + 1
  
  ## Determine selection of stations and umber of stations ns 
  if (!is.null(stid))  {
    if (!is.null(dim(x))) nstations <- dim(x)[2] else if (!is.null(x)) nstations <- 1 else nstations <- 0
    # if (append & (length(stid) != nstations)) {
    #   print(dim(x)); print(length(stid))
    #   stop('write2ncdf4.station: stid argument does not match x')
    # }
    ns <- length(stid) 
  } else if (length(dim(x))==2) {
    ns <- dim(x)[2] 
    stid <- 1:ns
  } else {
    ns <- 1
    dim(x) <- c(length(x),ns)
  }
  if (verbose) print(c(nt,ns))
  
  ## if (is.null(d)) d <- c(length(x),1)
  if (verbose) print(paste('Number of stations to save: ',paste(ns)))
  
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
  
  ## Compute summary statistics for the stations, e.g. mean, max, trend, etc.
  x0 <- x; missval0 <- missval; verbose0 <- verbose
  stats <- StationSumStats(x=x,missval=missval,ns=nspc,verbose=verbose,start=start)
  #list2env(StationSumStats(x=x,missval=missval,ns=nspc,verbose=verbose,start=start),envir=environment())
  verbose <- verbose[1]  ## There seemed to be several versions of 'verbose'
  x <- x0; missval <- missval0; rm('x0','missval0'); gc(reset=TRUE) ## REB in case something happened to x in the function call above...
  #>>>>>>> master
  if (verbose) print('Summary statistics computed')
  ## Only do summary statistics for stations with more than 30 years
  insufficient <- apply(coredata(x),2,nv) < nmin*365
  if (verbose) print(summary(insufficient))
  
  if (is.null(attr(x,'calendar'))) {
    attr(x,'calendar') <- 'standard'
  }
  
  if (!is.na(attr(x,'calendar'))) {
    calendar <- attr(x,'calendar')
  } else {
    calendar <- 'standard'
  }
  ## REB 2023-08-09: The lines for generating the netCDF file have been moved into a separate function for reuse...
  generate.station.ncfile(x,file,stats,missval,offset,scale,torg,prec,it,verbose,stid_unlim,
                          namelength,calendar,doi,namingauthority,processinglevel,creatortype,
                          creatoremail,institution,publishername,publisheremail,
                          publisherurl,project,insufficient)
}


#' Saves climate data as netCDF.
#' 
#' Method to save 'pca' data as netCDF, making sure to include the data
#' structure and meta-data (attributes). The code tries to follow the netCDf
#' 'CF' convention. The method is built on the \code{ncdf4} package.
#'
#' To save space, the values are saved as short (16-bit signed integer that
#' can hold values between -32768 and 32767).
#' (see NC_SHORT in \url{https://www.unidata.ucar.edu/software/netcdf/docs/data_type.html}).
#'
#' @seealso write2ncdf4
#' 
#' @param x data object
#' @param file file name
#' @param prec Precision: see \code{\link[ncdf4]{ncvar_def}}
#' @param scale Sets the atttribute 'scale_factor' which is used to scale
#' (multiply) the values stored (to save space may be represented as 'short').
#' @param offset Sets the attribute 'add_offset' which is added to the values
#' stored (to save space may be represented as 'short').
#' @param missval Missing value: see \code{\link[ncdf4]{ncvar_def}}
#' @param verbose TRUE - clutter the screen.
#' @param \dots additional arguments
#' 
#' @return None
#' 
#' @keywords netcdf ncdf4 save
#'
#' @exportS3Method
#' @export write2ncdf4.pca
write2ncdf4.pca <- function(x,...,file='esd.pca.nc',prec='short',verbose=FALSE,
                            scale=0.01,offset=0,missval=-32767) {
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
  ## KMP 2018-11-08: tim not defined
  tim <- ncvar_def("time", "year", dimtim, missval, 
                   longname='time', prec='float')
  lon <- ncvar_def("longitude", "degree_east", dimxy, missval, 
                   longname='longitude', prec='float')
  lat <- ncvar_def("latitude", "degree_north", dimxy, missval, 
                   longname='latitude', prec='float')
  alt <- ncvar_def("altitude", "m", dimxy, missval, 
                   longname='altitude', prec='float')
  stid <- ncvar_def("station_id", "number", dimxy, missval, 
                    longname='station ID', prec="integer")
  lambda <- ncvar_def("lambda", "number", dimpca, missval, 
                      longname='eigenvalues', prec="float")
  if (is.numeric(index(x))) index(x) <- year(x)
  dpca <- dim(pca); attributes(pca) <- NULL; dpca -> dim(pca)
  ## KMP 2018-11-08: nc has not been defined! should it be created or opened from a file?
  nc <- nc_create(file,vars=list(lon,lat,alt,stid,pca,pat,lambda))
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
  ncvar_put( nc, stid, stid(x) )
  ncvar_put( nc, lambda, attr(x,'eigenvalues') )
  ncatt_put( nc, pca, "history", paste(attr(x,'history'),collapse=';'), prec="char" )
  ncatt_put( nc, 0, 'class', class(x))
  ncatt_put( nc, 0, "esd-version", attr(x,'history')$session$esd.version)
  add_ADC_meta(nc,x,conventions=NA,title=NA,summary=NA,project=NA,license=NA,verbose=verbose)
}

#' Unfinished function that doesn't do anything.
#'
#' @param x input object of class 'dsensemble'
#' @param verbose if TRUE print progress
#' @param \dots additional arguments
#'
#' @seealso write2ncdf4
#' @exportS3Method
#' @export write2ncdf4.eof
write2ncdf4.eof <- function(x,...,verbose=FALSE){
  if(verbose) print("write2ncdf.eof")
  if(verbose) print("unfinished function that doesn't do anything")
}

## Recursive helping function to simplify the coding of write2ncdf4.station and alleviate problems with memory usage
## REB 2022-03-15
combinelist <- function(x,y,verbose=FALSE) {
  if (verbose) print('combinelist')
  ln <- names(x); ln2 <- length(names(y))
  if (verbose) {print(ln); print(ln2); print(length(x[[1]]))}
  for (i in 1:length(ln)) x[[ln[i]]] <- c(x[[ln[i]]],y[[ln[i]]])
  if (verbose) print(str(x))
  return(x)
}

daysold <- function(x,t) {
  ## This function returns the number of days between the last observation and the time now
  t <- t[is.finite(x)]
  lastdate.x <- max(t)
  y <- Sys.Date() - lastdate.x
  if (y < 0) print(lastdate.x)
  return(y)
}

StationSumStats <- function(x,missval=-32767,ns=300,verbose=FALSE,start='Jan') {
  if (verbose) print(paste('StationSumStats - precip?',is.precip(x)))
  
  if (is.null(dim(x))) {
    if (verbose) print('Single station - change dim(x)')
    dim(x) <- c(length(x),1)
  }
  
  d <- dim(x)
  ## REB 2022-03-31
  ## Make sure that this algorithm keeps track of the type of data (precip, temp, etc)
  if ( (length(esd::unit(x))!=d[2]) | (length(esd::varid(x))!=d[2]) ) {  
    attr(x,'unit') <- rep(esd::unit(x)[1],d[2])
    attr(x,'variable') <- rep(esd::varid(x)[1],d[2])
  }
  Y <- NULL
  if (d[2] > ns) {
    if (verbose) print(d)
    ## Split the station data into smaller chunks of ns stations to cope with large data object memory-wise
    for (is in seq(1,d[2],by=ns)) {
      if (verbose) print(paste('is =',is))
      js <- seq(is,min(c(is+ns-1,d[2])),by=1)
      if (verbose) {
        print(paste('is.precip(x) =',is.precip(x), paste(dim(x),collapse=' x '),
                    paste(range(js),collapse='-')))
        print(is.precip(subset(x,is=js))); str(subset(x,is=js,verbose=TRUE))
      }
      y <- StationSumStats(subset(x,is=js),missval=missval,verbose=verbose)
      if (verbose) print(names(y))
      if (is.null(Y)) Y <- y else Y <- combinelist(Y,y,verbose=verbose)
    }
    #if (verbose) print(str(Y))
    return(Y)
  } else { 
    
    ## Different summary statistics for precipitation, temperature and other variables.
    
    ## If there are more than one stations, use matrix oprations
    ## Estimate summary statistics for the station data
    if (verbose) print(paste('Estimate summary statistics - data dimension:',
                             paste(d,collapse='x'),' - precipitation?',is.precip(x)))
    #browser()
    mx <- apply(x,2,'max',na.rm=TRUE)  ## Maximum
    mn <- apply(x,2,'min',na.rm=TRUE)  ## Minimum
    # nhr <- apply(anomaly(x),2,'arec')  ## Record-high statistics: on anomalies
    # nlr <- apply(-anomaly(x),2,'arec') ## Record-low statistics
    ## REB 2023-06-20: the previous lines of codes crashed
    nhr <- apply(x,2,'arec')  ## Record-high statistics: on the absolute data
    daysold <- apply(x,2,'daysold',index(x))
    if (!is.precip(x)) nlr <- rep(NA,dim(x)[2]) else nlr <- apply(-x,2,'arec') ## Record-low statistics
    
    if (!is.precip(x)) {
      if (verbose) print('*** Not precipitation ***')
      if (verbose) print('Records')
      #nlr <- apply(-x,2,'arec') ## Record-low statistics
      ave <- apply(x,2,'mean',na.rm=TRUE)
      if (verbose) print('averages')
      ave.djf <- apply(subset(x,it='djf'),2,'mean',na.rm=TRUE)
      ave.mam <- apply(subset(x,it='mam'),2,'mean',na.rm=TRUE)
      ave.jja <- apply(subset(x,it='jja'),2,'mean',na.rm=TRUE)
      ave.son <- apply(subset(x,it='son'),2,'mean',na.rm=TRUE)
      if (verbose) print('standard deviations of anomalies')
      std <- apply(anomaly(x),2,'sd',na.rm=TRUE)
      std.djf <- apply(subset(anomaly(x),it='djf'),2,'sd',na.rm=TRUE)
      std.mam <- apply(subset(anomaly(x),it='mam'),2,'sd',na.rm=TRUE)
      std.jja <- apply(subset(anomaly(x),it='jja'),2,'sd',na.rm=TRUE)
      std.son <- apply(subset(anomaly(x),it='son'),2,'sd',na.rm=TRUE)
      if (verbose) print('trends')
      td <- apply(annual(x,start=start),2,'trend.coef',start=start)
      td.djf <- apply(annual(subset(x,it='djf'),'mean',nmin=75),2,'trend.coef')
      td.mam <- apply(annual(subset(x,it='mam'),'mean',nmin=75),2,'trend.coef')
      td.jja <- apply(annual(subset(x,it='jja'),'mean',nmin=75),2,'trend.coef')
      td.son <- apply(annual(subset(x,it='son'),'mean',nmin=75),2,'trend.coef')
    } else if (is.precip(x)) {
      if (verbose) print('ooo Precipitation ooo')
      if (verbose) print('sums')
      ave <- apply(annual(x,'sum',start=start),2,'mean',na.rm=TRUE)
      ave.djf <- apply(annual(subset(x,it='djf'),'sum',nmin=90),2,'mean',na.rm=TRUE)
      ave.mam <- apply(annual(subset(x,it='mam'),'sum',nmin=90),2,'mean',na.rm=TRUE)
      ave.jja <- apply(annual(subset(x,it='jja'),'sum',nmin=90),2,'mean',na.rm=TRUE)
      ave.son <- apply(annual(subset(x,it='son'),'sum',nmin=90),2,'mean',na.rm=TRUE)
      if (verbose) print('wet-day mean')
      mu <- apply(annual(x,'wetmean',start=start),2,'mean',na.rm=TRUE)
      mu.djf <- apply(annual(subset(x,it='djf'),'wetmean',nmin=75),2,'mean',na.rm=TRUE)
      mu.mam <- apply(annual(subset(x,it='mam'),'wetmean',nmin=75),2,'mean',na.rm=TRUE)
      mu.jja <- apply(annual(subset(x,it='jja'),'wetmean',nmin=75),2,'mean',na.rm=TRUE)
      mu.son <- apply(annual(subset(x,it='son'),'wetmean',nmin=75),2,'mean',na.rm=TRUE)
      if (verbose) print('wet-day frequency')
      fw <- apply(100*annual(x,'wetfreq',start=start),2,'mean',na.rm=TRUE)
      fw.djf <- apply(100*annual(subset(x,it='djf'),'wetfreq',nmin=75),2,'mean',na.rm=TRUE)
      fw.mam <- apply(100*annual(subset(x,it='mam'),'wetfreq',nmin=75),2,'mean',na.rm=TRUE)
      fw.jja <- apply(100*annual(subset(x,it='jja'),'wetfreq',nmin=75),2,'mean',na.rm=TRUE)
      fw.son <- apply(100*annual(subset(x,it='son'),'wetfreq',nmin=75),2,'mean',na.rm=TRUE)
      td <- apply(annual(x,FUN='sum',start=start),2,'trend.coef')
      td.djf <- apply(annual(subset(x,it='djf'),'sum',nmin=90),2,'trend.coef')
      td.mam <- apply(annual(subset(x,it='mam'),'sum',nmin=90),2,'trend.coef')
      td.jja <- apply(annual(subset(x,it='jja'),'sum',nmin=90),2,'trend.coef')
      td.son <- apply(annual(subset(x,it='son'),'sum',nmin=90),2,'trend.coef')
      tdfw <- apply(100*annual(x,FUN='wetfreq',start=start),2,'trend.coef')
      tdfw.djf <- apply(100*annual(subset(x,it='djf'),'wetfreq',nmin=75),2,'trend.coef')
      tdfw.mam <- apply(100*annual(subset(x,it='mam'),'wetfreq',nmin=75),2,'trend.coef')
      tdfw.jja <- apply(100*annual(subset(x,it='jja'),'wetfreq',nmin=75),2,'trend.coef')
      tdfw.son <- apply(100*annual(subset(x,it='son'),'wetfreq',nmin=75),2,'trend.coef')
      tdmu <- apply(annual(x,FUN='wetmean',start=start),2,'trend.coef')
      tdmu.djf <- apply(annual(subset(x,it='djf'),'wetmean',nmin=75),2,'trend.coef')
      tdmu.mam <- apply(annual(subset(x,it='mam'),'wetmean',nmin=75),2,'trend.coef')
      tdmu.jja <- apply(annual(subset(x,it='jja'),'wetmean',nmin=75),2,'trend.coef')
      tdmu.son <- apply(annual(subset(x,it='son'),'wetmean',nmin=75),2,'trend.coef')
      lr <- sapply(x,'lastrains')
      ld <- sapply(x,'lastdry')
      if (verbose) print('Sigma2')
      sigma2 <- apply(x,2,'rainvar')
      sigma2.djf <- apply(subset(x,it='djf'),2,'rainvar')
      sigma2.mam <- apply(subset(x,it='mam'),2,'rainvar')
      sigma2.jja <- apply(subset(x,it='jja'),2,'rainvar')
      sigma2.son <- apply(subset(x,it='son'),2,'rainvar')
      tsigma2 <- rainvartrend(x)
      tsigma2.djf <- rainvartrend(subset(x,it='djf'),nmin=90)
      tsigma2.mam <- rainvartrend(subset(x,it='mam'),nmin=90)
      tsigma2.jja <- rainvartrend(subset(x,it='jja'),nmin=90)
      tsigma2.son <- rainvartrend(subset(x,it='son'),nmin=90)
      ## Mean wet/dry-spell length
      if (verbose) print('Spell')
      t <- index(x)
      ss <- try(spell(x,threshold=1))
      ## If spell returns NULL, then return NAs...
      if ( (inherits(ss,'spell')) & (!is.null(ss)) ) { 
        mwsl <- colMeans(subset.station(ss,is=list(param='wet')),na.rm=TRUE)
        mdsl <- colMeans(subset.station(ss,is=list(param='dry')),na.rm=TRUE)
      } else {mwsl <- rep(missval,length(ave)); mdsl <- rep(missval,length(ave))}
    }
    
    if (verbose) print('lastelementrecord')
    ## Is the last element a high or low record?
    lehr <- lastelementrecord(x,verbose=verbose)
    lelr <- lastelementrecord(-x)
    
    rm('x','t','ss','d','is','ns','Y', 'missval')
    
    if (verbose) {ls(); print('StationSumStats returns all results in a list object...')}
    if (verbose) str(mget(ls()))
    ## Return all the objects defined in this environment:
    return(mget(ls()))
  }
}


generate.station.ncfile <- function(x,file,stats,missval,offset,scale,torg,prec='short',it=NULL,verbose=FALSE,
                                    stid_unlim=FALSE,namelength=24,calendar='standard',
                                    doi='NA',namingauthority='NA',processinglevel='NA',creatortype='NA',
                                    creatoremail='NA',institution='NA',publishername='NA',publisheremail='NA',
                                    publisherurl='NA',project='NA',insufficient) {
  ## Put the data into y in the form of a matrix with scaled values around zero and an offset: 
  if (verbose) print('<generate.station.ncfile>')
  src <- paste(rownames(table(attr(x,"source"))),collapse="/")
  if (verbose) print(src)
  y <- coredata(x)
  nt <- length(index(x))
  list2env(stats,envir=environment())
  verbose <- verbose[1]
  unitx <- attr(x,'unit')
  if (verbose) print(ls())
  if (!is.null(dim(y))) ns <- dim(y)[2] else ns <- 1
  if (verbose) {print('Scale the data after removing offset'); print(str(x)); print(summary(c(y))); 
    print(paste(sum(is.finite(y)),'good data &',sum(!is.finite(y)),'bad data')); print(c(offset,scale,missval))}
  
  y <- round((y - offset)/scale)
  if (verbose) {print('After scaling'); print(dim(y)); print(summary(c(y)))}
  y[!is.finite(y)] <- missval
  
  if (class(index(x))=='Date') {
    if (is.null(it)) {
      time <- julian(index(x)) - julian(as.Date(torg)) 
    } else {
      ## If it is provided, then ensure that the stations are saved for this time interval
      ## pad with NAs if necessary...
      if (verbose) print(paste('Use prescribed time coordinates',it[1],'-',it[2]))
      if (is.numeric(it)) {
        ## Assume its given as year - change to date
        it <- as.Date(c(paste0(it[1],'-01-01'),paste0(it[2],'-12-31')))
      }
      if (length(it)==2) it <- seq(min(it),max(it),by=1)
      y <- zoo(y,order.by=index(x))
      #print(dim(y))
      x2 <- merge(zoo(rep(0,nt),order.by=seq(min(it),max(it),length=nt)),zoo(y),all=TRUE)
      x2 <- window(x2[,-1],start=it[1],end=it[length(it)])
      x2 <- attrcp(x,x2); class(x2) <- class(x); y <- x2; #rm('x2')
      #print(dim(y))
      time <- julian(it) - julian(as.Date(torg))
      it <- index(y)
    }
    if (verbose) cat('index(x):',range(index(x)),' it:',range(it),' dim(x)',dim(x),'\n')
  } else if (inherits(x,'annual')) {
    time <- julian(as.Date(paste(year(x),'01-01',sep='-')))-julian(as.Date(torg))
  }
  #plot(y,new=FALSE); browser()
  
  if(verbose) print(paste('Period in data: ',min(firstyear(x)),' - ', max(lastyear(x)),' and time dimension: ',
                          paste(range(as.Date(time,origin=torg)),collapse=' - ')))
  
  # Attributes with same number of elements as stations are saved as variables
  if (is.null(it)) it <- index(y)
  
  # Define the dimensions: create a new file - do not append existing one
  #if (!append) {
  if (verbose) print('Define dimensions')
  if (verbose) print(stid(x))
  dimS <- ncdim_def( name="stid", units="number",vals=1:ns,unlim=stid_unlim)
  dimT <- ncdim_def( name="time", units=paste("days since",torg), vals=time, calendar=calendar,unlim=TRUE)
  dimnchar   <- ncdim_def("nchar",   "", 1:namelength, create_dimvar=FALSE )
  
  if (verbose) {
    print('Define variable')
    print(paste('create netCDF-file',file))
    print(summary(c(y)))
  }
  
  if (verbose) print('Define the netCDF structure')
  latid <- ncvar_def(name="lat",dim=list(dimS), units="degrees_north", missval=missval,longname="latitude", 
                     prec="float",verbose=verbose)
  
  lonid <- ncvar_def(name="lon",dim=list(dimS), units="degrees_east", missval=missval,longname="longitude", 
                     prec="float",verbose=verbose)
  altid <- ncvar_def(name="alt",dim=list(dimS), units="meters", missval=missval,longname="altitude", 
                     prec=prec,verbose=verbose)
  
  locid <- ncvar_def(name="loc",dim=list(dimnchar,dimS),units="name",prec="char",longname="location",
                     verbose=verbose)
  stid <- ncvar_def(name="stationID",dim=list(dimnchar,dimS),units="number",prec="char",longname="station_id",
                    verbose=verbose)
  cntrid <- ncvar_def(name="cntr",dim=list(dimnchar,dimS),units="name",prec="char",longname="country",
                      verbose=verbose)
  
  if (verbose) {
    print(paste('ncvar:',varid(x)[1]))
    print(unitx[1])
  }
  ## KMP 2018-11-02: devtools (run_examples) can only handle ASCII characters so I had to replace the 
  ## degree symbol with "\u00B0", but I'm not sure if it is going to work here.
  #    ncvar <- ncvar_def(name=varid(x)[1],dim=list(dimT,dimS), units=ifelse(unitx[1]=="\u00B0C", "degC",unitx[1]),
  #                       longname=attr(x,'longname')[1], prec=prec,compression=9,verbose=verbose)
  if (verbose) print(paste('write2ncdf4.station: longname=',attr(x,'longname')))
  options(encoding='UTF-8')
  longname <- attr(x,'longname')[1]
  if ( (is.null(longname)) | !is.character(longname) ) { 
    longname <- switch(varid(x)[1],'t2m'='daily mean temperature', 'prec'='24-hr precipitation',
                       'slp'='daily sea-level pressure','rr'='24-hr precipitation','tg'='daily mean temperature',
                       'tmax'='daily maximum temperature','tmin'='daily minimum temperature',
                       'tx'='daily maximum temperature','tn'='daily minimum temperature',
                       'sd'='snow depth','fg'='daily mean windspeed','fx'='daily maximum windspeed',
                       'cc'='cloud cover','ss'='daily sunshine duration','qq'='global radiation',
                       'dd'='daily wind direction','pp'='daily sea-level pressure')
    if (is.null(longname)) longname <- 'NA'
  }
  if (verbose) print(paste(varid(x)[1],longname))
  ncvar <- ncvar_def(name=varid(x)[1],dim=list(dimS,dimT), units=ifelse(unitx[1]=="\u00B0C", "degC",unitx[1]),
                     longname=longname, prec=prec,compression=9,verbose=verbose)
  
  if (verbose) print('The variables have been defined - now the summary statistics...')
  fyrid <- ncvar_def(name="first",dim=list(dimS), units="year", missval=missval,longname="first_year", 
                     prec="short",verbose=verbose)
  lyrid <- ncvar_def(name="last",dim=list(dimS), units="year", missval=missval,longname="last_year", 
                     prec="short",verbose=verbose)
  doid <- ncvar_def(name="days_old",dim=list(dimS), units="day", missval=missval,longname="days since last record and saving date", 
                    prec="integer",verbose=verbose)
  nvid <- ncvar_def(name="number",dim=list(dimS), units="count", missval=missval,longname="number_valid_data", 
                    prec="float",verbose=verbose)
  
  ## KMP 2018-11-02: devtools (run_examples) can only handle ASCII characters so I had to replace the 
  ## degree symbol with "\u00B0", but I'm not sure if it is going to work here.
  maxid <- ncvar_def(name="summary_max",dim=list(dimS), units= ifelse(attr(x,"unit")[1]=="\u00B0C", "degC",attr(x,"unit")[1]),
                     missval=missval,longname=varid(x)[1],prec="float",verbose=verbose)
  minid <- ncvar_def(name="summary_min",dim=list(dimS), ifelse(attr(x,"unit")[1]=="\u00B0C", "degC",attr(x,"unit")[1]), 
                     missval=missval,longname=varid(x)[1],prec="float",verbose=verbose)
  nhrid <- ncvar_def(name="summary_records",dim=list(dimS),
                     units=ifelse(attr(x,"unit")[1]=="\u00B0C", "degC",attr(x,"unit")[1]), 
                     missval=missval,longname="fraction_of_high_records",prec="float",verbose=verbose)
  ## REB 2024-09-24: check the unit for last_element (TRUE/FALSE)
  lehrid <- ncvar_def(name="last_element_highest",dim=list(dimS),
                      units=ifelse(attr(x,"unit")[1]=="\u00B0C", "degC",attr(x,"unit")[1]), 
                      missval=missval,longname="If_last_element_is_a_record",prec="short",verbose=verbose)
  
  if (is.T(x)) {
    meanid <- ncvar_def(name="summary_mean",dim=list(dimS), units="degC", 
                        missval=missval,longname="annual_mean_temperature",prec="float",verbose=verbose)
    meanid.djf <- ncvar_def(name="summary_mean_DJF",dim=list(dimS), units="degC", 
                            missval=missval,longname="seasonal_mean_temperature_Dec-Feb",prec="float",verbose=verbose)
    meanid.mam <- ncvar_def(name="summary_mean_MAM",dim=list(dimS), units="degC", 
                            missval=missval,longname="seasonal_mean_temperature_Mar-May",prec="float",verbose=verbose)
    meanid.jja <- ncvar_def(name="summary_mean_JJA",dim=list(dimS), units="degC", 
                            missval=missval,longname="seasonal_mean_temperature_Jun-Aug",prec="float",verbose=verbose)
    meanid.son <- ncvar_def(name="summary_mean_SON",dim=list(dimS), units="degC", 
                            missval=missval,longname="seasonal_mean_temperature_Sep-Nov",prec="float",verbose=verbose)
    sdid <- ncvar_def(name="summary_sd",dim=list(dimS), units="degC", 
                      missval=missval,longname="annual_temperature_anomaly_standard_deviation",prec="float",verbose=verbose)
    sdid.djf <- ncvar_def(name="summary_sd_DJF",dim=list(dimS), units="degC", 
                          missval=missval,longname="temperature_anomaly_standard_deviation_Dec-Feb",prec="float",verbose=verbose)
    sdid.mam <- ncvar_def(name="summary_sd_MAM",dim=list(dimS), units="degC", 
                          missval=missval,longname="temperature_anomaly_standard_deviation_Mar-May",prec="float",verbose=verbose)
    sdid.jja <- ncvar_def(name="summary_sd_JJA",dim=list(dimS), units="degC", 
                          missval=missval,longname="temperature_anomaly_standard_deviation_Jun-Aug",prec="float",verbose=verbose)
    sdid.son <- ncvar_def(name="summary_sd_SON",dim=list(dimS), units="degC", 
                          missval=missval,longname="temperature_anomaly_standard_deviation_Sep-Nov",prec="float",verbose=verbose)
    nlrid <- ncvar_def(name="summary_lows",dim=list(dimS), units="degC", 
                       missval=missval,longname="fraction_of_low_records",prec="float",verbose=verbose)
    tdid <- ncvar_def(name="summary_trend",dim=list(dimS), units="degC/decade", 
                      missval=missval,longname="annual_mean_temperature",prec="float",verbose=verbose)
    tdid.djf <- ncvar_def(name="summary_trend_DJF",dim=list(dimS), units="degC/decade", 
                          missval=missval,longname="seasonal_mean_temperature_Dec-Feb",prec="float",verbose=verbose)
    tdid.mam <- ncvar_def(name="summary_trend_MAM",dim=list(dimS), units="degC/decade", 
                          missval=missval,longname="seasonal_mean_temperature_Mar-May",prec="float",verbose=verbose)
    tdid.jja <- ncvar_def(name="summary_trend_JJA",dim=list(dimS), units="degC/decade", 
                          missval=missval,longname="seasonal_mean_temperature_Jun-Aug",prec="float",verbose=verbose)
    tdid.son <- ncvar_def(name="summary_trend_SON",dim=list(dimS), units="degC/decade", 
                          missval=missval,longname="seasonal_mean_temperature_Sep-Nov",prec="float",verbose=verbose)
    ## KMP 2018-11-02: devtools (run_examples) can only handle ASCII characters so I had to replace the 
    ## degree symbol with "\u00B0", but I'm not sure if it is going to work here.
    lelrid <- ncvar_def(name="last_element_lowest",dim=list(dimS), 
                        units=ifelse(attr(x,"unit")[1]=="\u00B0C", "degC",attr(x,"unit")[1]), 
                        missval=missval,longname="If_last_element_is_a_record",prec="short",verbose=verbose)
    
  } else if (is.precip(x)) {
    meanid <- ncvar_def(name="summary_mean",dim=list(dimS), units="mm/year", 
                        missval=missval,longname="mean_annual_precipitation_sum",prec="float",verbose=verbose)
    meanid.djf <- ncvar_def(name="summary_mean_DJF",dim=list(dimS), units="mm/season", 
                            missval=missval,longname="mean_seasonal_precip_sum_Dec-Feb",prec="float",verbose=verbose)
    meanid.mam <- ncvar_def(name="summary_mean_MAM",dim=list(dimS), units="mm/season", 
                            missval=missval,longname="mean_seasonal_precip_sum_Mar-May",prec="float",verbose=verbose)
    meanid.jja <- ncvar_def(name="summary_mean_JJA",dim=list(dimS), units="mm/season", 
                            missval=missval,longname="mean_seasonal_precip_sum_Jun-Aug",prec="float",verbose=verbose)
    meanid.son <- ncvar_def(name="summary_mean_SON",dim=list(dimS), units="mm/season", 
                            missval=missval,longname="mean_seasonal_precip_sum_Sep-Nov",prec="float",verbose=verbose)
    muid <- ncvar_def(name="summary_wetmean",dim=list(dimS), units="mm/day", 
                      missval=missval,longname="mean_annual_precipitation_wet_day_mean",prec="float",verbose=verbose)
    muid.djf <- ncvar_def(name="summary_wetmean_DJF",dim=list(dimS), units="mm/day", 
                          missval=missval,longname="mean_seasonal_precip_wetmean_Dec-Feb",prec="float",verbose=verbose)
    muid.mam <- ncvar_def(name="summary_wetmean_MAM",dim=list(dimS), units="mm/day", 
                          missval=missval,longname="mean_seasonal_precip_wetmean_Mar-May",prec="float",verbose=verbose)
    muid.jja <- ncvar_def(name="summary_wetmean_JJA",dim=list(dimS), units="mm/day", 
                          missval=missval,longname="mean_seasonal_precip_wetmean_Jun-Aug",prec="float",verbose=verbose)
    muid.son <- ncvar_def(name="summary_wetmean_SON",dim=list(dimS), units="mm/day", 
                          missval=missval,longname="mean_seasonal_precip_wetmean_Sep-Nov",prec="float",verbose=verbose)
    fwid <- ncvar_def(name="summary_wetfreq",dim=list(dimS), units="mm/day", 
                      missval=missval,longname="mean_annual_precipitation_wet_day_frequency",prec="float",verbose=verbose)
    fwid.djf <- ncvar_def(name="summary_wetfreq_DJF",dim=list(dimS), units="%", 
                          missval=missval,longname="mean_seasonal_precip_wetfreq_Dec-Feb",prec="float",verbose=verbose)
    fwid.mam <- ncvar_def(name="summary_wetfreq_MAM",dim=list(dimS), units="%", 
                          missval=missval,longname="mean_seasonal_precip_wetfreq_Mar-May",prec="float",verbose=verbose)
    fwid.jja <- ncvar_def(name="summary_wetfreq_JJA",dim=list(dimS), units="%", 
                          missval=missval,longname="mean_seasonal_precip_wetfreq_Jun-Aug",prec="float",verbose=verbose)
    fwid.son <- ncvar_def(name="summary_wetfreq_SON",dim=list(dimS), units="%", 
                          missval=missval,longname="mean_seasonal_precip_wetfreq_Sep-Nov",prec="float",verbose=verbose)
    tdid <- ncvar_def(name="summary_trend",dim=list(dimS), units="annual mm/decade", 
                      missval=missval,longname="trend_in_annual_precipitation_sum",prec="float",verbose=verbose)
    tdid.djf <- ncvar_def(name="summary_trend_DJF",dim=list(dimS), units="seasonal mm/decade", 
                          missval=missval,longname="trend_in_seasonal_precip_sum_Dec-Feb",prec="float",verbose=verbose)
    tdid.mam <- ncvar_def(name="summary_trend_MAM",dim=list(dimS), units="seasonal mm/decade", 
                          missval=missval,longname="trend_in_seasonal_precip_sum_Mar-May",prec="float",verbose=verbose)
    tdid.jja <- ncvar_def(name="summary_trend_JJA",dim=list(dimS), units="seasonal mm/decade", 
                          missval=missval,longname="trend_in_seasonal_precip_sum_Jun-Aug",prec="float",verbose=verbose)
    tdid.son <- ncvar_def(name="summary_trend_SON",dim=list(dimS), units="seasonal mm/decade", 
                          missval=missval,longname="trend_in_seasonal_precip_sum_Sep-Nov",prec="float",verbose=verbose)
    tdmuid <- ncvar_def(name="summary_trend_wetmean",dim=list(dimS), units="mm/day/decade", 
                        missval=missval,longname="annual_precipitation_wetday_mean",prec="float",verbose=verbose)
    tdmuid.djf <- ncvar_def(name="summary_trend_wetmean_DJF",dim=list(dimS), units="mm/day/decade", 
                            missval=missval,longname="trend_in_seasonal_precip_wetmean_Dec-Feb",prec="float",verbose=verbose)
    tdmuid.mam <- ncvar_def(name="summary_trend_wetmean_MAM",dim=list(dimS), units="mm/day/decade", 
                            missval=missval,longname="trend_in_seasonal_precip_wetmean_Mar-May",prec="float",verbose=verbose)
    tdmuid.jja <- ncvar_def(name="summary_trend_wetmean_JJA",dim=list(dimS), units="mm/day/decade", 
                            missval=missval,longname="trend_in_seasonal_precip_wetmean_Jun-Aug",prec="float",verbose=verbose)
    tdmuid.son <- ncvar_def(name="summary_trend_wetmean_SON",dim=list(dimS), units="mm/day/decade", 
                            missval=missval,longname="trend_in_seasonal_precip_wetmean_Sep-Nov",prec="float",verbose=verbose)
    tdfwid <- ncvar_def(name="summary_trend_wetfreq",dim=list(dimS), units="%/decade", 
                        missval=missval,longname="annual_precipitation_wet_day_frequency",prec="float",verbose=verbose)
    tdfwid.djf <- ncvar_def(name="summary_trend_wetfreq_DJF",dim=list(dimS), units="%/decade", 
                            missval=missval,longname="trend_in_seasonal_precip_wet_frequency_Dec-Feb",prec="float",verbose=verbose)
    tdfwid.mam <- ncvar_def(name="summary_trend_wetfreq_MAM",dim=list(dimS), units="%/decade", 
                            missval=missval,longname="trend_in_seasonal_precip_wet_frequency_Mar-May",prec="float",verbose=verbose)
    tdfwid.jja <- ncvar_def(name="summary_trend_wetfreq_JJA",dim=list(dimS), units="%/decade", 
                            missval=missval,longname="trend_in_seasonal_precip_wet_frequency_Jun-Aug",prec="float",verbose=verbose)
    tdfwid.son <- ncvar_def(name="summary_trend_wetfreq_SON",dim=list(dimS), units="%/decade", 
                            missval=missval,longname="trend_in_seasonal_precip_wet_frequency_Sep-Nov",prec="float",verbose=verbose)
    lrid <- ncvar_def(name="summary_lastrains",dim=list(dimS), units="days", 
                      missval=missval,longname="number_of_dry_days_at_end_of_record",prec="float",verbose=verbose)
    ldid <- ncvar_def(name="summary_lastdry",dim=list(dimS), units="days", 
                      missval=missval,longname="number_of_wet_days_at_end_of_record",prec="float",verbose=verbose)
    sigma2id <- ncvar_def(name="summary_sigma2",dim=list(dimS), units="mm^2", 
                          missval=missval,longname="variance_daily_precip",prec="float",verbose=verbose)
    sigma2id.djf <- ncvar_def(name="summary_sigma2_DJF",dim=list(dimS), units="mm^2", 
                              missval=missval,longname="variance_daily_precip_Dec-Feb",prec="float",verbose=verbose)
    sigma2id.mam <- ncvar_def(name="summary_sigma2_MAM",dim=list(dimS), units="mm^2", 
                              missval=missval,longname="variance_daily_precip_Mar-May",prec="float",verbose=verbose)
    sigma2id.jja <- ncvar_def(name="summary_sigma2_JJA",dim=list(dimS), units="mm^2", 
                              missval=missval,longname="variance_daily_precip_Jun-Aug",prec="float",verbose=verbose)
    sigma2id.son <- ncvar_def(name="summary_sigma2_SON",dim=list(dimS), units="mm^2", 
                              missval=missval,longname="variance_daily_precip_Sep-Nov",prec="float",verbose=verbose)
    tsigma2id <- ncvar_def(name="summary_trend_sigma2",dim=list(dimS), units="mm^2", 
                           missval=missval,longname="variance_daily_precip_trend",prec="float",verbose=verbose)
    tsigma2id.djf <- ncvar_def(name="summary_trend_sigma2_DJF",dim=list(dimS), units="mm^2", 
                               missval=missval,longname="variance_daily_precip_trend_Dec-Feb",prec="float",verbose=verbose)
    tsigma2id.mam <- ncvar_def(name="summary_trend_sigma2_MAM",dim=list(dimS), units="mm^2", 
                               missval=missval,longname="variance_daily_precip_trend_Mar-May",prec="float",verbose=verbose)
    tsigma2id.jja <- ncvar_def(name="summary_trend_sigma2_JJA",dim=list(dimS), units="mm^2", 
                               missval=missval,longname="variance_daily_precip_trend_Jun-Aug",prec="float",verbose=verbose)
    tsigma2id.son <- ncvar_def(name="summary_trend_sigma2_SON",dim=list(dimS), units="mm^2", 
                               missval=missval,longname="variance_daily_precip_trend_Sep-Nov",prec="float",verbose=verbose)
    mwslid <- ncvar_def(name="summary_mean_wetdur",dim=list(dimS), units="day", 
                        missval=missval,longname="mean_wet-day-spell_length",prec="float",verbose=verbose)
    mdslid <- ncvar_def(name="summary_mean_drydur",dim=list(dimS), units="day", 
                        missval=missval,longname="mean_dry-spell_length",prec="float",verbose=verbose)
  } else {
    meanid <- ncvar_def(name="summary_mean",dim=list(dimS), units=attr(x,"unit")[1], 
                        missval=missval,longname=paste("mean_annual",varid(x)[1],sep='_'),prec="float",verbose=verbose)
    meanid.djf <- ncvar_def(name="summary_mean_DJF",dim=list(dimS), units=attr(x,"unit")[1], 
                            missval=missval,longname=paste("Dec-Feb_mean",varid(x)[1],sep='_'),prec="float",verbose=verbose)
    meanid.mam <- ncvar_def(name="summary_mean_MAM",dim=list(dimS), units=attr(x,"unit")[1], 
                            missval=missval,longname=paste("Mar-May_mean",varid(x)[1],sep='_'),prec="float",verbose=verbose)
    meanid.jja <- ncvar_def(name="summary_mean_JJA",dim=list(dimS), units=attr(x,"unit")[1], 
                            missval=missval,longname=paste("Jun-Aug_mean",varid(x)[1],sep='_'),prec="float",verbose=verbose)
    meanid.son <- ncvar_def(name="summary_mean_SON",dim=list(dimS), units=attr(x,"unit")[1], 
                            missval=missval,longname=paste("Sep-Nov_mean",varid(x)[1],sep='_'),prec="float",verbose=verbose)
    sdid <- ncvar_def(name="summary_sd",dim=list(dimS), units="degC", 
                      missval=missval,longname="annual_temperature_anomaly_standard_deviation",prec="float",verbose=verbose)
    sdid.djf <- ncvar_def(name="summary_sd_DJF",dim=list(dimS), units="degC", 
                          missval=missval,longname="temperature_anomaly_standard_deviation_Dec-Feb",prec="float",verbose=verbose)
    sdid.mam <- ncvar_def(name="summary_sd_MAM",dim=list(dimS), units="degC", 
                          missval=missval,longname="temperature_anomaly_standard_deviation_Mar-May",prec="float",verbose=verbose)
    sdid.jja <- ncvar_def(name="summary_sd_JJA",dim=list(dimS), units="degC", 
                          missval=missval,longname="temperature_anomaly_standard_deviation_Jun-Aug",prec="float",verbose=verbose)
    sdid.son <- ncvar_def(name="summary_sd_SON",dim=list(dimS), units="degC", 
                          missval=missval,longname="temperature_anomaly_standard_deviation_Sep-Nov",prec="float",verbose=verbose)
    tdid <- ncvar_def(name="summary_trend",dim=list(dimS), units=attr(x,"unit")[1], 
                      missval=missval,longname="annual_mean_temperature",prec="float",verbose=verbose)
    tdid.djf <- ncvar_def(name="summary_trend_DJF",dim=list(dimS), units=attr(x,"unit")[1], 
                          missval=missval,longname=paste("Dec-Feb_mean",varid(x)[1],sep='_'),prec="float",verbose=verbose)
    tdid.mam <- ncvar_def(name="summary_trend_MAM",dim=list(dimS), units=attr(x,"unit")[1], 
                          missval=missval,longname=paste("Mar-May_mean",varid(x)[1],sep='_'),prec="float",verbose=verbose)
    tdid.jja <- ncvar_def(name="summary_trend_JJA",dim=list(dimS), units=attr(x,"unit")[1], 
                          missval=missval,longname=paste("Jun-Aug_mean",varid(x)[1],sep='_'),prec="float",verbose=verbose)
    tdid.son <- ncvar_def(name="summary_trend_SON",dim=list(dimS), units=attr(x,"unit")[1], 
                          missval=missval,longname=paste("Sep-Nov_mean",varid(x)[1],sep='_'),prec="float",verbose=verbose)
    ## KMP 2018-11-02: devtools (run_examples) can only handle ASCII characters so I had to replace the 
    ## degree symbol with "\u00B0", but I'm not sure if it is going to work here.
    lelrid <- ncvar_def(name="last_element_lowest",dim=list(dimS), 
                        units=ifelse(attr(x,"unit")[1]=="\u00B0C", "degC",attr(x,"unit")[1]), 
                        missval=missval,longname="If_last_element_is_a_record",prec="short",verbose=verbose)
    nlrid <- ncvar_def(name="summary_lows",dim=list(dimS), units="degC", 
                       missval=missval,longname="fraction_of_low_records",prec="float",verbose=verbose)
  }
  
  ## Set start and count: the time period it defined is already defined in y (padded NAs) 
  start <- c(1, 1)
  if (!is.null(dim(y))) count <- dim(t(y)) else count <- c(dim(y)[2],length(index(y)))
  
  if (verbose) {print(paste('Creating file',file)); str(ncvar)}
  if (is.T(x)) {
    ncid <- nc_create(file,vars=list(ncvar,lonid,latid,altid,locid,stid,cntrid, 
                                     fyrid,lyrid,nvid,meanid,meanid.djf,meanid.mam,meanid.jja,meanid.son,
                                     sdid,sdid.djf,sdid.mam,sdid.jja,sdid.son,maxid,minid,nhrid,nlrid,
                                     tdid,tdid.djf,tdid.mam,tdid.jja,tdid.son,lehrid,lelrid,doid))
  } else if (is.precip(x)) {
    ncid <- nc_create(file,vars=list(ncvar,lonid,latid,altid,locid,stid,cntrid, 
                                     fyrid,lyrid,nvid,meanid,meanid.djf,meanid.mam,meanid.jja,meanid.son,
                                     maxid,minid,nhrid,muid,muid.djf,muid.mam,muid.jja,muid.son,
                                     fwid,fwid.djf,fwid.mam,fwid.jja,fwid.son,
                                     tdid,tdid.djf,tdid.mam,tdid.jja,tdid.son,
                                     tdmuid,tdmuid.djf,tdmuid.mam,tdmuid.jja,tdmuid.son,
                                     tdfwid,tdfwid.djf,tdfwid.mam,tdfwid.jja,tdfwid.son,lrid,ldid,lehrid,
                                     sigma2id,sigma2id.djf,sigma2id.mam,sigma2id.jja,sigma2id.son,
                                     tsigma2id,tsigma2id.djf,tsigma2id.mam,tsigma2id.jja,tsigma2id.son,
                                     mwslid,mdslid,doid))
  } else {
    ncid <- nc_create(file,vars=list(ncvar,lonid,latid,altid,locid,stid,cntrid, 
                                     fyrid,lyrid,nvid,meanid,meanid.djf,meanid.mam,meanid.jja,meanid.son,
                                     sdid,sdid.djf,sdid.mam,sdid.jja,sdid.son,tdid,
                                     tdid.djf,tdid.mam,tdid.jja,tdid.son,maxid,minid,nhrid,nlrid,lehrid,
                                     lelrid,doid))
  }
  #}
  
  if (verbose) {
    print("start & count"); print(start); print(count); 
    print("dim(y)"); print(dim(y)); print(c(nt,length(y)))
    print("time"); print(range(time))
    print("netCDF dimensions"); print(c(nt,ns)); print(start+count-c(1,1))
    if (count[2]==0) browser('Detected suspect situation: count[0] = 0.')
  }
  
  if (verbose) {print('Saving the main data'); print(start); print(count); print(summary(c(y))); print(dim(t(y)))}
  ## Store the values -coredata - temporarily in yc for writing to netCDF
  yc <- coredata(y); bad <- !is.finite(yc)
  yc[bad] <- missval
  if (verbose) print(summary(c(yc)))
  
  ncvar_put( ncid, ncvar, t(yc),start=start,count=count); rm('yc')
  if (verbose) print('Saving attributes')
  ncatt_put( ncid, ncvar, 'add_offset',offset,prec='float')
  ncatt_put( ncid, ncvar, 'scale_factor',scale,prec='float')
  ncatt_put( ncid, ncvar, 'missing_value',missval,prec='float')
  if (verbose) print(paste('Saving metadata:',start[1],count[1]))
  ncvar_put( ncid, lonid, lon(x),start=start[1],count=count[1])
  ncvar_put( ncid, latid, lat(x),start=start[1],count=count[1])
  ncvar_put( ncid, altid, alt(x),start=start[1],count=count[1])
  if (verbose) print('First & last year')
  ncvar_put( ncid, fyrid, firstyear(x),start=start[1],count=count[1])
  ncvar_put( ncid, lyrid, lastyear(x),start=start[1],count=count[1])
  ncvar_put( ncid, doid, daysold,start=start[1],count=count[1])
  if (is.null(dim(x))) number <- sum(is.finite(coredata(x))) else
    if (length(dim(x))==2) number <- apply(coredata(x),2,FUN='nv') else number <- -1
  #ncvar_put( ncid, nvid, number,start=start[2],count=count[2])
  ncvar_put( ncid, nvid, number,start=start[1],count=count[1])
  
  if (verbose) print('Add summary statistics: mean')
  ave[insufficient] <- missval; ave.djf[insufficient] <- missval
  ave.mam[insufficient] <- missval; ave.jja[insufficient] <- missval
  ave.son[insufficient] <- missval; td[insufficient] <- missval
  td.djf[insufficient] <- missval; td.mam[insufficient] <- missval
  td.jja[insufficient] <- missval; td.son[insufficient] <- missval
  if (length(grep('mn',substr(ls(),1,2)))==0) mn <- rep(missval,length(mx)) ## REB: For precipitation, min (=0) is often not saved.
  mx[insufficient] <- missval; mn[insufficient] <- missval; 
  if (length(grep('lehr',substr(ls(),1,2)))==0) lehr <- rep(missval,length(mx)) ## REB: For precipitation, min (=0) is often not saved.
  nhr[insufficient] <- missval; lehr[insufficient] <- missval
  
  ncvar_put( ncid, meanid, ave,start=start[1],count=count[1])
  ncvar_put( ncid, meanid.djf, ave.djf,start=start[1],count=count[1])
  ncvar_put( ncid, meanid.mam, ave.mam,start=start[1],count=count[1])
  ncvar_put( ncid, meanid.jja, ave.jja,start=start[1],count=count[1])
  ncvar_put( ncid, meanid.son, ave.son,start=start[1],count=count[1])
  if (verbose) print('Add summary statistics: trend')
  ncvar_put( ncid, tdid, td  ,start=start[1],count=count[1])
  ncvar_put( ncid, tdid.djf, td.djf,start=start[1],count=count[1])
  ncvar_put( ncid, tdid.mam, td.mam,start=start[1],count=count[1])
  ncvar_put( ncid, tdid.jja, td.jja,start=start[1],count=count[1])
  ncvar_put( ncid, tdid.son, td.son,start=start[1],count=count[1])
  if (verbose) print('Add summary statistics: max, min')
  ncvar_put( ncid, maxid, mx, start=start[1],count=count[1])
  ncvar_put( ncid, minid, mn, start=start[1],count=count[1])
  if (verbose) print('Add summary statistics: records')
  ncvar_put( ncid, nhrid, nhr, start=start[1],count=count[1])
  ncvar_put( ncid, lehrid, lehr, start=start[1],count=count[1])
  
  if (is.T(x)) {
    if (verbose) print('extra for temperature')
    std[insufficient] <- missval; std.djf[insufficient] <- missval
    std.mam[insufficient] <- missval; std.jja[insufficient] <- missval
    std.son[insufficient] <- missval; nlr[insufficient] <- missval
    if (length(grep('lelr',substr(ls(),1,2)))==0) lelr <- rep(missval,length(mx)) ## REB: For precipitation, min (=0) is often not saved.
    lelr[insufficient] <- missval; 
    ncvar_put( ncid, sdid, std, start=start[1],count=count[1])
    ncvar_put( ncid, sdid.djf, std.djf, start=start[1],count=count[1])
    ncvar_put( ncid, sdid.mam, std.mam, start=start[1],count=count[1])
    ncvar_put( ncid, sdid.jja, std.jja, start=start[1],count=count[1])
    ncvar_put( ncid, sdid.son, std.son, start=start[1],count=count[1])
    ncvar_put( ncid, nlrid, nlr, start=start[1],count=count[1])
    ncvar_put( ncid, lelrid, lelr, start=start[1],count=count[1])
  }
  if (is.precip(x)) {
    if (verbose) print('extra for precipitation')
    mu[insufficient] <- missval; mu.djf[insufficient] <- missval
    mu.mam[insufficient] <- missval; mu.jja[insufficient] <- missval
    mu.son[insufficient] <- missval; fw[insufficient] <- missval
    fw.mam[insufficient] <- missval; fw.djf[insufficient] <- missval
    fw.jja[insufficient] <- missval; fw.son[insufficient] <- missval
    tdmu[insufficient] <- missval; tdmu.djf[insufficient] <- missval
    tdmu.mam[insufficient] <- missval; tdmu.jja[insufficient] <- missval
    tdmu.son[insufficient] <- missval; tdfw[insufficient] <- missval
    tdfw.mam[insufficient] <- missval; tdfw.djf[insufficient] <- missval
    tdfw.jja[insufficient] <- missval; tdfw.son[insufficient] <- missval
    lr[insufficient] <- missval; ld[insufficient] <- missval
    ncvar_put( ncid, muid, mu,start=start[1],count=count[1])
    ncvar_put( ncid, muid.djf, mu.djf,start=start[1],count=count[1])
    ncvar_put( ncid, muid.mam, mu.mam,start=start[1],count=count[1])
    ncvar_put( ncid, muid.jja, mu.jja,start=start[1],count=count[1])
    ncvar_put( ncid, muid.son, mu.son,start=start[1],count=count[1])
    ncvar_put( ncid, fwid, fw,start=start[1],count=count[1])
    ncvar_put( ncid, fwid.djf, fw.djf,start=start[1],count=count[1])
    ncvar_put( ncid, fwid.mam, fw.mam,start=start[1],count=count[1])
    ncvar_put( ncid, fwid.jja, fw.jja,start=start[1],count=count[1])
    ncvar_put( ncid, fwid.son, fw.son,start=start[1],count=count[1])
    ncvar_put( ncid, tdfwid, tdfw,start=start[1],count=count[1])
    ncvar_put( ncid, tdfwid.djf, tdfw.djf,start=start[1],count=count[1])
    ncvar_put( ncid, tdfwid.mam, tdfw.mam,start=start[1],count=count[1])
    ncvar_put( ncid, tdfwid.jja, tdfw.jja,start=start[1],count=count[1])
    ncvar_put( ncid, tdfwid.son, tdfw.son,start=start[1],count=count[1])
    ncvar_put( ncid, tdmuid, tdmu,start=start[1],count=count[1])
    ncvar_put( ncid, tdmuid.djf, tdmu.djf,start=start[1],count=count[1])
    ncvar_put( ncid, tdmuid.mam, tdmu.mam,start=start[1],count=count[1])
    ncvar_put( ncid, tdmuid.jja, tdmu.jja,start=start[1],count=count[1])
    ncvar_put( ncid, tdmuid.son, tdmu.son,start=start[1],count=count[1])
    ncvar_put( ncid, lrid, lr, start=start[1],count=count[1])
    ncvar_put( ncid, ldid, ld, start=start[1],count=count[1])
    if (verbose) print(paste('sigma2:',length(sigma2),start[1],count[1]))
    if (length(tsigma2)==count[1]) { 
      ncvar_put( ncid, sigma2id, sigma2,start=start[1],count=count[1])
      ncvar_put( ncid, sigma2id.djf, sigma2.djf,start=start[1],count=count[1])
      ncvar_put( ncid, sigma2id.mam, sigma2.mam,start=start[1],count=count[1])
      ncvar_put( ncid, sigma2id.jja, sigma2.jja,start=start[1],count=count[1])
      ncvar_put( ncid, sigma2id.son, sigma2.son,start=start[1],count=count[1])
      ncvar_put( ncid, tsigma2id, tsigma2,start=start[1],count=count[1])
      ncvar_put( ncid, tsigma2id.djf, tsigma2.djf,start=start[1],count=count[1])
      ncvar_put( ncid, tsigma2id.mam, tsigma2.mam,start=start[1],count=count[1])
      ncvar_put( ncid, tsigma2id.jja, tsigma2.jja,start=start[1],count=count[1])
      ncvar_put( ncid, tsigma2id.son, tsigma2.son,start=start[1],count=count[1])
    } else {
      print('write2ncdf4.station - precip: detected something wrong!'); 
      print(paste('sigma2:',length(sigma2),start[1],count[1]))
    }
    if (verbose) print('Mean spell length')
    ncvar_put( ncid, mwslid, mwsl,start=start[1],count=count[1])
    ncvar_put( ncid, mdslid, mdsl,start=start[1],count=count[1])
  } else {
    if (verbose) print(paste('extra for',varid(x)[1]))
    std[insufficient] <- missval; std.djf[insufficient] <- missval
    std.mam[insufficient] <- missval; std.jja[insufficient] <- missval
    std.son[insufficient] <- missval; nlr[insufficient] <- missval
    lelr[insufficient] <- missval; 
    ncvar_put( ncid, sdid, std, start=start[1],count=count[1])
    ncvar_put( ncid, sdid.djf, std.djf, start=start[1],count=count[1])
    ncvar_put( ncid, sdid.mam, std.mam, start=start[1],count=count[1])
    ncvar_put( ncid, sdid.jja, std.jja, start=start[1],count=count[1])
    ncvar_put( ncid, sdid.son, std.son, start=start[1],count=count[1])
    ncvar_put( ncid, nlrid, nlr, start=start[1],count=count[1])
    ncvar_put( ncid, lelrid, lelr, start=start[1],count=count[1])
  }
  
  ## There are some times problems saving text data, and there seems to be some 
  ## inconsistency in the ncdf4 package. To by-pass this problem, we had to make
  ## the following code more complicated. There seems to be a mix-up between the 
  ## dimensions sometimes.
  if (verbose) {print('Saving textual information'); print(dim(x)); print(dim(y))}
  test <- try(ncvar_put( ncid, locid, loc(x),start=c(1,start[1]),count=c(namelength,count[1])))
  if (inherits(test,'try-error'))
    try(ncvar_put( ncid, locid, loc(x),start=c(start[1],1),count=c(count[1],namelength)))
  if (verbose) {print('as.character(stid(x)):'); print(as.character(stid(x)))}
  test <- try(ncvar_put( ncid, stid, as.character(stid(x)),c(1,start[1]),count=c(namelength,count[1])))
  if (inherits(test,'try-error'))
    try(ncvar_put( ncid, stid, as.character(stid(x)),c(start[1],1),count=c(count[1],namelength)))
  test <- try(ncvar_put( ncid, cntrid, cntr(x),start=c(1,start[1]),count=c(namelength,count[1])))
  if (inherits(test,'try-error'))
    try(ncvar_put( ncid, cntrid, cntr(x),start=c(start[1],1),count=c(count[1],namelength)))
  #try(ncvar_put( ncid, cntrid, cntr(y),start=c(start[2],1),count=c(count[2],namelength)))
  if (verbose) print('textual data saved')
  
  #if (!append) {
  ## global attributes
  if (verbose) { print('global attributes'); print(paste(levels(factor(attr(x,"source"))),collapse="/")) }
  ncatt_put( ncid, 0, 'class', class(x))
  ncatt_put( ncid, 0, 'title', paste(levels(factor(attr(x,"info"))),collapse="/"))
  ncatt_put( ncid, 0, 'source', src)
  ncatt_put( ncid, 0, 'history', paste(unlist(attr(x,"history")),collapse="/"))
  ncatt_put( ncid, 0, 'references', paste(levels(factor(attr(x,"reference"))),collapse="/"))
  ncatt_put( ncid, 0, "esd-version", attr(x,'history')$session$esd.version)
  ## Add global attributes suggested in https://adc.met.no/node/4
  ncatt_put( ncid, 0, 'id', doi)
  ncatt_put( ncid, 0, 'naming_authority', namingauthority)
  ncatt_put( ncid, 0, 'title', 'netCDF-file for station data')
  ncatt_put( ncid, 0, 'summary', 'Daily data organised as station ID number and time')
  ncatt_put( ncid, 0, 'keywords', 'Observational record, summary statistics, station data')
  ncatt_put( ncid, 0, 'geospatial_lat_min', min(lat(x)))
  ncatt_put( ncid, 0, 'geospatial_lat_max', max(lat(x)))
  ncatt_put( ncid, 0, 'geospatial_lon_min', min(lon(x)))
  ncatt_put( ncid, 0, 'geospatial_lon_max', max(lon(x)))
  ncatt_put( ncid, 0, 'time_coverage_start', as.character(min(index(x))))
  ncatt_put( ncid, 0, 'time_coverage_end', as.character(max(index(x))))
  ncatt_put( ncid, 0, 'Conventions', 'CF-inspired')
  ncatt_put( ncid, 0, 'processing_level', processinglevel)
  ncatt_put( ncid, 0, 'date_created', date())
  ncatt_put( ncid, 0, 'creator_type', creatortype)
  ncatt_put( ncid, 0, 'creator_email', creatoremail)
  ncatt_put( ncid, 0, 'creator_url', 'github.com/metno/esd')
  ncatt_put( ncid, 0, 'institution', institution)
  ncatt_put( ncid, 0, 'publisher_name', publishername)
  ncatt_put( ncid, 0, 'publisher_email', publisheremail)
  ncatt_put( ncid, 0, 'publisher_url', publisherurl)
  ncatt_put( ncid, 0, 'project', project)
  #}
  nc_close(ncid)
  if (verbose) print('close')
}

#' @exportS3Method
#' @export write2ncdf4.eof
write2ncdf4.array <- function(x,file,...,verbose=FALSE){
  if (verbose) cat('write2ncdf4.array \n')
  args <- list(...)
  
  ## Get the attributes
  lons <- attr(x, "longitude")
  lats <- attr(x, "latitude")
  var_name <- attr(x, "variable")
  var_unit <- attr(x, "unit")
  
  
  
  dim_lon <- ncdim_def("lon", "degrees_east", as.double(lons))
  dim_lat <- ncdim_def("lat", "degrees_north", as.double(lats))
  if (length(dim(x))==3) { 
    times <- 1:4
    season_names <- "1: DJF, 2: MAM, 3: JJA, 4: SON"
    dim_time <- ncdim_def("time", "season_index", as.double(times), 
                          longname = "Season (1:DJF, 2:MAM, 3:JJA, 4:SON)")
    dims <- list(dim_lon, dim_lat, dim_time)
  } else dims <- list(dim_lon, dim_lat)
  
  
  fill_value <- -32767
  var_def <- ncvar_def(
    name = var_name,
    units = var_unit,
    dim = dims,
    missval = fill_value,
    longname = paste("Seasonal", var_name),
    prec = "float"
  )
  
  # 5. Opprett filen og skriv data
  nc_out <- nc_create(file, var_def, force_v4 = TRUE)
  
  # Skriv selve data-arrayen (X må ha dimensjonene n x m x t)
  ncvar_put(nc_out, var_def, x)
  ncatt_put(nc_out, 0, "title", "Seasonal Climate Data")
  if (length(dim(x))==3) ncatt_put(nc_out, 0, "seasons", season_names)
  ncatt_put(nc_out, 0, "history", paste("Created on", Sys.time(), "using R ncdf4"))
  
  nc_close(nc_out)
  message("Save the file as: ", file)
}
