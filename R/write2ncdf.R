## Author=? Date 

## https://www.unidata.ucar.edu/software/netcdf/docs/netcdf/CDF-Data-Types.html:
## short: 16-bit signed integers. The short type holds values between -32768 and 32767.

## Help functions 
firstyear <- function(x,na.rm=FALSE,verbose=FALSE) {
  if (verbose) print('firstyear')
  yrs <- year(x)
  if (verbose) print(range(as.numeric(yrs)))
  if (is.null(dim(x))) y <- min(yrs[is.finite(x)]) else { 
    nv <- apply(x,2,'nv')
    y <- rep(NA,length(nv))
    ok <- (1:length(nv))[nv > 0]
    y[ok] <- apply(x[,ok],2,function(x,yrs=yrs) min(yrs[is.finite(x)]),yrs)
    #for (i in ok) y[i] <- min(yrs[is.finite(x[,i])])
  }
  y[!is.finite(y)] <- NA
  if (verbose) print(table(as.numeric(y)))
  return(y)
}

lastyear <- function(x,na.rm=FALSE,verbose=FALSE) {
  if (verbose) print('lastyear')
  yrs <- year(x)
  if (verbose) print(range(as.numeric(yrs)))
  if (is.null(dim(x))) y <- max(yrs[is.finite(x)]) else { 
    nv <- apply(x,2,'nv')
    y <- rep(NA,length(nv))
    ok <- (1:length(nv))[nv > 0]
    y[ok] <- apply(x[,ok],2,function(x,yrs=yrs) max(yrs[is.finite(x)]),yrs)
  }     
  y[!is.finite(y)] <- NA
  if (verbose) print(table(as.numeric(y)))
  return(y)
}



#' Saves climate data as netCDF.
#' 
#' Method to save station data as netCDF, making sure to include the data
#' structure and meta-data (attributes). The code tries to follow the netCDf
#' 'CF' convention. The method is built on the \code{\link{ncdf4}} package.
#' 
#' 
#' @aliases write2ncdf4 write2ncdf4.station write2ncdf4.default
#' write2ncdf4.field write2ncdf4.list write2ncdf4.station write2ncdf4.eof
#' write2ncdf4.pca write2ncdf4.dsensemble
#' @param x data object
#' @param file file name
#' @param prec Precision: see \code{\link[ncdf4]{ncvar_def}}
#' @param missval Missing value: see \code{\link[ncdf4]{ncvar_def}}
#' @param offset Sets the attribute 'add_offset' which is added to the values
#' stored (to save space may be represented as 'short').
#' @param scale Sets the atttribute 'scale_factor' which is used to scale
#' (multiply) the values stored (to save space may be represented as 'short').
#' @param torg Time origin
#' @param verbose TRUE - clutter the screen.
#' @return A "zoo" "field" object with additional attributes used for further
#' processing.
#' @author R.E. Benestad
#' @seealso \code{\link{test.retrieve.ncdf4}}.
#' @keywords netcdf ncdf4 save
#' @examples
#' 
#' nacd <- station(src='nacd')
#' X <- annual(nacd)
#' write2ncdf4(X,file='test.nc')
#' 
#' @export write2ncdf4
write2ncdf4 <- function(x,...) UseMethod("write2ncdf4")

write2ncdf4.default <- function(x,...) {
}

write2ncdf4.list <- function(x,file='field.nc',prec='short',scale=0.1,offset=NULL,
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
  units <- unlist(lapply(x,function(x) attr(x,"unit")[1]))
  if (verbose) {print(varids); print(units); print(n)}
  x4nc <- list()
  for (i in 1:n) {
    if (verbose) print(paste(i,'ncvar_def',varids[i]))
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

write2ncdf4.field <- function(x,file='field.nc',prec='short',scale=NULL,offset=NULL,
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
  x4nc <- ncvar_def(varid(x)[1], attr(x,"unit")[1], list(dimlon,dimlat,dimtim), -1, 
                    longname=attr(x,'longname'), prec=prec)
     
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
  ncatt_put( ncnew, 0, 'class', class(x))
  ncatt_put( ncnew, 0, "description", 
             paste("Saved from esd using write2ncdf4",date()))
  if (verbose) print(attr(x,'history'))
  ncatt_put( ncnew, 0, "esd-version", attr(x,'history')$session$esd.version)
  nc_close(ncnew)
}



# https://www.unidata.ucar.edu/software/netcdf/docs/netcdf/CDL-Data-Types.html:
# short: 16-bit signed integers. The short type holds values between -32768 and 32767. 

write2ncdf4.station <- function(x,file='station.nc',prec='short',offset=0, missval=-99,it=NULL,stid=NULL,append=FALSE,
                                scale=0.1,torg='1899-12-31',stid_unlim=FALSE,namelength=24,verbose=FALSE) {
  
  if (!inherits(x,"station")) stop('x argument must be a station object') 
  unitx <- attr(x,'unit')

    if (verbose) {print('write2ncdf4.station'); print(range(index(x))); print(range(c(coredata(x)),na.rm=TRUE))}
  ## Quality check - remove obviously unrealistic values
  if (is.precip(x)) {
    cx <- coredata(x); cx[cx < 0] <- NA; cx[cx > 1500] <- NA; cx -> coredata(x); rm('cx') 
  }
  if (is.T(x)) {
    cx <- coredata(x); cx[cx < -100] <- NA; cx[cx > 100] <- NA; cx -> coredata(x); rm('cx')
  }
  
  cx <- coredata(x); cx[cx < -100] <- NA;
  if (prec=='short') cx[cx > 3200] <- NA;
  cx -> coredata(x); rm('cx') 
  
  ## Don't save empty space:
  x0 <- x
  if (length(dim(x))==2) {
    cx <- coredata(x); cx[cx <= missval] <- NA; coredata(x) <- cx; rm('cx')
  }
  if (sum(is.finite(x))==0) {
    print('write2ncdf4.station: No valid data after weeding missing values')
    print(summary(x0)); print(summary(x))
    return()
  }
  rm('x0'); gc(reset=TRUE)
  ## Weed out stations with short time series:
  if (length(dim(x))==2) {
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
  
  ## Get time 
  if (is.null(it)) nt <- dim(x)[1] else nt <- length(it)
  if (!is.null(stid))  {
    if (!is.null(dim(x))) nstations <- dim(x)[2] else if (!is.null(x)) nstations <- 1 else nstations <- 0
    if (append & (length(stid) != nstations)) {
      print(dim(x)); print(length(stid))
      stop('write2ncdf4.station: stid argument does not match x')
    }
    ns <- length(stid) 
  } else if (length(dim(x))==2) {
    ns <- dim(x)[2] 
    stid <- 1:ns
  } else ns <- 1
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
  
  ## For single stations, we need to fix the dimensions
  if (is.null(dim(x))) dim(x) <- c(length(x),1)
  
  ## Esimate summary statistics for the station data
  if (verbose) print('Estimate summary statistics')
  mx <- apply(x,2,'max',na.rm=TRUE)
  mn <- apply(x,2,'min',na.rm=TRUE)
  nhr <- apply(anomaly(x),2,'arec')
  nlr <- apply(-anomaly(x),2,'arec')
  
  ## Is the last element a high or low record?
  lehr <- lastelementrecord(x,verbose=verbose)
  lelr <- lastelementrecord(-x)
  if (is.T(x)) {
    if (verbose) print('Temperature')
    ## Maximum temperature
    ave <- apply(x,2,'mean',na.rm=TRUE)
    ave.djf <- apply(subset(x,it='djf'),2,'mean',na.rm=TRUE)
    ave.mam <- apply(subset(x,it='mam'),2,'mean',na.rm=TRUE)
    ave.jja <- apply(subset(x,it='jja'),2,'mean',na.rm=TRUE)
    ave.son <- apply(subset(x,it='son'),2,'mean',na.rm=TRUE)
    std <- apply(anomaly(x),2,'sd',na.rm=TRUE)
    std.djf <- apply(subset(anomaly(x),it='djf'),2,'sd',na.rm=TRUE)
    std.mam <- apply(subset(anomaly(x),it='mam'),2,'sd',na.rm=TRUE)
    std.jja <- apply(subset(anomaly(x),it='jja'),2,'sd',na.rm=TRUE)
    std.son <- apply(subset(anomaly(x),it='son'),2,'sd',na.rm=TRUE)
    td <- apply(annual(x),2,'trend.coef')
    td.djf <- apply(annual(subset(x,it='djf'),'mean',nmin=75),2,'trend.coef')
    td.mam <- apply(annual(subset(x,it='mam'),'mean',nmin=75),2,'trend.coef')
    td.jja <- apply(annual(subset(x,it='jja'),'mean',nmin=75),2,'trend.coef')
    td.son <- apply(annual(subset(x,it='son'),'mean',nmin=75),2,'trend.coef')
  } else if (is.precip(x)) {
    if (verbose) print('Precipitation')
    ave <- apply(annual(x,'sum'),2,'mean',na.rm=TRUE)
    ave.djf <- apply(annual(subset(x,it='djf'),'sum',nmin=90),2,'mean',na.rm=TRUE)
    ave.mam <- apply(annual(subset(x,it='mam'),'sum',nmin=90),2,'mean',na.rm=TRUE)
    ave.jja <- apply(annual(subset(x,it='jja'),'sum',nmin=90),2,'mean',na.rm=TRUE)
    ave.son <- apply(annual(subset(x,it='son'),'sum',nmin=90),2,'mean',na.rm=TRUE)
    mu <- apply(annual(x,'wetmean'),2,'mean',na.rm=TRUE)
    mu.djf <- apply(annual(subset(x,it='djf'),'wetmean',nmin=75),2,'mean',na.rm=TRUE)
    mu.mam <- apply(annual(subset(x,it='mam'),'wetmean',nmin=75),2,'mean',na.rm=TRUE)
    mu.jja <- apply(annual(subset(x,it='jja'),'wetmean',nmin=75),2,'mean',na.rm=TRUE)
    mu.son <- apply(annual(subset(x,it='son'),'wetmean',nmin=75),2,'mean',na.rm=TRUE)
    fw <- apply(100*annual(x,'wetfreq'),2,'mean',na.rm=TRUE)
    fw.djf <- apply(100*annual(subset(x,it='djf'),'wetfreq',nmin=75),2,'mean',na.rm=TRUE)
    fw.mam <- apply(100*annual(subset(x,it='mam'),'wetfreq',nmin=75),2,'mean',na.rm=TRUE)
    fw.jja <- apply(100*annual(subset(x,it='jja'),'wetfreq',nmin=75),2,'mean',na.rm=TRUE)
    fw.son <- apply(100*annual(subset(x,it='son'),'wetfreq',nmin=75),2,'mean',na.rm=TRUE)
    td <- apply(annual(x,FUN='sum'),2,'trend.coef')
    td.djf <- apply(annual(subset(x,it='djf'),'sum',nmin=90),2,'trend.coef')
    td.mam <- apply(annual(subset(x,it='mam'),'sum',nmin=90),2,'trend.coef')
    td.jja <- apply(annual(subset(x,it='jja'),'sum',nmin=90),2,'trend.coef')
    td.son <- apply(annual(subset(x,it='son'),'sum',nmin=90),2,'trend.coef')
    tdfw <- apply(100*annual(x,FUN='wetfreq'),2,'trend.coef')
    tdfw.djf <- apply(100*annual(subset(x,it='djf'),'wetfreq',nmin=75),2,'trend.coef')
    tdfw.mam <- apply(100*annual(subset(x,it='mam'),'wetfreq',nmin=75),2,'trend.coef')
    tdfw.jja <- apply(100*annual(subset(x,it='jja'),'wetfreq',nmin=75),2,'trend.coef')
    tdfw.son <- apply(100*annual(subset(x,it='son'),'wetfreq',nmin=75),2,'trend.coef')
    tdmu <- apply(annual(x,FUN='wetmean'),2,'trend.coef')
    tdmu.djf <- apply(annual(subset(x,it='djf'),'wetmean',nmin=75),2,'trend.coef')
    tdmu.mam <- apply(annual(subset(x,it='mam'),'wetmean',nmin=75),2,'trend.coef')
    tdmu.jja <- apply(annual(subset(x,it='jja'),'wetmean',nmin=75),2,'trend.coef')
    tdmu.son <- apply(annual(subset(x,it='son'),'wetmean',nmin=75),2,'trend.coef')
    lr <- sapply(x,'lastrains')
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
    ss <- spell(x,1)
    mwsl <- colMeans(subset.station(ss,is=list(param='wet')),na.rm=TRUE)
    mdsl <- colMeans(subset.station(ss,is=list(param='dry')),na.rm=TRUE)
  } else {
    ave <- apply(x,2,'mean',na.rm=TRUE)
    ave.djf <- apply(subset(x,it='djf'),2,'mean',na.rm=TRUE)
    ave.mam <- apply(subset(x,it='mam'),2,'mean',na.rm=TRUE)
    ave.jja <- apply(subset(x,it='jja'),2,'mean',na.rm=TRUE)
    ave.son <- apply(subset(x,it='son'),2,'mean',na.rm=TRUE)
    std <- apply(anomaly(x),2,'sd',na.rm=TRUE)
    std.djf <- apply(subset(anomaly(x),it='djf'),2,'sd',na.rm=TRUE)
    std.mam <- apply(subset(anomaly(x),it='mam'),2,'sd',na.rm=TRUE)
    std.jja <- apply(subset(anomaly(x),it='jja'),2,'sd',na.rm=TRUE)
    std.son <- apply(subset(anomaly(x),it='son'),2,'sd',na.rm=TRUE)
    td <- apply(annual(x),2,'trend.coef')
    td.djf <- apply(annual(subset(x,it='djf'),'mean',nmin=75),2,'trend.coef')
    td.mam <- apply(annual(subset(x,it='mam'),'mean',nmin=75),2,'trend.coef')
    td.jja <- apply(annual(subset(x,it='jja'),'mean',nmin=75),2,'trend.coef')
    td.son <- apply(annual(subset(x,it='son'),'mean',nmin=75),2,'trend.coef')
  }
  if (verbose) print('Summary statistics computed')
  ## Only do summary statistics for stations with more than 30 years
  insufficient <- apply(coredata(x),2,nv) < 30*365
  if (verbose) print(nv)
  
  y <- coredata(x)
  y[!is.finite(y)] <- missval
  y <- round((y - offset)/scale)
  
  if (is.null(attr(x,'calendar'))) {
    attr(x,'calendar') <- 'standard'
  }
  
  if (!is.na(attr(x,'calendar'))) {
    calendar <- attr(x,'calendar')
  } else {
    calendar <- 'standard'
  }
  
  if (class(index(x))=='Date') {
    if (is.null(it)) {
      time <- julian(index(x)) - julian(as.Date(torg)) 
    } else {
      if (verbose) print('Use prescribed time coordinates')
      y <- zoo(y,order.by=index(x))
      x2 <- merge(zoo(rep(0,nt),order.by=it),zoo(y),all=FALSE)
      x2 <- window(x2[,-1],start=it[1],end=it[length(it)])
      x2 <- attrcp(x,x2); class(x2) <- class(x); y <- x2; rm('x2')
      time <- julian(it) - julian(as.Date(torg))
      if (verbose) {
        print(range(index(x)))
        print(range(it))
        print(dim(x))
      }
    }
  } else if (inherits(x,'annual')) {
    time <- julian(as.Date(paste(year(x),'01-01',sep='-')))-julian(as.Date(torg))
  }
  if(verbose) print(paste('Period in data: ',min(firstyear(x)),' - ', max(lastyear(x)),' and time dimension: ',
                          paste(range(as.Date(time,origin=torg)),collapse=' - ')))
  
  # Attributes with same number of elements as stations are saved as variables
  if (is.null(it)) it <- index(y)
  start <- c( (1:length(it))[is.element(it,index(y)[1])],stid[1] )
  if (length(start)==1) start <- c(start,1)
  
  if (!is.null(dim(y))) count <- dim(y) else count <- c(length(y),1)
  if (verbose) {
    print("start & count"); print(start); print(count); 
    print("dim(y)"); print(dim(y))
    print("time"); print(range(time))
    print("netCDF dimensions"); print(c(nt,ns)); print(start+count-c(1,1))
  }
  
  # Define the dimensions
  if (!append) {
    if (verbose) print('Define dimensions')
    if (verbose) print(stid(x))
    dimS <- ncdim_def( name="stid", units="number",vals=1:ns,unlim=stid_unlim)
    #dimT <- ncdim_def( name="time", units=paste("days since",torg), vals=1:nt, calendar=calendar,unlim=TRUE)
    dimT <- ncdim_def( name="time", units=paste("days since",torg), vals=time, calendar=calendar,unlim=TRUE)
    dimnchar   <- ncdim_def("nchar",   "", 1:namelength, create_dimvar=FALSE )
    #dimstation <- ncdim_def("station", "", 1:ns, create_dimvar=FALSE )
    
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
    ncvar <- ncvar_def(name=varid(x)[1],dim=list(dimT,dimS), units=ifelse(unitx[1]=="\u00B0C", "degC",unitx[1]),
                       longname=attr(x,'longname')[1], prec=prec,compression=9,verbose=verbose)
    
    if (verbose) print('The variables have been defined - now the summary statistics...')
    
    fyrid <- ncvar_def(name="first",dim=list(dimS), units="year", missval=missval,longname="first_year", 
                       prec="short",verbose=verbose)
    lyrid <- ncvar_def(name="last",dim=list(dimS), units="year", missval=missval,longname="last_year", 
                       prec="short",verbose=verbose)
    nvid <- ncvar_def(name="number",dim=list(dimS), units="count", missval=missval,longname="number_valid_data", 
                      prec="float",verbose=verbose)
    
    ## KMP 2018-11-02: devtools (run_examples) can only handle ASCII characters so I had to replace the 
    ## degree symbol with "\u00B0", but I'm not sure if it is going to work here.
    maxid <- ncvar_def(name="summary_max",dim=list(dimS), units= ifelse(attr(x,"unit")[1]=="\u00B0C", "degC",attr(x,"unit")[1]),
                       missval=missval,longname=varid(x),prec="float",verbose=verbose)
    minid <- ncvar_def(name="summary_min",dim=list(dimS), ifelse(attr(x,"unit")[1]=="\u00B0C", "degC",attr(x,"unit")[1]), 
                       missval=missval,longname=varid(x),prec="float",verbose=verbose)
    nhrid <- ncvar_def(name="summary_records",dim=list(dimS),
                       units=ifelse(attr(x,"unit")[1]=="\u00B0C", "degC",attr(x,"unit")[1]), 
                       missval=missval,longname="fraction_of_high_records",prec="float",verbose=verbose)
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
                          missval=missval,longname=paste("mean_annual",varid(x),sep='_'),prec="float",verbose=verbose)
      meanid.djf <- ncvar_def(name="summary_mean_DJF",dim=list(dimS), units=attr(x,"unit")[1], 
                              missval=missval,longname=paste("Dec-Feb_mean",varid(x),sep='_'),prec="float",verbose=verbose)
      meanid.mam <- ncvar_def(name="summary_mean_MAM",dim=list(dimS), units=attr(x,"unit")[1], 
                              missval=missval,longname=paste("Mar-May_mean",varid(x),sep='_'),prec="float",verbose=verbose)
      meanid.jja <- ncvar_def(name="summary_mean_JJA",dim=list(dimS), units=attr(x,"unit")[1], 
                              missval=missval,longname=paste("Jun-Aug_mean",varid(x),sep='_'),prec="float",verbose=verbose)
      meanid.son <- ncvar_def(name="summary_mean_SON",dim=list(dimS), units=attr(x,"unit")[1], 
                              missval=missval,longname=paste("Sep-Nov_mean",varid(x),sep='_'),prec="float",verbose=verbose)
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
                            missval=missval,longname=paste("Dec-Feb_mean",varid(x),sep='_'),prec="float",verbose=verbose)
      tdid.mam <- ncvar_def(name="summary_trend_MAM",dim=list(dimS), units=attr(x,"unit")[1], 
                            missval=missval,longname=paste("Mar-May_mean",varid(x),sep='_'),prec="float",verbose=verbose)
      tdid.jja <- ncvar_def(name="summary_trend_JJA",dim=list(dimS), units=attr(x,"unit")[1], 
                            missval=missval,longname=paste("Jun-Aug_mean",varid(x),sep='_'),prec="float",verbose=verbose)
      tdid.son <- ncvar_def(name="summary_trend_SON",dim=list(dimS), units=attr(x,"unit")[1], 
                            missval=missval,longname=paste("Sep-Nov_mean",varid(x),sep='_'),prec="float",verbose=verbose)
      ## KMP 2018-11-02: devtools (run_examples) can only handle ASCII characters so I had to replace the 
      ## degree symbol with "\u00B0", but I'm not sure if it is going to work here.
      lelrid <- ncvar_def(name="last_element_lowest",dim=list(dimS), 
                          units=ifelse(attr(x,"unit")[1]=="\u00B0C", "degC",attr(x,"unit")[1]), 
                          missval=missval,longname="If_last_element_is_a_record",prec="short",verbose=verbose)
    }
  } 
  
  if (append & file.exists(file)) {
    if (verbose) print(paste('Appending',file))
    ncid <- nc_open(file, write=TRUE)
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
    meanid <- ncid$var[["summary_mean"]]
    meanid.djf <- ncid$var[["summary_mean_DJF"]]
    meanid.mam <- ncid$var[["summary_mean_MAM"]]
    meanid.jja <- ncid$var[["summary_mean_JJA"]]
    meanid.son <- ncid$var[["summary_mean_SON"]]
    tdid <- ncid$var[["summary_trend"]]
    tdid.djf <- ncid$var[["summary_trend_DJF"]]
    tdid.mam <- ncid$var[["summary_trend_MAM"]]
    tdid.jja <- ncid$var[["summary_trend_JJA"]]
    tdid.son <- ncid$var[["summary_trend_SON"]]
    maxid <- ncid$var[["summary_max"]]
    minid <- ncid$var[["summary_min"]]
    nhrid <- ncid$var[["summary_records"]]
    lehrid <- ncid$var[["last_element_highest"]]
    
    if (is.T(x)) {
     sdid <- ncid$var[["summary_sd"]]
     sdid.djf <- ncid$var[["summary_sd_DJF"]]
     sdid.mam <- ncid$var[["summary_sd_MAM"]]
     sdid.jja <- ncid$var[["summary_sd_JJA"]]
     sdid.son <- ncid$var[["summary_sd_SON"]]
     nlrid <- ncid$var[["summary_lows"]]
     lelrid <- ncid$var[["last_element_lowest"]]
    } else if (is.precip(x)) {
      muid <- ncid$var[["summary_wetmean"]]
      muid.djf <- ncid$var[["summary_wetmean_DJF"]]
      muid.mam <- ncid$var[["summary_wetmean_MAM"]]
      muid.jja <- ncid$var[["summary_wetmean_JJA"]]
      muid.son <- ncid$var[["summary_wetmean_SON"]]
      fwid <- ncid$var[["summary_wetfreq"]]
      fwid.djf <- ncid$var[["summary_wetfreq_DJF"]]
      fwid.mam <- ncid$var[["summary_wetfreq_MAM"]]
      fwid.jja <- ncid$var[["summary_wetfreq_JJA"]]
      fwid.son <- ncid$var[["summary_wetfreq_SON"]]
      tdmuid <- ncid$var[["summary_trend_wetmean"]]
      tdmuid.djf <- ncid$var[["summary_trend_wetmean_DJF"]]
      tdmuid.mam <- ncid$var[["summary_trend_wetmean_MAM"]]
      tdmuid.jja <- ncid$var[["summary_trend_wetmean_JJA"]]
      tdmuid.son <- ncid$var[["summary_trend_wetmean_SON"]]
      tdfwid <- ncid$var[["summary_trend_wetfreq"]]
      tdfwid.djf <- ncid$var[["summary_trend_wetfreq_DJF"]]
      tdfwid.mam <- ncid$var[["summary_trend_wetfreq_MAM"]]
      tdfwid.jja <- ncid$var[["summary_trend_wetfreq_JJA"]]
      tdfwid.son <- ncid$var[["summary_trend_wetfreq_SON"]]
      lrid <- ncid$var[["summary_lastrains"]]
      sigma2id <- ncid$var[["summary_sigma2"]]
      sigma2id.djf <- ncid$var[["summary_sigma2_DJF"]]
      sigma2id.mam <- ncid$var[["summary_sigma2_MAM"]]
      sigma2id.jja <- ncid$var[["summary_sigma2_JJA"]]
      sigma2id.son <- ncid$var[["summary_sigma2_SON"]]
      tsigma2id <- ncid$var[["summary_trend_sigma2"]]
      tsigma2id.djf <- ncid$var[["summary_trend_sigma2_DJF"]]
      tsigma2id.mam <- ncid$var[["summary_trend_sigma2_MAM"]]
      tsigma2id.jja <- ncid$var[["summary_trend_sigma2_JJA"]]
      tsigma2id.son <- ncid$var[["summary_trend_sigma2_SON"]]
      mwslid <- ncid$var[["summary_mean_wetdur"]]
      mdslid <- ncid$var[["summary_mean_drydur"]]
    } else {
      sdid <- ncid$var[["summary_sd"]]
      sdid.djf <- ncid$var[["summary_sd_DJF"]]
      sdid.mam <- ncid$var[["summary_sd_MAM"]]
      sdid.jja <- ncid$var[["summary_sd_JJA"]]
      sdid.son <- ncid$var[["summary_sd_SON"]]
      lelrid <- ncid$var[["last_element_lowest"]]
    }
    
    ## Appending the data after those that already exist:
    if (verbose) print(paste('Adjust start[2] so tht data is added after',ns,'stations'))
    start[2] <- start[2] + ns
    
  } else {
    if (verbose) print(paste('Creating file',file))
    if (is.T(x)) ncid <- nc_create(file,vars=list(ncvar,lonid,latid,altid,locid,stid,cntrid, 
                                                   fyrid,lyrid,nvid,meanid,meanid.djf,meanid.mam,meanid.jja,meanid.son,
                                                   sdid,sdid.djf,sdid.mam,sdid.jja,sdid.son,maxid,minid,nhrid,nlrid,
                                                   tdid,tdid.djf,tdid.mam,tdid.jja,tdid.son,lehrid,lelrid)) else
        if (is.precip(x)) ncid <- nc_create(file,vars=list(ncvar,lonid,latid,altid,locid,stid,cntrid, 
                                            fyrid,lyrid,nvid,meanid,meanid.djf,meanid.mam,meanid.jja,meanid.son,
                                            maxid,minid,nhrid,muid,muid.djf,muid.mam,muid.jja,muid.son,
                                            fwid,fwid.djf,fwid.mam,fwid.jja,fwid.son,
                                            tdid,tdid.djf,tdid.mam,tdid.jja,tdid.son,
                                            tdmuid,tdmuid.djf,tdmuid.mam,tdmuid.jja,tdmuid.son,
                                            tdfwid,tdfwid.djf,tdfwid.mam,tdfwid.jja,tdfwid.son,lrid,lehrid,
                                            sigma2id,sigma2id.djf,sigma2id.mam,sigma2id.jja,sigma2id.son,
                                            tsigma2id,tsigma2id.djf,tsigma2id.mam,tsigma2id.jja,tsigma2id.son,
                                            mwslid,mdslid)) else 
        ncid <- nc_create(file,vars=list(ncvar,lonid,latid,altid,locid,stid,cntrid, 
                                          fyrid,lyrid,nvid,meanid,meanid.djf,meanid.mam,meanid.jja,meanid.son,
                                          sdid,sdid.djf,sdid.mam,sdid.jja,sdid.son,tdid,
                                          tdid.djf,tdid.mam,tdid.jja,tdid.son,maxid,minid,nhrid,lehrid,lelrid))
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
  
  if (verbose) print('Add summary statistics: mean')
  ave[insufficient] <- missval; ave.djf[insufficient] <- missval
  ave.mam[insufficient] <- missval; ave.jja[insufficient] <- missval
  ave.son[insufficient] <- missval; td[insufficient] <- missval
  td.djf[insufficient] <- missval; td.mam[insufficient] <- missval
  td.jja[insufficient] <- missval; td.son[insufficient] <- missval
  mx[insufficient] <- missval; mn[insufficient] <- missval; 
  nhr[insufficient] <- missval; lehr[insufficient] <- missval
  
  ncvar_put( ncid, meanid, ave,start=start[2],count=count[2])
  ncvar_put( ncid, meanid.djf, ave.djf,start=start[2],count=count[2])
  ncvar_put( ncid, meanid.mam, ave.mam,start=start[2],count=count[2])
  ncvar_put( ncid, meanid.jja, ave.jja,start=start[2],count=count[2])
  ncvar_put( ncid, meanid.son, ave.son,start=start[2],count=count[2])
  if (verbose) print('Add summary statistics: trend')
  ncvar_put( ncid, tdid, td  ,start=start[2],count=count[2])
  ncvar_put( ncid, tdid.djf, td.djf,start=start[2],count=count[2])
  ncvar_put( ncid, tdid.mam, td.mam,start=start[2],count=count[2])
  ncvar_put( ncid, tdid.jja, td.jja,start=start[2],count=count[2])
  ncvar_put( ncid, tdid.son, td.son,start=start[2],count=count[2])
  if (verbose) print('Add summary statistics: max, min')
  ncvar_put( ncid, maxid, mx, start=start[2],count=count[2])
  ncvar_put( ncid, minid, mn, start=start[2],count=count[2])
  if (verbose) print('Add summary statistics: records')
  ncvar_put( ncid, nhrid, nhr, start=start[2],count=count[2])
  ncvar_put( ncid, lehrid, lehr, start=start[2],count=count[2])
  if (is.T(x)) {
    if (verbose) print('extra for temperature')
    std[insufficient] <- missval; std.djf[insufficient] <- missval
    std.mam[insufficient] <- missval; std.jja[insufficient] <- missval
    std.son[insufficient] <- missval; nlr[insufficient] <- missval
    lelr[insufficient] <- missval; 
    ncvar_put( ncid, sdid, std, start=start[2],count=count[2])
    ncvar_put( ncid, sdid.djf, std.djf, start=start[2],count=count[2])
    ncvar_put( ncid, sdid.mam, std.mam, start=start[2],count=count[2])
    ncvar_put( ncid, sdid.jja, std.jja, start=start[2],count=count[2])
    ncvar_put( ncid, sdid.son, std.son, start=start[2],count=count[2])
    ncvar_put( ncid, nlrid, nlr, start=start[2],count=count[2])
    ncvar_put( ncid, lelrid, lelr, start=start[2],count=count[2])
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
    lr[insufficient] <- missval
    ncvar_put( ncid, muid, mu,start=start[2],count=count[2])
    ncvar_put( ncid, muid.djf, mu.djf,start=start[2],count=count[2])
    ncvar_put( ncid, muid.mam, mu.mam,start=start[2],count=count[2])
    ncvar_put( ncid, muid.jja, mu.jja,start=start[2],count=count[2])
    ncvar_put( ncid, muid.son, mu.son,start=start[2],count=count[2])
    ncvar_put( ncid, fwid, fw,start=start[2],count=count[2])
    ncvar_put( ncid, fwid.djf, fw.djf,start=start[2],count=count[2])
    ncvar_put( ncid, fwid.mam, fw.mam,start=start[2],count=count[2])
    ncvar_put( ncid, fwid.jja, fw.jja,start=start[2],count=count[2])
    ncvar_put( ncid, fwid.son, fw.son,start=start[2],count=count[2])
    ncvar_put( ncid, tdfwid, tdfw,start=start[2],count=count[2])
    ncvar_put( ncid, tdfwid.djf, tdfw.djf,start=start[2],count=count[2])
    ncvar_put( ncid, tdfwid.mam, tdfw.mam,start=start[2],count=count[2])
    ncvar_put( ncid, tdfwid.jja, tdfw.jja,start=start[2],count=count[2])
    ncvar_put( ncid, tdfwid.son, tdfw.son,start=start[2],count=count[2])
    ncvar_put( ncid, tdmuid, tdmu,start=start[2],count=count[2])
    ncvar_put( ncid, tdmuid.djf, tdmu.djf,start=start[2],count=count[2])
    ncvar_put( ncid, tdmuid.mam, tdmu.mam,start=start[2],count=count[2])
    ncvar_put( ncid, tdmuid.jja, tdmu.jja,start=start[2],count=count[2])
    ncvar_put( ncid, tdmuid.son, tdmu.son,start=start[2],count=count[2])
    ncvar_put( ncid, lrid, lr, start=start[2],count=count[2])
    ncvar_put( ncid, sigma2id, sigma2,start=start[2],count=count[2])
    ncvar_put( ncid, sigma2id.djf, sigma2.djf,start=start[2],count=count[2])
    ncvar_put( ncid, sigma2id.mam, sigma2.mam,start=start[2],count=count[2])
    ncvar_put( ncid, sigma2id.jja, sigma2.jja,start=start[2],count=count[2])
    ncvar_put( ncid, sigma2id.son, sigma2.son,start=start[2],count=count[2])
    ncvar_put( ncid, tsigma2id, tsigma2,start=start[2],count=count[2])
    ncvar_put( ncid, tsigma2id.djf, tsigma2.djf,start=start[2],count=count[2])
    ncvar_put( ncid, tsigma2id.mam, tsigma2.mam,start=start[2],count=count[2])
    ncvar_put( ncid, tsigma2id.jja, tsigma2.jja,start=start[2],count=count[2])
    ncvar_put( ncid, tsigma2id.son, tsigma2.son,start=start[2],count=count[2])
    if (verbose) print('Mean spell length')
    ncvar_put( ncid, mwslid, mwsl,start=start[2],count=count[2])
    ncvar_put( ncid, mdslid, mdsl,start=start[2],count=count[2])
  } else {
    if (verbose) print(paste('extra for',varid(x)[1]))
    std[insufficient] <- missval; std.djf[insufficient] <- missval
    std.mam[insufficient] <- missval; std.jja[insufficient] <- missval
    std.son[insufficient] <- missval; nlr[insufficient] <- missval
    lelr[insufficient] <- missval; 
    ncvar_put( ncid, sdid, std, start=start[2],count=count[2])
    ncvar_put( ncid, sdid.djf, std.djf, start=start[2],count=count[2])
    ncvar_put( ncid, sdid.mam, std.mam, start=start[2],count=count[2])
    ncvar_put( ncid, sdid.jja, std.jja, start=start[2],count=count[2])
    ncvar_put( ncid, sdid.son, std.son, start=start[2],count=count[2])
    ncvar_put( ncid, lelrid, lelr, start=start[2],count=count[2])
  }
  
  ## There are some times problems saving text data, and there seems to be some 
  ## inconsistency in the ncdf4 package. To by-pass this problem, we had to make
  ## the following code more complicated. There seems to be a mix-up between the 
  ## dimensions sometimes.
  if (verbose) print('Saving textual information')
  test <- try(ncvar_put( ncid, locid, loc(y),start=c(1,start[2]),count=c(namelength,count[2])))
  if (inherits(test,'try-error'))
    try(ncvar_put( ncid, locid, loc(y),start=c(start[2],1),count=c(count[2],namelength)))
  test <- try(ncvar_put( ncid, stid, as.character(stid(y)),c(1,start[2]),count=c(namelength,count[2])))
  if (inherits(test,'try-error'))
    try(ncvar_put( ncid, stid, as.character(stid(y)),c(start[2],1),count=c(count[2],namelength)))
  test <- try(ncvar_put( ncid, cntrid, cntr(y),start=c(1,start[2]),count=c(namelength,count[2])))
  if (inherits(test,'try-error'))
    try(ncvar_put( ncid, cntrid, cntr(y),start=c(start[2],1),count=c(count[2],namelength)))
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
write2ncdf4.pca <- function(x,file='esd.pca.nc',prec='short',verbose=FALSE,scale=0.01,offset=0,missval=-99) {
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
}

write2ncdf4.eof <- function(x,file='eof.nc',prec='short',scale=10,offset=NULL,torg="1970-01-01",missval=-999) {
}

  
write2ncdf4.dsensemble <- function(x,file='esd.dsensemble.nc',prec='short',offset=0,scale=0.1,
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

  ncnew <- nc_create( file, varlist,verbose=verbose)
  
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
  if (verbose) print(paste('Finished successfully - file', file))  
}
