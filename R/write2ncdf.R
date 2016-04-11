
## https://www.unidata.ucar.edu/software/netcdf/docs/netcdf/CDF-Data-Types.html:
## short: 16-bit signed integers. The short type holds values between -32768 and 32767.

write2ncdf4 <- function(x,...) UseMethod("write2ncdf4")

write2ncdf4.default <- function(x,...) {
}

write2ncdf4.field <- function(x,fname='field.nc',prec='short',scale=0.1,offset=NULL,
                              torg="1970-01-01",missval=-999,verbose=FALSE) {
  if (verbose) print('write2ncdf4.field')

  y <- coredata(x)
  if (is.null(offset)) offset <- mean(y,na.rm=TRUE)
  if (is.null(scale)) scale <- 1
  y <- t(y)
  y[!is.finite(y)] <- missval
  y <- round((y-offset)/scale)
  if (verbose) print(attr(y,'dimensions'))
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
  ncvar_put( ncnew, x4nc, round(y) )
  ncatt_put( ncnew, x4nc, "add_offset", offset, prec="float" )
  ncatt_put( ncnew, x4nc, "scale_factor", scale, prec="float" ) 
  ncatt_put( ncnew, x4nc, "_FillValue", missval, prec="float" ) 
  ncatt_put( ncnew, x4nc, "missing_value", missval, prec="float" ) 
  ncatt_put( ncnew, 0, "description", 
             "Saved from esd using write2ncdf4")
  nc_close(ncnew)
}


write2ncdf4.eof <- function(x,fname='eof.nc',prec='short',scale=10,offset=NULL,torg="1970-01-01",missval=-999) {
}

write2ncdf4.pca <- function(x,fname='pca.nc',prec='short',scale=10,offset=NULL,torg="1970-01-01",missval=-999) {
}

write2ncdf4.dsensemble <- function(x,fname='dsensemble.nc',prec='short',scale=10,offset=NULL,torg="1970-01-01",missval=-999) {
}

# https://www.unidata.ucar.edu/software/netcdf/docs/netcdf/CDL-Data-Types.html:
# short: 16-bit signed integers. The short type holds values between -32768 and 32767. 

write2ncdf4.station <- function(x,fname,prec='short',offset=0, missval=-999,
                                scale=0.1,torg='1899-12-31',verbose=FALSE) {
  #require(ncdf4)

  if (!inherits(x,"station")) stop('x argument must be a station object') 
  
  if (verbose) print('write.ncdf')
  
  ## Write a station object as a netCDF file using the short-type combined with add_offsetet and scale_factor
  ## to reduce the size.   
  ## Examine the station object: dimensions and attributes  

  ## Get time 
  nt <- dim(x)[1]
  ns <- dim(x)[2]
  
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
  
  y <- t(coredata(x))
  y[!is.finite(y)] <- missval
  y <- round((y - offset)/scale)
  dim(y) <- c(ns,nt)
  
  if (is.null(attr(x,'calendar')))
      attr(x,'calendar') <- 'standard'

  if (!is.na(attr(x,'calendar')))
      calendar <- attr(x,'calendar')
  else calendar <- 'standard'

  if (class(index(x))=='Date')
      time <- julian(index(x)) - julian(as.Date(torg))
  else if (inherits(x,'annual'))
      time <- julian(as.Date(paste(year(x),'01-01',sep='-')))-julian(as.Date(torg))
  
# Attributes with same number of elements as stations are saved as variables
  
# Define the dimensions
  if (verbose) print('Define dimensions')
  if (verbose) print(stid(x))
  dimS <- ncdim_def( name="stid", units="number",vals=c(1:ns))
  dimT <- ncdim_def( name="time", units=paste("days since",torg), vals=time, calendar=calendar)
  
  if (verbose) print('Define variable')

  lon <- lon(x)
  lat <- lat(x)
  alt <- alt(x)
  
  if (verbose) print(paste('create netCDF-file',fname))
     
  if (verbose) {str(y); print(summary(c(y)))}
  latid <- ncvar_def(name="lat",dim=list(dimS), units="degrees_north", missval=missval,longname="latitude", prec=prec,verbose=verbose)

  lonid <- ncvar_def(name="lon",dim=list(dimS), units="degrees_east", missval=missval,longname="longitude", prec=prec,verbose=verbose)
   altid <- ncvar_def(name="alt",dim=list(dimS), units="meters", missval=missval,longname="altitude", prec=prec,verbose=verbose)

  locid <- ncvar_def(name="loc",dim=list(dimS),units="strings",prec="char",longname="location",verbose=verbose)
  
  ncvar <- ncvar_def(name=varid(x)[1],dim=list(dimT,dimS), units=ifelse(unit(x)[1]=="°C", "degC",unit(x)[1]),longname=attr(x,'longname')[1], prec="float",compression=9,verbose=verbose)

  ncid <- nc_create(fname,vars=list(ncvar,lonid,latid,altid,locid)) ## vars)
  ncvar_put( ncid, ncvar, y)
  ncatt_put( ncid, ncvar, 'add_offset',offset,prec='float')
  ncatt_put( ncid, ncvar, 'scale_factor',scale,prec='float')
  ncatt_put( ncid, ncvar, 'missing_value',missval,prec='float')
  ncatt_put( ncid, ncvar, 'location',paste(loc(x),collapse=", "),prec='char')
  ncatt_put( ncid, ncvar, 'country',paste(cntr(x),collapse=", "),prec='character')
  ncvar_put( ncid, lonid, attr(x,"longitude"))
  ncvar_put( ncid, latid, attr(x,"latitude"))
  ncvar_put( ncid, altid, attr(x,"altitude"))
  ncvar_put( ncid, locid, as.array(attr(x,"location")))
  
  ## global attributes
  ncatt_put( ncid, 0, 'title', paste(levels(factor(attr(x,"info"))),collapse="/"))
  ncatt_put( ncid, 0, 'source', paste(levels(factor(attr(x,"source"))),collapse="/"))
  ncatt_put( ncid, 0, 'history', paste(unlist(attr(tmax,"history")),collapse="/"))
  ncatt_put( ncid, 0, 'references', paste(levels(factor(attr(tmax,"reference"))),collapse="/"))
  
  nc_close(ncid)
  if (verbose) print('close')
}
