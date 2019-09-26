# Export esd-objects to netCDF files
# 
# \code{esd2ncdf4} saves the data and its main structure (metadata) as netCDF
# files according to the CF-convention \url{http://cfconventions.org/}.
# 
# 
# @aliases esd2ncdf4 esd2ncdf4.field esd2ncdf4.station
# esd2ncdf4.eof esd2ncdf4.pca esd2ncdf4.dsensemble
# 
# @param x esd object
# @param fname file name
# @param scale scaling factor
# @param offset offset
# @param torg time origin
# @param \dots additional arguments
# 
esd2ncdf4 <- function(x,...) UseMethod("esd2ncdf4")

esd2ncdf4.default <- function(x,...) {
}

esd2ncdf4.field <- function(x,fname='field.nc',scale=10,offset=NULL,torg="1970-01-01",verbose=FALSE) {
  if(verbose) print('esd2ncdf4.field')
  dimlon <- ncdim_def( "longitude", "degree_east", lon(x) )
  dimlat <- ncdim_def( "latitude", "degree_north", lat(x) )
  if (inherits(index(x),c('numeric','integer')))
      index(x) <- as.Date(paste(index(x),'-01-01',sep=''))
  
  dimtim <- ncdim_def( "time", paste("days since",torg),
                      as.numeric(as.Date(index(x),origin=torg)) )
  x4nc <- ncvar_def(varid(x)[1], unit(x)[1], list(dimlon,dimlat,dimtim), -1, 
                    longname=attr(x,'longname'), prec="short")
     
     # Create a netCDF file with this variable
  ncnew <- nc_create( fname, x4nc )

  y <- coredata(x); attributes(y) <- NULL
  if (is.null(offset)) offset <- mean(y,na.rm=TRUE)
  y <- round(scale*(y-offset))
  # Write some values to this variable on disk.
  ncvar_put( ncnew, x4nc, round(y) )
  ncvar_put( ncnew, x4nc, round(y) )
  ncatt_put( ncnew, x4nc, "add_offset", mean(y,na.rm=TRUE), prec="float" )
  ncatt_put( ncnew, x4nc, "scale_factor", 0.01, prec="float" ) 
  ncatt_put( ncnew, x4nc, "_FillValue", -99, prec="float" ) 
  ncatt_put( ncnew, x4nc, "missing_value", -99, prec="float" ) 
  ncatt_put( ncnew, 0, "description", 
             "Saved from esd using esd2ncdf4")
  nc_close(ncnew)
}

esd2ncdf4.station <- function(x,fname='station.nc',scale=10,offset=NULL,torg="1970-01-01") {
}


esd2ncdf4.eof <- function(x,fname='eof.nc',scale=10,offset=NULL,torg="1970-01-01") {
}

esd2ncdf4.pca <- function(x,fname='pca.nc',scale=10,offset=NULL,torg="1970-01-01") {
}

esd2ncdf4.dsensemble <- function(x,fname='dsensemble.nc',scale=10,offset=NULL,torg="1970-01-01") {
}
