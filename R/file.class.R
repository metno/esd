# file.class.R

## Used to check the contents in netCDF file - to use in retrieve to call retrieve.dsenemble,
## retrieve.eof or retrieve.station rather than the standard form to read field objects.
## Assumes that empty class attribute means a field object
file.class <- function(ncfile) {
  nc <- nc_open(ncfile)
  dimnames <- names(nc$dim)
  class.x <- ncatt_get(nc,0,'class')
  nc_close(nc)
  class.x$dimnames <- dimnames
  return(class.x)
}
