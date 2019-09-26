#' Used to check the contents in netCDF file.
#'
#' Assumes that empty class attribute means a field object
#'
#' @param ncfile filename of netcdf file
#'
#' @seealso check.ncdf4 retrieve
#'
#' @export
file.class <- function(ncfile) {
  nc <- nc_open(ncfile)
  dimnames <- names(nc$dim)
  class.x <- ncatt_get(nc,0,'class')
  nc_close(nc)
  class.x$dimnames <- dimnames
  return(class.x)
}
