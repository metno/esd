#' Reads and extracts the attribute information in a netCDF files
#' 
#' @param filename Name of netCDF file including path
#' 
#' @export getatt
getatt <- function(filename) {
  ncid <- ncdf4::nc_open(filename)
  ncdf4::nc_close(ncid)
  return(ncid)
}