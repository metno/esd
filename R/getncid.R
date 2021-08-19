#' Generic function to retrieve information from a netcdf file
#'
#' @param filename Names of files
#' @param path Path to files, if not provided in filename
#' @param verbose if TRUE print information about progress
#'
#' @export
getncid <- function(filename,path=NULL,verbose=FALSE) {
  if(verbose) print("getncid")
  if(!is.null(path)) filename <- file.path(path,filename)
  if(verbose) print("Collect information stored in the netCDF header")
  cid <- getatt(filename)
  ncid <- ncdf4::nc_open(filename)
  ncid2 <- check.ncdf4(ncid,param="auto")
  ncdf4::nc_close(ncid)
  cid$dates <- paste(range(ncid2$time$vdate),collapse=",")
  cid$model <- ncid2$model
  cid$project_id <- cid$model$project_id
  if(is.null(cid$project_id) & !is.null(cid$model$title)) {
    if(grepl("CMIP5",cid$model$title)) cid$project_id <- "CMIP5"
  }
  return(cid)
}