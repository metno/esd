#' @exportS3Method
#' @export retrieve.dsensemble
retrieve.dsensemble <- function(ncfile,...) {
  ## REB 2026-01-22 
  args <- list(...)
  if (!is.null(args$param)) param <- args$param else param <-"auto"
  if (!is.null(args$verbose)) verbose <- args$verbose else verbose <-FALSE
  path <- args$path
    
  if(verbose) print("retrieve.dsensemble")
  nc <- nc_open(file.path(path,ncfile))
  #browser()
  dimnames <- names(nc$dim)
  lon <- ncvar_get(nc,dimnames[grep("lon|x",tolower(dimnames))])
  lat <- ncvar_get(nc,dimnames[grep("lat|y",tolower(dimnames))])
  nc_close(nc)
  if ( (length(dim(lon))==1) & (length(dim(lat))==1) )  {
    if (verbose) print('Regular grid field found')
    X <- retrieve.ncdf4(ncfile,path=path,param=param,verbose=verbose,...)
  } else {
    if (verbose) print('Irregular grid field found')
    X <- retrieve.rcm(ncfile,path=path,param=param,verbose=verbose,...) 
  }
}

