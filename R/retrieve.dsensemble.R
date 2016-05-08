retrieve.dsensemble <- function(ncfile,param="auto",type="ncdf4",
                             path=NULL,verbose=FALSE,...) {
   if ((type=="ncdf") | (class(ncfile)=="ncdf")) { ##(library("ncdf",logical.return=TRUE)) {
        nc <- open.ncdf(file.path(path,ncfile))
        dimnames <- names(nc$dim)
        lon <- get.var.ncdf(nc,dimnames[grep("lon|x",tolower(dimnames))])
        lat <- get.var.ncdf(nc,dimnames[grep("lat|y",tolower(dimnames))])
        close.ncdf(nc)
        if ( (length(dim(lon))==1) & (length(dim(lat))==1) ) {
            if (verbose) print('Regular grid field found')
            X <- retrieve.ncdf(ncfile,path=path,param=param,verbose=verbose,...)
        } else {
            if (verbose) print('Irregular grid field found')
            X <- retrieve.rcm(ncfile,path=path,param=param,verbose=verbose,...) 
        }
    } else if ((type=="ncdf4") | (class(ncfile)=="ncdf4")) {##(library("ncdf4",logical.return=TRUE)) {
        nc <- nc_open(file.path(path,ncfile))
        dimnames <- names(nc$dim)
        lon <- ncvar_get(nc,dimnames[grep("lon|x",tolower(dimnames))])
        lat <- ncvar_get(nc,dimnames[grep("lat|y",tolower(dimnames))])
        nc_close(nc)
        if ( (length(dim(lon))==1) & (length(dim(lat))==1) )  {
            if (verbose) print('Regular grid field found')
            X <- retrieve.ncdf4(ncfile,path=path,param=param,verbose=verbose,...)
        }
        else {
            if (verbose) print('Irregular grid field found')
            X <- retrieve.rcm(ncfile,path=path,param=param,verbose=verbose,...) 
        }
    } else {
      print("No suitable ncdf or ncdf4 libraries found to read your file or data")
    }
}

