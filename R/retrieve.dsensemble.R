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

## These small functions are common code that simplify saving data as netCDF 
prepare.pca4nc <- function(x,verbose=FALSE,scale=0.01,offset=0,missing=-99,prec='short') {
  if (verbose) print('prepare.pca4nc')
  pcaatts <- names(attributes(x))
  pattern <- attr(x,'pattern')
  if (verbose) print(pcaatts)
  ## set up the dimensions of the PCA
  dimpca <- ncdim_def( "pca", "index", 1:dim(x)[2] )
  dimtim <- ncdim_def( "time_pca", "year", index(x) )
  dimxy <- ncdim_def( "space_pca", "index", 1:dim(pattern)[1] )
  dimsea <- ncdim_def( "season", "index", 1:4 )
  pca <- ncvar_def("pca", "weights", list(dimsea,dimtim,dimpca), missing, 
                    longname='principal components', prec=prec)
  pat <- ncvar_def("pattern_pca", "weights", list(dimsea,dimxy,dimpca), missing, 
                    longname='principal component analysis patterns', prec=prec)
  lon <- ncvar_def("longitude", "degree_east", dimxy,miss, 
                    longname='longitude', prec='float')
  lat <- ncvar_def("latitude", "degree_north", dimxy,miss, 
                    longname='latitude', prec='float')
  alt <- ncvar_def("altitude", "m", dimxy,miss, 
                    longname='altitude', prec='float')
  stid <- ncvar_def("station_id", "number", dimxy,miss, 
                    longname='station ID', prec="char")
  loc <- ncvar_def("location", "name", dimxy,miss, 
                    longname='location name', prec="char")
  cntr <- ncvar_def("country", "name", dimxy,miss, 
                    longname='country name', prec="char")
  src <- ncvar_def("src", "name", dimxy,miss, 
                    longname='source', prec="char")
  lambda <- ncvar_def("lambda", "number", dimpca,miss, 
                    longname='eigenvalues', prec="float")
  ncstr <- list(pca,pat,tim,lon,lat,alt,stid,lambda,loc,cntr,src)
  return(ncstr)
  }

add.pca2nc <- function(x,nc,ncstr,verbose=FALSE,scale=0.01,offset=0,missing=-99,prec='short') {
  if (verbose) print('add.pca2nc')
  pcaatts <- names(attributes(x))
  pattern <- attr(x,'pattern')
  pattern[!is.finite(pattern)] <- missing
  dpat <- dim(pattern); attributes(pattern) <- NULL; dpat -> dim(pattern)
  pca <- coredata(x)
  if (is.numeric(index(x))) index(x) <- year(x)
  dpca <- dim(pca); attributes(pca) <- NULL; dpca -> dim(pca)
  ncvar_put( nc, ncstr$pca, round((pca - offset)/scale) )
  ncatt_put( nc, ncstr$pca, "add_offset", offset, prec="float" )
  ncatt_put( nc, ncstr$pca, "scale_factor", scale, prec="float" ) 
  ncatt_put( nc, ncstr$pca, "_FillValue", missval, prec="float" ) 
  ncatt_put( nc, ncstr$pca, "missing_value", missval, prec="float" ) 
  ncatt_put( nc, ncstr$pca, "dimensions_pca", attr(pattern,'dimensions'), prec="char" ) 
  ncvar_put( nc, ncstr$pat, round((pattern - offset)/scale) )
  ncatt_put( nc, ncstr$pat, "add_offset", offset, prec="float" )
  ncatt_put( nc, ncstr$pat, "scale_factor", scale, prec="float" ) 
  ncatt_put( nc, ncstr$pat, "_FillValue", missval, prec="float" ) 
  ncatt_put( nc, ncstr$pat, "missing_value", missval, prec="float" ) 
  ncvar_put( nc, ncstr$tim, index(x) )
  ncvar_put( nc, ncstr$lon, lon(x) )
  ncvar_put( nc, ncstr$lat, lat(x) )
  ncvar_put( nc, ncstr$alt, alt(x) )
  ncvar_put( nc, ncstr$loc, loc(x) )
  ncvar_put( nc, ncstr$stid, stid(x) )
  ncvar_put( nc, ncstr$cntr, cntr(x) )
  ncvar_put( nc, ncstr$src, src(x) )
  ncvar_put( nc, ncstr$lambde, attr(x,'eigenvalues') )
  ncatt_put( nc, ncstr$pca, "history", paste(attr(x,'history'),collapse=';'), prec="char" )
}

prepare.eof4nc <- function(x,nc) {
}


  
write2ncdf4.dsensemble <- function(x,fname='field.nc',prec='short',scale=0.1,offset=NULL,
                              torg="1970-01-01",missval=-999,verbose=FALSE) {
  if (verbose) print('write2ncdf4.field')
  names.x <- names(x)
  class.x <- class(x)
  ngcms <- length(x) - 2
  ## Get the two first elements of the list which contain information common to the rest
  info <- x$info
  if (!is.null(x$pca)) pca <- x$pca else pca <- NULL
  if (!is.null(x$eof)) eof <- x$eof else eof <- NULL
  
  ## Clear these so that the rest are zoo objects describing the model runs
  x$info <- NULL
  x$pca <- NULL
  names.x <- names.x[-c(1:2)]
  if (verbose) print(names.x)
  
  ## Define the variables and dimensions associated with the PCA/EOFs
  if (!is.null(pca)) varlist <- prepare.pca4nc(pca)
  if (!is.null(eof)) varlist <- prepare.eof4nc(eof)
  if (is.numeric(index(x))) index(x) <- year(x)
  dimtim <- ncdim_def( "time", paste("days since",as.Date(0)), index(x) )
  dimens <- ncdim_def( "ensemble_member", 'index', 1:length(x) )
  varlist$gcm <- ncvar_def("pca", "weights", list(dimens,dimtim,varlist$dimpca), missing, 
                    longname='principal components', prec=prec)
  varlist$names <- ncvar_def("names", "text", list(dimens,dimtim,varlist$dimpca), missing, 
                    longname='names of GCM and run', prec="char")
     
     # Create a netCDF file with this variable
  ncnew <- nc_create( fname, varlist )
  add.pca2nc(x,ncnew,varlist)
  X <- unlist(x); dim(X) <- c(length(x),dim(x[[1]]))
  ncvar_put( ncnew, varlist$gcm, round((X - offset)/scale) )
  # Write some values to this variable on disk.
  ncatt_put( ncnew, varlist$gcm, "add_offset", offset, prec="float" )
  ncatt_put( ncnew, varlist$gcm, "scale_factor", scale, prec="float" ) 
  ncatt_put( ncnew, varlist$gcm, "_FillValue", missval, prec="float" ) 
  ncatt_put( ncnew, varlist$gcm, "missing_value", missval, prec="float" ) 
  ncvar_put( ncnew, varlist$names, names.x )
  ncatt_put( ncnew, 0, "description", "Saved from esd using write2ncdf4.dsensemble")
  ncatt_put( ncnew, 0, "class", paste(class.x,collapse='-'))
  nc_close(ncnew)
}

## Used to check the contents in netCDF file - to use in retrieve to call retrieve.dsenemble,
## retrieve.eof or retrieve.station rather than the standard form to read field objects.
## Assumes that empty class attribute means a field object
file.class <- function(ncfile,path=NULL,default='field',type="ncdf4") {
  if (type=='ncdf4') {
    nc <- nc_open(file.path(path,ncfile))
    dimnames <- names(nc$dim)
    class <- ncatt_get(nc,0,'class')
    close.ncdf(nc)
  } else {
    nc <- open.ncdf(file.path(path,ncfile))
    dimnames <- names(nc$dim)
    lon <- get.att.ncdf(nc,0,'class')
    close.ncdf(nc)
  }
  attr(class,'dimnames') <- dimnames
  return(class)
}
