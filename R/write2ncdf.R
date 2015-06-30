save2ncdf <- function(x,file,...)
    UseMethod("save2ncdf")

## https://www.unidata.ucar.edu/software/netcdf/docs/netcdf/CDL-Data-Types.html:
## short: 16-bit signed integers. The short type holds values between -32768 and 32767. 

save2ncdf.station <- function(x,file,prec='short',missval=-99.9,offs=0,
                              scalf=0.1,torg='1899-12-31',verbose=FALSE) {
    ##require(ncdf4)
    if (verbose) print('save2ncdf.station')
    
    ## Save a station object as a netCDF file using the short-type combined with add_offset and scale_factor
    ## to reduce the size.   
    ## Examine the station object: dimensions and attributes  
    d <- dim(x)

    if (is.null(d)) d <- c(length(x),1)
    if (verbose) print(paste('dimension of station: ',paste(d, collapse= ' x ')))
    atts <- names(attributes(x))
    if (verbose) print(atts)
    iattr2ncdf <- !is.element(atts,c('dim','dimnames','index',
                                     'station_id','variable','unit',
                                     'longname','history','location','country'))
    atts <- atts[iattr2ncdf]
    na <- length(atts); la <- rep(0,na)
    attrprec <- rep('character',na)
    for (i in 1:na) {
        ##if (verbose) print(atts[i])
        la[i] <- length(attr(x,atts[i]))
        attrprec[i] <- switch(tolower(atts[i]),
                              'longitude'='float',latitude='float',altitude='float',
                              'wmo_id'='integer','class'='character','location'='character',
                              'country'='character','aspect'='character','source'='character',
                              'quality'='character','url'='character','reference'='character',
                              'info'='character')
    }
    attr(x,'quality') <- as.character(attr(x,'quality'))
    if (verbose) print(paste('attributes:', paste(atts, collapse=', '),
                             '; types:',paste(attrprec, collapse=', ')))
    
    y <- t(coredata(x))
    attributes(y) <- NULL
    y[!is.finite(y)] <- missval
    
    y <- round((y - offs)/scalf)
    dim(y) <- c(d[2],d[1])
    
    if (is.null(attr(x,'calendar'))) attr(x,'calendar') <- 'standard'
    if (!is.na(attr(x,'calendar'))) calendar <- attr(x,'calendar') else
    calendar <- 'standard'
    if (class(index(x))=='Date')
        time <- julian(index(x)) - julian(as.Date(torg))
    else if (inherits(x,'annual'))
        time <- julian(as.Date(paste(year(x),'01-01',sep='-')))-julian(as.Date(torg))
    
    ## Attributes with same number of elements as stations are saved as variables
    
    ## Define the dimensions
    if (verbose) print('define dimensions')
    ##browser()
    if (verbose) print(stid(x))
    dimSt <- ncdim_def( "stid", "number", as.numeric(stid(x)),longname='station_number')
    dimTi <- ncdim_def( "time", paste("days since",torg), time, calendar=calendar )
    
    if (verbose) print('define variable')
    ncvar <- ncvar_def(varid(x)[1],unit(x)[1], list(dimSt,dimTi),missval=missval, attr(x,'longname')[1], prec=prec)
    
    if (verbose) print(paste('create netCDF-file',file))
    ncid <- nc_create(file,ncvar)
    if (verbose) {str(y); print(summary(c(y)))}
    ncvar_put( ncid, ncvar, y)
    ncatt_put( ncid, ncvar, 'add_offset',offs,prec='float')
    ncatt_put( ncid, ncvar, 'scale_factor',scalf,prec='float')
    ncatt_put( ncid, ncvar, 'missing_value',missval,prec='float')
    ncatt_put( ncid, ncvar, 'location',paste(loc(x),collapse=", "),prec='character')
    ncatt_put( ncid, ncvar, 'country',paste(cntr(x),collapse=", "),prec='character')
    
    if (verbose) print('save attributes')
    for (i in 1:na) {
        if ( (la[i]==d[2]) & !is.na(attr(x,atts[i])[1]) ) {
            if (verbose) print(paste(atts[i],attrprec[i],sep=': '))
            ncatt_put( ncid, ncvar, atts[i], attr(x,atts[i]), prec=attrprec[i] )
        }
    }
                                        # global attributes:
    if (verbose) print('save global attributes')
    for (i in 1:na) {
        if (la[i]!=d[2]) {
            if (verbose) print(paste(atts[i],attrprec[i],sep=': '))
            ncatt_put( ncid, 0, atts[i], attr(x,atts[i]), prec=attrprec[i] )
        }
    }
    
    nc_close(ncid)
    if (verbose) print('close')
}
