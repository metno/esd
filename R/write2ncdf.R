save2ncdf <- function(x,file,...) UseMethod("save2ncdf")

## https://www.unidata.ucar.edu/software/netcdf/docs/netcdf/CDF-Data-Types.html:
## short: 16-bit signed integers. The short type holds values between -32768 and 32767.

save2ncdf.station <- function(x,file,prec='short',missval=-99.9,offs=0,
                              scalf=0.1,torg='1899-12-31',verbose=FALSE) {

    require(netcdf4)
    if(verbose) print("save2ncdf.station")

    ## Save a station object as a netCDF file using the short-type combined with add_offsett and scale_factor to reduce the size.
    ## Examine the station object: dimensions and attributes
    d <- dim(x)

    if (is.null(d)) d <- c(length(x),1)
    if (verbose) print(paste("dimension of station",paste(d,collapse=" x ")))
    atts <- names(attributes(x))
    if (verbose) print(atts)
    iattr2ncdf <- !is.element(atts,c('dim','dimnames','index',
                                     'station_id','variable','unit',
                                     'longname','history','location','country'))
    attr <- attr[iattr2ncdf]
    na <- length(atts); la <- rep(0,na)
    attrprec <- rep("character",na)
    for (i in 1:na) {
        la[i] <- length(attr(x,atts[i]))
        attrprec[i] <- switch(tolower(atts[i]),
                              'longitude'='float',latitude='float',altitude='float',
                              'wmo_id'='integer','class'='character','location'='character',
                              'country'='character','aspect'='character','source'='character',
                              'quality'='character','url'='character','reference'='character',
                              'info'='character')
    }
    attr(x,'quality') <- as.character(attr(x,'quality'))
    if (verbose) print(paste('attributes:', paste(atts,collapse=', '),
                             '; types:',paste(attrprec,collapse=', ')))
    y <- t(coredata(x))
    attributes(y) <- NULL
    y[!is.finite(y)] <- missval
}
