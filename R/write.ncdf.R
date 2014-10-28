write.ncdf <- function(x,file,...)
  UseMethod("write2ncdf")

# https://www.unidata.ucar.edu/software/netcdf/docs/netcdf/CDL-Data-Types.html:
# short: 16-bit signed integers. The short type holds values between -32768 and 32767. 

write.ncdf.station <- function(x,file,prec='short',offs=0, missval=NA,
                               scalf=0.1,torg='1899-12-31',verbose=FALSE) {
  #require(ncdf4)

  if (!inherits(x,"station")) stop('x argument must be a station object') 
  
  if (verbose) print('write.ncdf')
  
  ## Write a station object as a netCDF file using the short-type combined with add_offset and scale_factor
  ## to reduce the size.   
  ## Examine the station object: dimensions and attributes  

  ## Get time 
  nt <- dim(x)[1]
  ns <- dim(x)[2]
  
  ## if (is.null(d)) d <- c(length(x),1)
  if (verbose) print(paste('Number of stations: ',paste(ns)))
  
  atts <- names(attributes(x))
 
  if (verbose) print(atts)
  ## attr2attr <- is.element(atts,c('dim','dimnames','index',
  ##                                 'station_id','variable','unit',
  ##                                'longname','history','location','country'))
  attr2attr <- is.element(atts,c('station_id','variable','unit','longname','location','country'))
  ##atts <- atts[iattr2ncdf]
  attr2var <- is.element(atts,c('longitude','latitude','altitude'))
  na <- length(atts); la <- rep(0,na)
  attrprec <- rep('character',na)
  ## for (i in 1:na) {
  ##  #if (verbose) print(atts[i])
  ##  la[i] <- length(attr(x,atts[i]))
  ##  attrprec[i] <- switch(tolower(atts[i]),
  ##                     'longitude'='float','latitude'='float','altitude'='float',
  ##                     'wmo_id'='integer','class'='character','location'='character',
  ##                     'country'='character','aspect'='character','source'='character',
  ##                     'quality'='character','url'='character','reference'='character',
  ##                     'info'='character')
  ##}
  attr(x,'quality') <- as.character(attr(x,'quality'))
  if (verbose) print(paste('attributes:', paste(atts, collapse=', '),
                           '; types:',paste(attrprec, collapse=', ')))
  
  y <- t(coredata(x))
  attributes(y) <- NULL
  y[!is.finite(y)] <- missval
    
  y <- round((y - offs)/scalf)
  dim(y) <- c(ns,nt)
  
  if (is.null(attr(x,'calendar')))
      attr(x,'calendar') <- 'standard'

  if (!is.na(attr(x,'calendar')))
      calendar <- attr(x,'calendar')
  else calendar <- 'standard'

  if (class(index(x))=='Date')
      time <- julian(index(x)) - julian(as.Date(torg))
  else if (inherits(x,'annual'))
      time <- julian(as.Date(paste(year(x),'01-01',sep='-')))-julian(as.Date(torg))
  
# Attributes with same number of elements as stations are saved as variables
  
# Define the dimensions
  if (verbose) print('Define dimensions')
  ##browser()
  if (verbose) print(stid(x))
  dimS <- ncdim_def( name="stid", units="number",vals=c(1:ns))
  dimT <- ncdim_def( name="time", units=paste("days since",torg), vals=time, calendar=calendar)
  
  if (verbose) print('Define variable')

### identify main variables in case of several variables in rda like tmin, tmax - must be implemented
  ## varname <- as.character(levels(factor(varid(x))))
  ## split variables
  ## for (i in 1:length(varname))
  ##    eval(parse(text=paste('var',i,'<-')))

  lon <- lon(x)
  lat <- lat(x)
  alt <- alt(x)
  
  ## Define variables 
  ##vars <- NA
  ##for (i in 1:na) {
  ##    if (la[i]==ns) 
  ##        eval(parse(text=paste("var",as.character(i)," <- ncvar_def(atts[",as.character(i),"],units=attrprec[",as.character(i),"],list(dimS),longname=atts[",as.character(i),"])",sep="")))
  ##}
  if (verbose) print(paste('create netCDF-file',file))
  ## browser()
     
  if (verbose) {str(y); print(summary(c(y)))}
  latid <- ncvar_def(name="lat",dim=list(dimS), units="degrees_north", missval=missval,longname="latitude", prec=prec,verbose=verbose)

  lonid <- ncvar_def(name="lon",dim=list(dimS), units="degrees_east", missval=missval,longname="longitude", prec=prec,verbose=verbose)
   altid <- ncvar_def(name="alt",dim=list(dimS), units="meters", missval=missval,longname="altitude", prec=prec,verbose=verbose)

  locid <- ncvar_def(name="loc",dim=list(dimS),units="strings",prec="char",longname="location",verbose=verbose)
  
  ncvar <- ncvar_def(name=varid(x)[1],dim=list(dimT,dimS), units=ifelse(unit(x)[1]=="°C", "degC",unit(x)[1]),longname=attr(x,'longname')[1], prec="float",compression=9,verbose=verbose)

  ## eval(parse(text= paste("vars <- list(ncvar,",paste(paste('var',c(1:na)[la==ns],sep=""),collapse=","),")",collapse="")))
  
  ncid <- nc_create(file,vars=list(ncvar,lonid,latid,altid,locid)) ## vars)
  ncvar_put( ncid, ncvar, y)
  ncatt_put( ncid, ncvar, 'add_offset',offs,prec='float')
  ncatt_put( ncid, ncvar, 'scale_factor',scalf,prec='float')
  ncatt_put( ncid, ncvar, 'missing_value',missval,prec='float')
  ncatt_put( ncid, ncvar, 'location',paste(loc(x),collapse=", "),prec='char')
  ncatt_put( ncid, ncvar, 'country',paste(cntr(x),collapse=", "),prec='character')
  browser()
  ncvar_put( ncid, lonid, attr(x,"longitude"))
  ncvar_put( ncid, latid, attr(x,"latitude"))
  ncvar_put( ncid, altid, attr(x,"altitude"))
  ncvar_put( ncid, locid, as.array(attr(x,"location")))
  
  ## for (i in 1:na) {
  ##    if (la[i]==ns)
  ##        eval(parse(text=paste("ncvar_put(ncid,var",as.character(i),",vals=attr(x,atts[",as.character(i),"]))",sep="")))
  ##    }
  
  ## if (verbose) print('save attributes')
  ## for (i in 1:na) {
  ##   if ( (la[i]==ns) & !is.na(attr(x,atts[i])[1]) ){
  ##    if (verbose) print(paste(atts[i],attrprec[i],sep=': '))
  ##     ncatt_put( ncid, ncvar, atts[i], attr(x,atts[i]), prec=attrprec[i] )
  ##   }
  ##}
 
  ##if (verbose) print('save global attributes')
  ##for (i in 1:na) {
  ##  if (la[i]!=ns) {
  ##    if (verbose) print(paste(atts[i],attrprec[i],sep=': '))
  ##    ncatt_put( ncid, 0, atts[i], attr(x,atts[i]), prec=attrprec[i] )
  ##  }
  ##}

  
  ## global attributes
  ncatt_put( ncid, 0, 'title', paste(levels(factor(attr(x,"info"))),collapse="/"))
  ncatt_put( ncid, 0, 'source', paste(levels(factor(attr(x,"source"))),collapse="/"))
  ncatt_put( ncid, 0, 'history', paste(unlist(attr(tmax,"history")),collapse="/"))
  ncatt_put( ncid, 0, 'references', paste(levels(factor(attr(tmax,"reference"))),collapse="/"))
  
  nc_close(ncid)
  if (verbose) print('close')
}
