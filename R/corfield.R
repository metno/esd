# R.E. Benestad, met.no, Oslo, Norway 12.04.2013
# rasmus.benestad@met.no
#------------------------------------------------------------------------


corfield<-function(x,y,...) UseMethod("corfield")


corfield.default <- function(x,y,...) {
  cor(x,y)
}


corfield.zoo <- function(x,y,plot=TRUE,use='pairwise.complete.obs',verbose=FALSE,...) {
  if (verbose) {print("corfield.zoo:"); print('station against field')}

  # Keep track of which is an eof object and which is a station record:
  swapped <- FALSE
  if ( inherits(x,c("field"))) {
    yy <- x
    x <- y
    y <- yy
    swapped <- TRUE
  }
  
  stopifnot(inherits(y,'field'))
  # If the station data is daily, aggigate as monthly mean
  if ( inherits(y,"month") ) {
    x <- aggregate(x, as.yearmon, mean)
    print("Warning in corfield.station: aggregated y to monthly values.")
  }

  #print("HERE")
  if (verbose) {print(dim(x)); print(dim(y)); print(class(x)); print(class(y))}
  yx <- merge(zoo(x),zoo(y),all=FALSE)
  #print("OK so far?")

  d <- dim(yx)
  r <- apply(coredata(yx[,2:d[2]]),2,cor,coredata(yx[,1]),use=use,...)
  #print(length(r))
  
  r <- attrcp(y,r)

  #print(names(attributes(x))); print(attr(y,'dimensions'))
  attr(r,'dimensions') <- attr(y,'dimensions')[1:2]
  attr(r,'longitude') <- attr(y,'longitude')
  attr(r,'latitude') <- attr(y,'latitude')
  attr(r,'variable') <- attr(y,'variable')
  attr(r,'source') <- attr(y,'source')
  attr(r,'unit') <- attr(y,'unit')
  attr(r,'time') <- range(index(y))
  class(r) <- 'corfield'

  #print("map")
  if (plot) map(r)
  return(r)
}

corfield.field <- function(x,y,plot=TRUE,use='pairwise.complete.obs',verbose=FALSE,...) {
 
  if (verbose) {print('corfield.field'); print('field against field')}
  cor2s <- function(x,use,...) {
    n <- length(x);
    #print(n); print(length(x[1:(n/2)])); print(length(x[(n/2+1):n]))
    r <- cor(x[1:(n/2)],x[(n/2+1):n],use,...)
    return(r)
  }

  #print("corfield.field")
  #print(dim(x)); print(attr(x,"dimensions"))
  #print(dim(y)); print(attr(y,"dimensions"))
 
  if (inherits(y,'station')) {
    r <- corfield.field.station(y,x,use=use,...)
    return(r)
  }
  #print("synchonise")
  yx <- combine(x,y,dimension="space",all=FALSE)
  #str(yx)
  #print("..."); print(summary(yx))
  y <- yx$y; x <- yx$X
  #print(dim(x)); print(dim(y))
  #print(names(attributes(yx$X)))
  #print(names(attributes(yx$y)))
  #print("After combine")
  #print(dim(x)); print(attr(x,"dimensions"))
  #print(dim(y)); print(attr(y,"dimensions"))
  #map(y)
  
  if (sum(is.na(c( MATCH(attr(x,'longitude'), attr(y,'longitude')),
                   MATCH(attr(x,'latitude'), attr(y,'latitude')) )) ) > 0) {
    #print("regrid y:")
    #print(class(x)); print(class(y))
    #print(attr(x,'dimensions')); print(attr(y,'dimensions'))
    #print(dim(x)); print(dim(y))
    #map(x); map(y)

    #greenwich <- (max()*min() < 0) | (max()*min() < 0)
    x <- g2dl(x,greenwich=FALSE); y <- g2dl(y,greenwich=FALSE)
    #print(attr(x,'longitude')); print(attr(y,'longitude'))
    #print(attr(x,'latitude')); print(attr(y,'latitude'))
    
    lon.rng <- c(max(min(attr(x,'longitude')), min(attr(y,'longitude'))),
                 min(max(attr(x,'longitude')), max(attr(y,'longitude'))))
    lat.rng <- c(max(min(attr(x,'latitude')), min(attr(y,'latitude'))),
                 min(max(attr(x,'latitude')), max(attr(y,'latitude'))))
    #print("Region:"); print(c(lon.rng,lat.rng))
    x <- subset(x,is=list(lon.rng,lat.rng))
    y <- subset(y,is=list(lon.rng,lat.rng))
    #print(dim(x)); print(dim(y))
    #print("HERE")
    y <- regrid(y,is=x)
    #print("HERE")
  } else y <- regrid(y,is=x)
  #map(y)
  
  #print(dim(y)); print(dim(x)); print(dim(rbind(coredata(x),coredata(y))))
  r <- apply(rbind(coredata(x),coredata(y)),2,cor2s,use=use,...)
  #print(length(r))

  #print(names(attributes(x)))
  #nattr <- softattr(x)
  #for (i in 1:length(nattr))
  #  attr(r,nattr[i]) <- attr(x,nattr[i])
  r <- attrcp(x,r)
  
  attr(r,'dimensions') <- attr(x,'dimensions')[1:2]
  attr(r,'longitude') <- attr(x,'longitude')
  attr(r,'latitude') <- attr(x,'latitude')
  attr(r,'time') <- range(index(x))
  if (attr(x,'variable')==attr(y,'variable'))
    attr(r,'variable') <- attr(x,'variable') else
    attr(r,'variable') <- attr(x,'variable')   # Need to finalise this...
#    attr(r,'variable') <- paste(attr(x,'variable'),attr(y,'variable'),sep="/")
  attr(r,'source') <- paste(attr(x,'source'),attr(y,'source'),sep="/")
  attr(r,'location') <- attr(y,'location')
  if (attr(x,'unit')==attr(y,'unit')) 
    attr(r,'unit') <- attr(x,'unit') else
    attr(r,'unit') <- attr(x,'unit')   # Need to finalise this...
#    attr(r,'unit') <- paste(attr(x,'unit'),attr(y,'unit'),sep="/")
  class(r) <- 'corfield'
  if (plot) map(r)
  return(r)
}


corfield.field.station <- function(x,y,plot=TRUE,
                                   use='pairwise.complete.obs',...) {
  r <- corfield.station(y,x,plot=plot,use=use,...)
  return(r)
}


corfield.station <- function(x,y,plot=TRUE,use='pairwise.complete.obs',na.action='na.omit',...) {
  #print("corfield.station:")
  # Keep track of which is an eof object and which is a station record:
  nval <- function(x) sum(is.finite(x))
  
  swapped <- FALSE
  if ( inherits(y,c("station")) & inherits(x,c("field"))) {
    yy <- x
    x <- y
    y <- yy
    swapped <- TRUE
  }
  
  stopifnot(inherits(x,'station'),inherits(y,'field'))
  # If the station data is daily, aggigate as monthly mean
  if ( (inherits(x,"day")) & (inherits(y,"month")) ) {
    x <- aggregate(x, as.yearmon, mean)
    print("Warning in corfield.station: aggregated y to monthly values.")
  }

#  arg <- list(...)
#  ina <- grep('na.action',names(arg))
#  if (length(ina)==1) {
  if (na.action=='na.omit') {
      ok <- is.finite(rowSums(y))
      y <- subset(y,it=range(year(y)[ok]))
    }
#    arg[[-ina]] -> list(...)
#  }
  
  #print("HERE")
  #print(length(x)); print(dim(y))
  #yx <- combine(x,y,all=FALSE)
    #y <- yx$y; x <- yx$X
  ngood <- apply(coredata(y),2,nval)
  ok <- ngood == length(index(y))
  if (sum(ok)==0) {
    print('Problem with missing data. Try with the following argument')
    print('na.action="na.omit"')
    stop()
  }
  #print(table(ngood)); print(sum(ok))
  yok <- y[,ok]
  #print(table(apply(coredata(yok),2,nval)))
  #print(dim(yok)); image(coredata(yok))
  merge(zoo(x),yok,all=FALSE) -> xy
  d <- dim(xy)
  Y <- xy[,1]
  X <- xy[,2:d[2]]
  
  #browser()
  #print(d);print(dim(X)); print(length(y)); print(class(x)); print(class(y))
  rok <- apply(coredata(X),2,cor,coredata(Y),use=use,...)
  r <- rep(NA,dim(y)[2])
  r[ok] <- rok
  #print(length(r))
  
  r <- attrcp(y,r)

  #print(names(attributes(x))); print(attr(y,'dimensions'))
  attr(r,'dimensions') <- attr(y,'dimensions')[1:2]
  attr(r,'longitude') <- attr(y,'longitude')
  attr(r,'latitude') <- attr(y,'latitude')
  attr(r,'time') <- range(index(X))
  if (attr(x,'variable')==attr(y,'variable'))
    attr(r,'variable') <- attr(x,'variable') else
    attr(r,'variable') <- attr(x,'variable')   # Need to finalise this...
  attr(r,'source') <- paste(attr(x,'source'),attr(y,'source'),sep="/")
  attr(r,'location') <- attr(x,'location')
  attr(r,'x.longitude') <- attr(x,'longitude')
  attr(r,'x.latitude') <- attr(x,'latitude')
  attr(r,'x.altitude') <- attr(x,'altitude')
  attr(r,'stadion_id') <- attr(x,'station_id')
  attr(r,'country') <- attr(x,'country')
  if (attr(x,'unit')==attr(y,'unit')) 
    attr(r,'unit') <- attr(x,'unit') else
    attr(r,'unit') <- attr(y,'unit') # Need to finalise this...
  class(r) <- 'corfield'

  #print("map")
  if (plot) map(r)
  invisible(r)
}

corfield.eof <- function(x,y,pattern=1,plot=TRUE,
                         use='pairwise.complete.obs',na.action='na.omit',...) {
  stopifnot(inherits(x,'eof'),inherits(y,'field'))
  z <- as.station(x[,pattern],loc=paste('eof',pattern),param='PC',unit='dimensionless')
  r <- corfield(z,y,plot=plot,use=use,na.action=na.action)
  invisible(r)
}

