#' Correlation
#' 
#' Compute the correlation between field objects and station/field.
#' 
#' @aliases corfield corfield.default corfield.zoo corfield.field
#' corfield.field.station corfield.station corfield.eof corfield.trajectory
#'
#' @param x data object
#' @param y data object
#' @param plot TRUE: plot the results
#' @param use see \code{\link{cor}}
#' @param new see \code{link{map}}
#' @param ip index EOF pattern
#' @param \dots additional arguments
#'
#' @return Map of correlation
#'
#' @keywords manip
#'
#' @examples
#' x <- t2m.DNMI(lon=c(-40,30),lat=c(0,50))
#' y <- t2m.NorESM.M(lon=c(-40,30),lat=c(0,50))
#' r <- corfield(annual(x),annual(y))
#' 
#' data(Oslo)
#' t2m <- t2m.DNMI()
#' x <- subset(Oslo,it='january')
#' y <- subset(t2m,it='january')
#' r <- corfield(x,y)
#' 
#' @export corfield
corfield<-function(x,y,...) UseMethod("corfield")

#' @exportS3Method
#' @export corfield.default
corfield.default <- function(x,y,...) {
  cor(x,y)
}

#' @exportS3Method
#' @export corfield.zoo
corfield.zoo <- function(x,y,...,plot=TRUE,use='pairwise.complete.obs',verbose=FALSE,new=TRUE,
                         colbar=list(breaks=seq(-1,1,by=0.05),rev=TRUE)) {
  if (verbose) print("corfield.zoo:")
  stopifnot(inherits(y,'field'))

  # Keep track of which is an eof object and which is a station record:
  swapped <- FALSE
  if ( inherits(x,c("field"))) {
    yy <- x
    x <- y
    y <- yy
    swapped <- TRUE
  }

  # If the station data is daily, aggigate as monthly mean
  if ( inherits(y,"month") ) {
    x <- aggregate(x, as.yearmon, mean)
    print("Warning in corfield.station: aggregated y to monthly values.")
  }

  #print("HERE")
  if (verbose) {
    print(dim(x))
    print(dim(y))
    print(class(x))
    print(class(y))
  }
  yx <- merge(zoo(x),zoo(y),all=FALSE)
  #print("OK so far?")

  d <- dim(yx)
  r <- apply(coredata(yx[,2:d[2]]),2,cor,coredata(yx[,1]),use=use)
  #print(length(r))
  
  r <- attrcp(y,r)
  if (is.null(src(x))) attr(x,'source') <- NA
  if (is.null(src(y))) attr(y,'source') <- NA
  
  #print(names(attributes(x))); print(attr(y,'dimensions'))
  attr(r,'dimensions') <- attr(y,'dimensions')[1:2]
  attr(r,'longitude') <- attr(y,'longitude')
  attr(r,'latitude') <- attr(y,'latitude')
  attr(r,'variable') <- attr(y,'variable')
  attr(r,'source') <- attr(y,'source')
  attr(r,'unit') <- attr(y,'unit')
  attr(r,'time') <- range(index(y))
  attr(r,'n_data_points') <- d[1]
  class(r) <- 'corfield'

  #print("map")
  if (plot) map(r,colbar=colbar,new=new)
  invisible(r)
}

#' @exportS3Method
#' @export corfield.field
corfield.field <- function(x,y,...,plot=TRUE,use='pairwise.complete.obs',verbose=FALSE,new=TRUE,
                           colbar=list(breaks=seq(-1,1,by=0.05),rev=TRUE)) {
 
  if (verbose) {print('corfield.field'); print('field against field')}
  if (class(index(x))!=class(index(y))) warning('corfield.field: Different time index classes!')
  cor2s <- function(x,use,...) {
    n <- length(x);
    #print(n); print(length(x[1:(n/2)])); print(length(x[(n/2+1):n]))
    r <- cor(x[1:(n/2)],x[(n/2+1):n],use,...)
    invisible(r)
  }

  if (verbose) {print(dim(x)); print(attr(x,"dimensions"))
                print(dim(y)); print(attr(y,"dimensions"))}
 
  if (inherits(y,'station')) {
    r <- corfield.field.station(y,x,use=use,verbose=verbose,colbar=colbar,...)
    invisible(r)
  }
  if (verbose) print("synchonise")
  yx <- combine(x,y,dimension="space",all=FALSE,verbose=verbose)
  y <- yx$y; x <- yx$X
  if (verbose) {
    print(names(attributes(yx$X)))
    print(names(attributes(yx$y)))
    print("After combine")
    print(dim(x)); print(attr(x,"dimensions"))
    print(dim(y)); print(attr(y,"dimensions"))
  }
  
  if (sum(is.na(c( MATCH(attr(x,'longitude'), attr(y,'longitude')),
                   MATCH(attr(x,'latitude'), attr(y,'latitude')) )) ) > 0) {
    if (verbose) {print("regrid y:")
                  print(class(x)); print(class(y))
                  print(attr(x,'dimensions')); print(attr(y,'dimensions'))
                  print(dim(x)); print(dim(y))
                }

    #greenwich <- (max()*min() < 0) | (max()*min() < 0)
    x <- g2dl(x,greenwich=FALSE); y <- g2dl(y,greenwich=FALSE)
    #print(attr(x,'longitude')); print(attr(y,'longitude'))
    #print(attr(x,'latitude')); print(attr(y,'latitude'))

    lon.rng <- c(max(min(lon(x),lon(y),na.rm=TRUE)),
                 min(max(lon(x),lon(y),na.rm=TRUE)))
    lat.rng <- c(max(min(lat(x),lat(y),na.rm=TRUE)),
                 min(max(lat(x),lat(y),na.rm=TRUE)))
    #print("Region:"); print(c(lon.rng,lat.rng))
    x <- subset(x,is=list(lon.rng,lat.rng))
    y <- subset(y,is=list(lon.rng,lat.rng))
    y <- regrid(y,is=x)
  } else y <- regrid(y,is=x)
  
  #print(dim(y)); print(dim(x)); print(dim(rbind(coredata(x),coredata(y))))
  r <- apply(rbind(coredata(x),coredata(y)),2,cor2s,use=use,...)

  if (verbose) print(names(attributes(x)))
  r <- attrcp(x,r)
  
  attr(r,'dimensions') <- attr(x,'dimensions')[1:2]
  attr(r,'longitude') <- attr(x,'longitude')
  attr(r,'latitude') <- attr(x,'latitude')
  attr(r,'time') <- range(index(x))
  if (is.null(attr(x,'variable'))) attr(x,'variable') <- '.'
  if (is.null(attr(y,'variable'))) attr(y,'variable') <- '.'
  if (attr(x,'variable')[1]==attr(y,'variable')[1])
    attr(r,'variable') <- attr(x,'variable')[1] else
    attr(r,'variable') <- c(attr(x,'variable')[1], attr(y,'variable')[1])
  attr(r,'source') <- paste(attr(x,'source'),attr(y,'source'),sep="/")
  attr(r,'location') <- attr(y,'location')
  attr(r,'n_data_points') <- dim(x)[1]
  if (attr(x,'unit')[1]==attr(y,'unit')[1]) 
    attr(r,'unit') <- attr(x,'unit')[1] else
    attr(r,'unit') <- c(attr(x,'unit')[1],attr(y,'unit')[1])   
  attr(r,'history') <- history.stamp(x)
  attr(r,'history') <- history.stamp(y)
  class(r) <- 'corfield'
  if (plot) map(r,colbar=colbar,new=new)
  invisible(r)
}

#' @exportS3Method
#' @export corfield.field.station
corfield.field.station <- function(x,y,...,plot=TRUE,verbose=FALSE,new=TRUE,
                                   colbar=list(breaks=seq(-1,1,by=0.05),rev=TRUE),
                                   use='pairwise.complete.obs') {
  r <- corfield.station(y,x,plot=plot,verbose=verbose,use=use,colbar=colbar,new=new,...)
  invisible(r)
}

#' @exportS3Method
#' @export corfield.station
corfield.station <- function(x,y,...,plot=TRUE,verbose=FALSE,new=TRUE,
                             use='pairwise.complete.obs',
                             na.action='na.omit',colbar=list(breaks=seq(-1,1,by=0.05),rev=TRUE)) {
  if (verbose) print("corfield.station:")
  
  # Keep track of which is an eof object and which is a station record:
  nval <- function(x) sum(is.finite(x))
 
  swapped <- FALSE
  if ( inherits(y,c("station")) & inherits(x,c("field"))) {
    if (verbose) print('swap station and field')
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
  # ngood <- apply(coredata(y),2,nval)
  # ok <- ngood == length(index(y))
  # 
  # if (sum(ok)==0) {
  #   print('Problem with missing data. Try with the following argument')
  #   print('na.action="na.omit"')
  #   stop()
  # }
  # if (verbose) { print(paste(sum(ok),'good data points')); print(table(ngood)) }
  # yok <- y[,ok]
  # if (verbose) print(table(apply(coredata(yok),2,nval)))
  # #print(dim(yok)); image(coredata(yok))
  # cbind(zoo(x),yok,all=FALSE) -> xy
  # d <- dim(xy)
  # Y <- xy[,1]
  # X <- xy[,2:d[2]]
  # 
  # browser()
  # if (verbose) {print(d);print(dim(X)); print(length(y)); print(class(x)); print(class(y))}
  # rok <- apply(coredata(X),2,cor,coredata(Y),use=use,...)
  # r <- rep(NA,dim(y)[2])
  # r[ok] <- rok
  # if (verbose) {print(length(r))}
  if (verbose) {print('Before date mathing:');print(dim(y)); print(length(x))}
  y <- matchdate(y,it=x)
  x <- matchdate(x,it=y)
  if (verbose) {print('After date mathing:');print(dim(y)); print(length(x))}
  r <- apply(coredata(y),2,cor,coredata(x),use=use,...)
  r <- attrcp(y,r)

  if (verbose) {print(names(attributes(x))); print(attr(y,'dimensions'))}
  dim(r) <- attr(y,'dimensions')[1:2]
  attr(r,'dimensions') <- attr(y,'dimensions')[1:2]
  attr(r,'longitude') <- attr(y,'longitude')
  attr(r,'latitude') <- attr(y,'latitude')
  attr(r,'time') <- range(index(x))
  if (attr(x,'variable')[1]==attr(y,'variable')[1])
    attr(r,'variable') <- attr(x,'variable') else
    attr(r,'variable') <- c(varid(x)[1],varid(y)[1])
  attr(r,'source') <- paste(attr(x,'source'),attr(y,'source'),sep="/")
  attr(r,'location') <- attr(x,'location')
  attr(r,'x.longitude') <- attr(x,'longitude')
  attr(r,'x.latitude') <- attr(x,'latitude')
  attr(r,'x.altitude') <- attr(x,'altitude')
  attr(r,'stadion_id') <- attr(x,'station_id')
  attr(r,'country') <- attr(x,'country')
  attr(r,'n_data_points') <- dim(y)[1]
  if (attr(x,'unit')[1]==attr(y,'unit')[1]) 
    attr(r,'unit') <- attr(x,'unit')[1] else
    attr(r,'unit') <- c(unit(x)[1],unit(y)[1])
  attr(r,'history') <- history.stamp(x)
  attr(r,'history') <- history.stamp(y)
  class(r) <- 'corfield'

  #print("map")
  if (plot) map(r,verbose=verbose,colbar=colbar,new=new,...)
  invisible(r)
}

#' @exportS3Method
#' @export corfield.eof
corfield.eof <- function(x,y,...,ip=1,plot=TRUE,new=TRUE,
                         colbar=list(breaks=seq(-1,1,by=0.05),rev=TRUE),
                         use='pairwise.complete.obs',na.action='na.omit') {
  stopifnot(inherits(x,'eof'),inherits(y,'field'))
  z <- as.station(x[,ip],loc=paste('eof',ip),param='PC',unit='dimensionless')
  r <- corfield(z,y,plot=plot,use=use,na.action=na.action,colbar=colbar,new=new)
  invisible(r)
}

#' @exportS3Method
#' @export corfield.trajectory
corfield.trajectory <- function(x,y,...,it=NULL,is=NULL,param=NULL,FUN="count",
                                unit=NULL,longname=NULL,loc=NULL,
                                use="pairwise.complete.obs",
                                colbar=list(breaks=seq(-1,1,by=0.05),rev=TRUE),
                                method="pearson") {

  swapped <- FALSE
  if ( inherits(y,c("trajectory")) & inherits(x,c("field"))) {
    yy <- x
    x <- y
    y <- yy
    swapped <- TRUE
  }
  stopifnot(inherits(x,'trajectory'),inherits(y,'field'))

  x <- subset(x,it=it,is=is)
  y <- subset(y,it=it)
    
  xs <- trajectory2station(x,param=param,FUN=FUN,
                          unit=unit,longname=longname,loc=loc)
  if(any("season" %in% class(y))) {
    xs <- as.4seasons(xs)
  } else if (any("month" %in% class(y))) {
    xs <- as.monthly(xs)
  }
  xs <- subset(xs,it=y)
  r <- corfield(xs,y,use=use,method=method,colbar=colbar)
  invisible(r)
}
