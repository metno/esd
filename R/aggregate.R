#' aggregate
#' 
#' The aggregation functions are based on the S3 method for \code{zoo} objects,
#' but takes care of extra house keeping, such as attributes with meta data.
#'
#' \code{aggregate} calculates a time aggregate of an input object,
#' e.g, the mean seasonal cycle (if \code{by}=month and \code{FUN}="mean") or the
#' annual sum (if \code{by}=year and \code{FUN}="sum").
#'
#' \code{aggregate.area} is used for aggregating spatial statistics, such as
#' the global mean or the global area of some phenomenon.
#' \code{aggregate.size} is similar to \code{aggregate.area}, but returns the size statistics (square
#' meters) for individual events (defined as gridboxes touching each other).
#' 
#' @aliases aggregate aggregate.comb aggregate.field
#' @seealso aggregate.area aggregate.size aggregate.size.matrix aggregate.size.field
#' 
#' @param x A \code{\link{station}} object
#' @param by see \code{\link{aggregate.zoo}}
#' @param FUN see \code{\link{aggregate.zoo}}. Additional options: 'area','exceedance','lessthan'.
#' @param regular see \code{\link{aggregate.zoo}}
#' @param frequency see \code{\link{aggregate.zoo}}
#' @param na.rm TRUE: ignore NA - see see \code{\link{mean}}
#'
#' @return The call returns a station object
#'
#' @author R.E. Benestad
#' @keywords utilities
#' @examples
#' 
#' ## S3 method for class 'station'
#' data(Svalbard)
#' x <- aggregate(Svalbard, month, FUN='mean', na.rm=TRUE)
#' plot(x)
#'
#' ## S3 method for class 'field'
#' slp <- slp.DNMI()
#' y <- aggregate(slp, year, FUN='mean', na.rm=FALSE)
#'
#' ## Aggregate area
#' w <- aggregate.area(y)
#' plot(w)
#'
#' @export
aggregate.station <- function(x, by, FUN='mean', na.rm=TRUE, regular = NULL, ...,
                              frequency = NULL, verbose=FALSE, threshold=1) {

  if(verbose) print("aggregate.station")
  class(x) -> cls
  class(x) <- "zoo"

  if (by=='yearmon') {
    by <- as.Date(as.yearmon(index(x)))
  }
  
  if ( (sum(is.element(names(formals(FUN)),'na.rm')==1)) |
       (sum(is.element(FUN,c('mean','min','max','sum','quantile')))>0 ) ) {
    y <- aggregate(x, by, FUN, na.rm=na.rm, ...,
                   regular = regular, frequency = frequency)
  } else {
    y <- aggregate(x, by, FUN, ..., regular = regular, frequency = frequency)
  }

 # if (inherits(by[1],'character')) index(y) <- as.Date(index(y))
    
  if (class(index(y))=="Date") {
  dy <- day(y); mo <- month(y); yr <- year(y)
    if (dy[2] - dy[1] > 0) cls[length(cls) - 1] <- "day" else
    if (mo[2] - mo[1] == 1) cls[length(cls) - 1] <- "month" else
    if (mo[2] - mo[1] == 3) cls[length(cls) - 1] <- "season" else
    if (yr[2] - yr[1] > 0) cls[length(cls) - 1] <- "annual"
   } else
  if (class(index(y))=="yearmon") cls[length(cls) - 1] <- "month" else
  if (class(index(y))=="yearqtr") cls[length(cls) - 1] <- "qtr" else
  if (class(index(y))=="numeric") cls[length(cls) - 1] <- "annual" else
  if (class(index(y))=="character") cls[length(cls) - 1] <- "annual"
  if ( (length(index(y)) <= 12) & (class(index(y))=="numeric") & 
       (min(index(y)) >= 1) & (max(index(y)) <= 12) ) cls[length(cls) - 1] <- "seasonalcycle"
  class(y) <- cls
  y <- attrcp(x,y)

   #print(FUN)
  if (substr(FUN,1,4)=="count")  {
    #print("Count")
    attr(y,'unit') <- paste("counts | X >",threshold," * ",attr(x,'unit'))
  } else if (substr(FUN,1,4)=="freq") {
    #print("Frequency")
    attr(y,'variable') <- 'f'
    attr(y,'unit') <- paste("frequency | X >",threshold," * ",attr(x,'unit'))
  } else if (FUN=="wetfreq") {
    #print("Wet-day frequency")
    attr(y,'variable')[] <- 'f[w]'
    attr(y,'longname')[] <- 'Wet-day frequency'
    attr(y,'unit')[] <- paste("frequency | X >",threshold," * ",attr(x,'unit'))
  } else if (FUN=="wetmean") {
    #print("Wet-day mean")
    attr(y,'variable')[] <- 'mu'
    attr(y,'longname')[] <- 'Wet-day mean precipitation'
    attr(y,'unit')[] <- 'mm/day'
    n <- aggregate(x,by,FUN='nv', ...,
                   regular = regular, frequency = frequency)
    std.err <- 2*coredata(y)/sqrt(coredata(n)-1)
    attributes(std.err) <- NULL
    dim(std.err) <- dim(y)
    attr(y,'standard.error') <- zoo(std.err,order.by=index(y))
  } else if (FUN=="mean") {
    #print("Mean")
    sigma <- aggregate(x, by, FUN='sd', na.rm=na.rm, ...,
                       regular = regular, frequency = frequency)
    n <- aggregate(x, by, FUN='nv', na.rm=na.rm, ...,
                   regular = regular, frequency = frequency)
    #n <- count(x,threshold=threshold)
    std.err <- 2*coredata(sigma)/sqrt(coredata(n)-1)
    attributes(std.err) <- NULL
    dim(std.err) <- dim(y)
    #Finite population correction factor
    #http://en.wikipedia.org/wiki/Standard_error#Standard_error_of_the_mean
    #N <- length(n)
    #fpc <- sqrt((N-n)/(N-1))
    #std.err <- fpc*std.err
    attr(y,'standard.error') <- zoo(std.err,order.by=index(y))
  } else if (FUN=="HDD") {
    attr(y,'variable') <- 'HDD'
    attr(y,'unit') <- 'degree-days'
  } else if (FUN=="CDD") {
    attr(y,'variable') <- 'CDD'
    attr(y,'unit') <- 'degree-days'
  } else if (FUN=="GDD") {
    attr(y,'variable') <- 'GDD'
    attr(y,'unit') <- 'degree-days'
  } else attr(y,'unit') <- attr(x,'unit')

  attr(y,'history') <- history.stamp(y)
  return(y)
}

# Aggregate S3 method for a 'comb' object
#' @export
aggregate.comb <- function(x,by,FUN = 'mean', ...,
                              regular = NULL, frequency = NULL) {
  # Also need to apply the aggregation to the appended fields
  #if (verbose) print("aggregate.comb")
  #print(class(x))

  if (inherits(FUN,'function')) FUN <- deparse(substitute(FUN)) # REB140314
  if (deparse(substitute(by))=="year")
    by <- as.Date(strptime(paste(year(x),1,1,sep='-'),'%Y-%m-%d'))
  
  x <- aggregate.field(x,by=by,FUN=FUN, ...,
                       regular = regular, frequency = frequency)
  n <- attr(x,'n.apps')

  print("appended fields...")
  for (i in 1:n) {
    eval(parse(text=paste("z <- attr(x,'appendix.",i,"')",sep="")))
    #print(class(z)); print(by); print(index(z))
    z <- aggregate(z,by=by,FUN=FUN, ...,
                   regular = regular, frequency = frequency)
    print("update appended field")
    eval(parse(text=paste("z -> attr(x,'appendix.",i,"')",sep="")))
  }
  return(z)
}

# Aggregate S3 method for a field object
#' @export
aggregate.field <- function(x,by,FUN = 'mean', ...,
                            regular = NULL, frequency = NULL,
			    threshold=0, verbose=FALSE) {
  
  if (verbose) print('aggregate.field')
  class(x) -> cls
  #print(class(index(x)))
  class(x) <- "zoo"
  if (inherits(FUN,'function')) FUN <- deparse(substitute(FUN)) # REB140314
  #str(X); plot(X)
  #print("here...")
  #print(class(by))
  
  if (!is.list(by)) {
  # Temporal aggregation:
    #print("HERE")
    #print(deparse(substitute(by)))
    #clsy2 <- switch(deparse(substitute(by)),
    #                     "as.yearmon"="month",
    #                     "as.yearqtr"="quarter",
    #                     "as.annual"="annual",
    #                     "year"="annual",
    #                     "by" = "by")
    #if (is.null(clsy2)) clsy2 <- deparse(substitute(by))
    #print(clsy2)
    #if (deparse(substitute(by))[1]=="year") {
    #  ## KMP 2017-05-07: annual mean should have year as index, not date
    #  #by <- as.Date(strptime(paste(year(x),1,1,sep='-'),'%Y-%m-%d'))
    #  by <- year(x)
    #  index(x) <- year(x)
    #}
    ## REB - 'what do the following lines do?'year' changed to 'month' in the if-statement
    #if (deparse(substitute(by))[1]=="month") {
    #  print('fixed bug')
    #  by <- month(x)
    #  x <- as.data.frame(x)
    #  index(x) <- month(x)
    #}
    ## BER
    #print(deparse(substitute(by)))
    #print(class(x))
    #print(class(index(x)))
    #print('aggregate')
    ## y <- aggregate(x, by, match.fun(FUN), ...) ## AM quick fix replaced by
    y <- aggregate(x, by, FUN, ...)
    class(x) <- cls; 
    
    if (class(index(y))=="Date") {
      dy <- day(y); mo <- month(y); yr <- year(y)
      if (dy[2] - dy[1] > 0) cls[length(cls) - 1] <- "day" else
        if (mo[2] - mo[1] == 1) cls[length(cls) - 1] <- "month" else
          if (mo[2] - mo[1] == 3) cls[length(cls) - 1] <- "season" else
            if (yr[2] - yr[1] > 0) cls[length(cls) - 1] <- "annual"
    } else
      if (class(index(y))=="yearmon") cls[length(cls) - 1] <- "month" else
        if (class(index(y))=="yearqtr") cls[length(cls) - 1] <- "qtr" else
          if (class(index(y))=="numeric") cls[length(cls) - 1] <- "annual" else
            if (class(index(y))=="character") cls[length(cls) - 1] <- "annual"
    if ( (length(index(y)) <= 12) & (class(index(y))=="numeric") & 
         (min(index(y)) >= 1) & (max(index(y)) <= 12) ) cls[length(cls) - 1] <- "seasonalcycle"
    class(y) <- cls
    #class(y)[2] <- clsy2
    
    y <- attrcp(x,y)
    #nattr <- softattr(x)
    #for (i in 1:length(nattr))
    #  attr(y,nattr[i]) <- attr(x,nattr[i])
    attr(y,'dimensions') <- c(attr(x,'dimensions')[-3],length(index(y)))
    #attr(y,'history') <- c(attr(x,'history'),'aggregate')
    #attr(y,'date-stamp') <- date()
    #attr(y,'call') <- match.call()
    attr(y,'history') <- history.stamp(x)
    attr(y,'dimnames') <- NULL
    return(y)
  } else {
    #print("there")
  # spatial aggregation
    Lon <- by[[1]]; Lat=by[[2]]
    dx <- diff(Lon)[1]; dy <- diff(Lat)[1]
    t <- index(x)
    d <- attr(x,'dimensions')
    lon <- dx*floor(attr(x,'longitude')/dx)
    lat <- dy*floor(attr(x,'latitude')/dy)
    xy <- paste(rep(lon,d[2]),sort(rep(lat,d[1])),sep="/")
#    xy <- paste(rep(lat,d[1]),sort(rep(lon,d[2])),sep="/")
#    xy <- paste(sort(rep(lon,d[2])),rep(lat,d[1]),sep="/")
    lon <- as.numeric(rownames(table(lon)))
    lat <- as.numeric(rownames(table(lat)))
    D <-  c(length(lon),length(lat),length(t))
    print(paste("spatial aggregation:",d[1],"x",d[2]," ->",D[1],"x",D[2]))
    Z <- t(coredata(x)); attributes(Z) <- NULL
    dim(Z) <- c(d[1]*d[2],d[3]); rownames(Z) <- xy   
    ## z0 <- aggregate(Z,by=list(xy), match.fun(FUN),simplify=TRUE) ## AM 14-04-2015 replaced by
    z0 <- aggregate(Z,by=list(xy), FUN,simplify=TRUE)

    # The aggregate function rearranges the order of lon-lat:
    lonlat <- z0$Group.1
    ll <- strsplit(lonlat,split='/')
    lxy <- rep(NA,length(ll))
    for (i in 1:length(ll))
      lxy[i] <- as.numeric(ll[[i]][1]) + as.numeric(ll[[i]][2])*100
    #print(xy[1:100]); print(lonlat[1:100]); print(lonlat[order(lxy)][1:100])
    z <- as.matrix(z0[,2:(length(t)+1)])
    dim(z) <- c(D[1]*D[2],D[3])
    z <- z[order(lxy),]
    y <- zoo(t(as.matrix(z)),order.by=t)
    y <- attrcp(x,y)
    #nattr <- softattr(x)
    #for (i in 1:length(nattr))
    #  attr(y,nattr[i]) <- attr(x,nattr[i])
    attr(y,'dimensions') <- D
    attr(y,'latitude') <- lat
    attr(y,'longitude') <-lon
    #attr(y,'history') <- c(attr(x,'history'),'aggrgate')
    #attr(y,'date-stamp') <- date()
    #attr(y,'call') <- match.call()
    attr(y,'history') <- history.stamp(x)
    attr(y,'dimnames') <- NULL
    class(y) <- cls
    return(y)
    }
}
