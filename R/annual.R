#' Conversion to esd objects.
#' 
#' \code{annual} aggregates time series into annual values (e.g. means).
#' 
#' @aliases annual annual.zoo annual.default annual.dsensemble annual.station
#' annual.spell annual.field
#'
#' @seealso 
#' 
#' @param x an input object of, e.g., class 'station', 'field'
#' @param FUN a function, see \code{\link{aggregate.zoo}}
#' @param nmin Minimum number of data points (e.g. days or months) with valid
#' data accepted for annual estimate. NULL demands complete years.
#' @param format 'numeric' or 'character'
#' @param na.rm TRUE: ignore NA - see see \code{\link{mean}}
#' @param l length of window
#' @param it0 Value of specific year to include in the aggregated series
#' (it0=1970 and l=1 will give year 1970, 1975, 1980, ...)
#' 
#' @return Same class as x
#'
#' @seealso \code{\link{as.annual}} \code{\link{aggregate.station}} \code{\link{pentad}}
#' @keywords utilities
#'
#' @examples
#' data(ferder)
#' plot(annual(ferder,FUN="min"))
#' plot(annual(ferder,FUN="IQR",na.rm=TRUE))
#' 
#' data(bjornholt)
#' plot(annual(bjornholt,FUN="exceedance",fun="count"))
#' plot(annual(bjornholt,FUN="exceedance",fun="freq"))
#' plot(annual(bjornholt,FUN="exceedance"))
#' 
#' @export annual
annual <- function(x, ...) UseMethod("annual")

#' @export
annual.zoo <- function(x,FUN='mean',na.rm=TRUE,nmin=NULL, verbose=FALSE,...) {
  if (verbose) print("annual.zoo")
  if (inherits(x,'annual')) return(x)
  attr(x,'names') <- NULL
#  yr <- year(x)  REB: 08.09.2014
  class(x) <- 'zoo'
  ## Update the units for annual sums:
  if (FUN=='sum') {
    attr(x,'unit') <- sub('day','year',attr(x,'unit'))
    attr(x,'unit') <- sub('month','year',attr(x,'unit'))
    attr(x,'unit') <- sub('season','year',attr(x,'unit'))
  }
  
#  y <- aggregate(x,yr,FUN=match.fun(FUN),...,na.rm=na.rm)
  if ( (sum(is.element(names(formals(FUN)),'na.rm')==1)) |
       (sum(is.element(FUN,c('mean','min','max','sum','quantile')))>0) )
#    y <- aggregate(x,yr,FUN=FUN,...,na.rm=na.rm) else
#    y <- aggregate(x,yr,FUN=FUN,...)
      y <- aggregate(x,year,FUN=FUN,...,na.rm=na.rm)
  else
      y <- aggregate(x,year,FUN=FUN,...)
  ## replace infinite values by NA
  y[which(is.infinite(y))] <- NA
  ## Check minimum valid data points:
  d <- dim(y)
  if (!is.null(nmin)) {
    ## Need to account for both multiple and single series
    if (verbose) print(paste('nmin=',nmin))
    nok <- coredata(aggregate(x,year,FUN='nv'))
    ycd <- coredata(y) 
    ycd[nok < nmin] <-  NA
    ## If multivariate/matrix: reset dimensions
    if (!is.null(d)) dim(ycd) <- d
    coredata(y) <- ycd
  }
  attr(y,'dimnames') <- NULL
  invisible(y)
}


#' @export
annual.default <- function(x,FUN='mean',na.rm=TRUE, nmin=NULL,...,
                           threshold=NULL,regular=NULL,frequency=NULL,
                           verbose=FALSE) { ## 

  if (verbose) print('annual.default')
  
  ## Case when subsetting one specific season / in this case nmin =1
  
  ## If already annual, then return
  if (inherits(x,'annual')) return(x)
  ## Update the units for annual sums:
  if (FUN=='sum') {
    attr(x,'unit') <- sub('day','year',attr(x,'unit'))
    attr(x,'unit') <- sub('month','year',attr(x,'unit'))
    attr(x,'unit') <- sub('season','year',attr(x,'unit'))
  }
  
  ## This line to make the function more robust.
  if (length(grep('nmin',ls()))==0) nmin <- NULL
  
  if (inherits(FUN,'function')) FUN <- deparse(substitute(FUN)) # REB110314
  attr(x,'names') <- NULL
  yr <- year(x)
  nmo <- length(levels(factor(month(x))))
  d <- dim(x)
  if(is.null(d)) d <- c(length(x),1)
  #print(table(yr))
  
  YR <- as.numeric(rownames(table(yr)))
  nyr <- as.numeric(table(yr))

  # Need to accomodate for the possibility of more than one station series.
  if (inherits(x,'day')) {
    if (is.null(nmin)) nmin <- 30*nmo
  } else if (inherits(x,'month')) {
    if (is.null(nmin)) nmin <- 12
  } else if (inherits(x,'season')) {
    if (is.null(nmin)) nmin <-  length(levels(factor(month(x))))
  } else {
    nmin <- NA
  }
  if (verbose) {print(paste('nmin=',nmin)); print(class(x))}
  
  ## Convert x to a zoo-object:
  if (verbose) print('Number of valid data points')
  X <- zoo(coredata(x),order.by=index(x))
  attr(X,'units') <- unit(x)
  attr(X,'variable') <- varid(x)
  
  ## Check how many valid data points)
  nok <- aggregate(X,year,FUN='nv')

  if (FUN == 'sum') na.rm <- FALSE ## AM
  if (verbose) print(paste('aggregate: FUN=',FUN))

  if (verbose) str(X)
  
  if (sum(is.element(names(formals(FUN)),'threshold')==1)) {
    ## If threshold needed - set a default:
    if (is.null(threshold) & inherits(x,'station')) {
      threshold <- 1 ## AM added 20-05-2015
      if (verbose) print('Warning : threshold value not found and set to 1')
    }
    y <- aggregate(X,year,FUN=FUN,...,threshold=threshold) ## AM 20-05-2015
  } else if ((sum(is.element(names(formals(FUN)),'na.rm')==1)) |
           (sum(is.element(FUN,c('mean','min','max','sum','quantile')))>0)) {
    if (verbose) print('Function has na.rm-argument')
    y <- aggregate(X,year,FUN=FUN,...,na.rm=na.rm)
  } else {
    if (verbose) print('Function has no treshold nor na.rm arguments')
    y <- aggregate(X,year,FUN=FUN,...) # REB
  }
  y[!is.finite(y)] <- NA ## AM

  if (verbose) print('check for incomplete sampling')
  ## Flag the data with incomplete sampling as NA
  if (!is.na(nmin)) {
    ## Need to account for both multiple and single series
    if (verbose) {print(paste('nmin=',nmin)); print(nok)}
    #nok <- coredata(aggregate(x,year,FUN='nv'))
    ycd <- coredata(y) 
    ycd[coredata(nok) < nmin] <-  NA
    ## If multivariate/matrix: reset dimensions
    if (!is.null(dim(x))) dim(ycd) <- dim(y)
    y <- zoo(ycd,order.by=index(y))
    if (verbose) print(paste('mask',sum(nok < nmin),'years with nv <',nmin))
  }

  ## Copy the old attributes and reset as the original class:
  y <- attrcp(x,y,ignore="names")
  args <- list(...)
  if (verbose) print(names(args))

  ## Set appropriate units and variable names:
  if (verbose) print('Set appropriate units and variable names')
  if (FUN=="count")  {
    if (verbose) print("Count")
    attr(y,'unit') <-
      rep(paste("counts | X >",threshold," * ",attr(x,'unit')),d[2])
  } else if (FUN=="freq") {
    if (verbose) print("Frequency")
    attr(y,'variable') <- rep('f',d[2])
    attr(y,'unit') <- rep('fraction',d[2])
#    attr(y,'unit') <- rep(paste("frequency | X >",threshold," * ",attr(x,'unit')),d[2])
  } else if (FUN=="wetfreq") {
    if (verbose) print("Wet-day frequency")
    attr(y,'variable') <- rep('f[w]',d[2])
    attr(y,'unit') <- rep('fraction',d[2])
    attr(y,'longname')[] <- 'Wet-day frequency'
#    attr(y,'unit') <- rep(paste("frequency | X >",threshold," * ",attr(x,'unit')),d[2])
  } else if (FUN=="wetmean") {
    if (verbose) print("Wet-day mean")
    attr(y,'variable') <- rep('mu',d[2])
    attr(y,'longname')[] <- 'Wet-day mean precipitation'
    attr(y,'unit') <- rep('mm/day',d[2])
#    n <- count(X,threshold=threshold) # REB
#    n <- aggregate(X,year,FUN='count', threshold=threshold,...,
#                   regular = regular, frequency = frequency)
    n <- nok  # Not the count above threshold byut number of valid data points
    bad <- coredata(n)==0
    coredata(n)[bad] <- 1
    std.err <- 2*coredata(y)/sqrt(coredata(n)-1)
    std.err[bad] <- NA
    attributes(std.err) <- NULL
    dim(std.err) <- dim(y)
    attr(y,'standard.error') <- zoo(std.err,order.by=index(y))
  } else if (FUN=="mean") {
    if (verbose) print("mean")
    sigma <- aggregate(X, year, FUN='sd', ...,
                       regular = regular, frequency = frequency)
#    n <- count(x,threshold=threshold)
    n <- aggregate(X,year,FUN='count', threshold=threshold,...,
                   regular = regular, frequency = frequency)
    bad <- coredata(n)==0
    coredata(n)[bad] <- 1
    std.err <- 2*coredata(sigma)/sqrt(coredata(n)-1)
    std.err[bad] <- NA
    attributes(std.err) <- NULL
    dim(std.err) <- dim(sigma)
    attr(y,'standard.error') <- zoo(std.err,order.by=index(sigma))
  } else if (FUN=="HDD") {
    attr(y,'variable') <- rep('HDD',d[2])
    attr(y,'unit') <- rep('degree-days',d[2])
  } else if (FUN=="CDD") {
    attr(y,'variable') <- rep('CDD',d[2])
    attr(y,'unit') <- rep('degree-days',d[2])
  } else if (FUN=="GDD") {
    attr(y,'variable') <- rep('GDD',d[2])
    attr(y,'unit') <- rep('degree-days',d[2])
  } else attr(y,'unit') <- attr(x,'unit')

  attr(y,'history') <- history.stamp(x)
  class(y) <- class(x)
  class(y)[length(class(y))-1] <- "annual"
  if (class(y)[1]=="spell") class(y) <- class(y)[-1]
  #print(class(y)); print(class(x))
  invisible(y)
}

#' @export
annual.station <- function(x,FUN='mean',nmin=NULL,threshold=NULL,verbose=FALSE,...) {
  if (verbose) print('annual.station')
  attr(x,'names') <- NULL
  y <- annual.default(x,FUN=FUN,nmin=nmin,threshold=threshold,verbose=verbose,...)
  y[which(is.infinite(y))] <- NA
  invisible(y)
}

#' @export
annual.spell <- function(x,FUN='mean',nmin=0,threshold=NULL,verbose=FALSE,...) {
  attr(x,'names') <- NULL
  if ( (inherits(x,'mon'))  & is.null(nmin) ) {
    iy <- year(x)
    nmy <- as.numeric(table(iy))
    full <- nmy[nmy==12]
    x[is.element(iy,!full),] <- NA
    na.rm=FALSE
  }
  #y <- annual.default(x,FUN=match.fun(FUN),...)
  y <- annual.default(x,FUN=FUN,nmin=nmin,threshold=threshold,...) 
  invisible(y)
}

#' @export
annual.dsensemble <- function(x,FUN='mean',verbose=FALSE,...) {
  if (verbose) print("annual.dsensemble")
  clsx <- class(x)
  clss <- class(attr(x,'station'))
  if (!inherits(x,c('day','month','annual','season')))
      class(x) <- c(clsx[1],clss[2],clsx[2])
  if (inherits(x,'season'))
      y <- subset(x,it=0,verbose=verbose)
  else
      y <- annual.default(x,FUN=FUN,verbose=verbose,...)

  attr(y,'station') <- annual.station(attr(x,'station'),...)
  names(y) <- names(x)
  invisible(y)
}

#' @export
annual.field <- function(x,FUN='mean',na.rm=TRUE,nmin=NULL,verbose=FALSE, ...) {
  if (verbose) print('annual.field')
  attr(x,'names') <- NULL
  if (inherits(FUN,'function')) FUN <- deparse(substitute(FUN)) # REB110314
  yr <- year(x)
  cls <- class(x)
#  class(x) <- "zoo"
  if ( (inherits(x,'mon'))  & is.null(nmin) ) {
    iy <- year(x)
    nmy <- as.numeric(table(iy))
    full <- nmy[nmy==12]
    x[is.element(iy,!full),] <- NA
    na.rm=FALSE
  }
#  y <- aggregate(x,yr,FUN=match.fun(FUN),...,na.rm=na.rm)
#  y <- aggregate(x,yr,FUN=FUN,...,na.rm=na.rm)
  y <- annual.default(x,FUN=FUN,nmin=nmin,verbose=verbose,...) 
  y <- attrcp(x,y)
  attr(y,'history') <- history.stamp(x)
  attr(y,'dimensions') <- c(attr(x,'dimensions')[1:2],length(index(y)))
  class(y) <- cls
  # KMP 2017-10-25: time resolution is not always second in class vector, e.g., for an object 
  # of class [eof, comb, field, month, zoo] the following line will replace comb instead of month
  #class(y)[2] <- "annual"
  class(y)[length(cls)-1] <- "annual"
  y[which(is.infinite(y))] <- NA
  invisible(y)
}

#' @export
annual.eof <- function(x,FUN='mean',na.rm=TRUE,nmin=NULL,verbose=FALSE, ...) {
  if (verbose) print('annual.eof')
  attr(x,'names') <- NULL
  if (inherits(FUN,'function')) FUN <- deparse(substitute(FUN)) # REB110314
  yr <- year(x)
  cls <- class(x)
  #  class(x) <- "zoo"
  if ( (inherits(x,'mon'))  & is.null(nmin) ) {
    iy <- year(x)
    nmy <- as.numeric(table(iy))
    full <- nmy[nmy==12]
    x[is.element(iy,!full),] <- NA
    na.rm=FALSE
  }
  
  y <- annual.default(x,FUN=FUN,nmin=nmin,verbose=verbose,...) 
  y <- attrcp(x,y)
  attr(y,'history') <- history.stamp(x)
  attr(y,'dimensions') <- c(attr(x,'dimensions')[1:2],length(index(y)))
  class(y) <- cls
  # KMP 2017-10-25: time resolution is not always second in class vector, e.g., for an object 
  # of class [eof, comb, field, month, zoo] the following line will replace comb instead of month
  #class(y)[2] <- "annual"
  class(y)[length(cls)-1] <- "annual"
  y[which(is.infinite(y))] <- NA
  invisible(y)
}


