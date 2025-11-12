#' Conversion to esd objects.
#' 
#' \code{as.annual} and \code{annual} aggregates time series into annual values (e.g. means).
#'
#' \code{as.monthly} aggregates time series into monthly values (e.g. means).
#'
#' \code{as.daily} aggregates time series into daily values (e.g. means).
#'
#' \code{as.4seasons} aggregates to four seasons ('djf': December-February, 'mam': March-May, 'jja': June-August, 'son': September-November)
#'
#' \code{as.seasons} aggregates to a user defined season with input arguments 'start' and 'end' giving the dates. To select march to september, use either start='03-01' and end='09-30' or start = 3 and end = 9.
#'
#' \code{as.OctMar} aggregates to the season October to March, which is the rainy season in parts of Africa. 
#'
#' @aliases as.annual as.annual.default as.annual.numeric as.annual.integer as.annual.yearqtr as.annual.station as.annual.spell
#' annual annual.zoo annual.default annual.dsensemble annual.station annual.spell annual.field annual.eof
#' as.monthly as.monthly.default as.monthly.station as.monthly.field
#' as.4seasons as.4seasons.default as.4seasons.day as.4seasons.station as.4seasons.spell as.4seasons.field as.4seasons.dsensemble as.seasons as.daily as.OctMar
#'
#' @seealso aggregate
#' 
#' @param x an input object of, e.g., class 'station', 'field'
#' @param FUN a function, see \code{\link{aggregate.zoo}}
#' @param nmin Minimum number of data points (e.g. days or months) with valid
#' data accepted for annual estimate. NULL demands complete years.
#' @param start makes it possible to estimate annual aggregated statistics that start into the year. 
#' E.g. start='Sep' starts the year on September 1st. Other allowed formats are start='MM-DD'. 
#' @param format 'numeric' or 'character'
#' @param na.rm a boolean; if TRUE, ignore NA - see see \code{\link{mean}}
#' @param verbose a boolean; if TRUE print information about progress
#' @param \dots additional argument
#' 
#' @return Same class as x
#'
#' @keywords utilities
#'
#' @examples
#' data(ferder)
#' plot(annual(ferder,FUN="min"))
#' plot(annual(ferder,FUN="IQR",na.rm=TRUE))
#' plot(annual(ferder))
#' lines(annual(ferder,start='Jul'))
#' 
#' data(bjornholt)
#' plot(as.4seasons(bjornholt,threshold=1,FUN="exceedance"))
#' 
#' @export annual
annual <- function(x, ...) UseMethod("annual")

#' @exportS3Method
#' @export annual.zoo
annual.zoo <- function(x,FUN='mean',na.rm=TRUE,nmin=NULL, start = NULL, verbose=FALSE,...) {
  if (verbose) print("annual.zoo")
  if (inherits(x,c('annual','year'))) return(x)
  attr(x,'names') <- NULL
  #  yr <- year(x)  REB: 08.09.2014
  class(x) <- 'zoo'
  ## If start specified, then use lag to make the series start to estimate annual aggregate
  ## starting from any random day in the year
  if (!is.null(start)) x <- shiftyear(x,start,verbose)
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


#' @exportS3Method
#' @export
annual.default <- function(x,FUN='mean',na.rm=TRUE, nmin=NULL,start=NULL,...,
                           minlen=NULL, threshold=NULL,regular=NULL,frequency=NULL,
                           verbose=FALSE) { ## 
  
  if (verbose) print(paste('annual.default',FUN))
  ## Case when subsetting one specific season / in this case nmin =1
  if (is.null(minlen) & !is.null(nmin)) minlen <- nmin 
  ## If already annual, then return
  if (inherits(x,'annual')) return(x)
  ## If start specified, then use lag to make the series start to estimate annual aggregate
  ## starting from any random day in the year
  if (!is.null(start)) x <- shiftyear(x,start,verbose)
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
  
  # Need to accommodate for the possibility of more than one station series.
  if (inherits(x,'day')) {
    if (is.null(nmin)) nmin <- 30*nmo
  } else if (inherits(x,'month')) {
    if (is.null(nmin)) nmin <- 12
  } else if (inherits(x,'season')) {
    if (is.null(nmin)) nmin <-  length(levels(factor(month(x))))
  } else {
    nmin <- NULL
  }
  if (verbose) {print(paste('nmin=',nmin)); print(class(x))}
  
  ## Convert x to a zoo-object:
  if (verbose) print('Number of valid data points')
  X <- zoo(coredata(x),order.by=index(x))
  attr(X,'units') <- esd::unit(x)
  attr(X,'variable') <- varid(x)
  
  ## Check how many valid data points)
  nok <- aggregate(X,year,FUN='nv')
  
  ## Check how many data points, valid or not
  nlen <- aggregate(X,year,FUN=length)
  
  ## REB 2024-08-27: The argument 'nmin' takes care of NAs or by setting 'na.rm=FALSE' as argument.
  #if (FUN == 'sum') na.rm <- FALSE ## AM
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
    if (verbose) print(paste(FUN,'has na.rm-argument:',na.rm))
    y <- aggregate(X,year,FUN=FUN,...,na.rm=na.rm)
  } else {
    if (verbose) print('Function has no treshold nor na.rm arguments')
    y <- aggregate(X,year,FUN=FUN,...) # REB
  }
  y[!is.finite(y)] <- NA ## AM
  if (verbose) print(y)
  
  if (verbose) print('check for incomplete sampling')
  ## Need to account for both multiple and single series
  ycd <- coredata(y)

  ## Mask values with few data points
  if(is.null(minlen)) {
    if(inherits(x, "day")) minlen <- 360 else 
      if (inherits(x, "month")) minlen <- 12 else 
        if (inherits(x, "season")) minlen <- 4 else minlen <- 1
  }
  few <- nlen < minlen
  if (verbose) cat('Years with less than',minlen,'data points',sum(few))
  ycd[few] <-  NA
  if (verbose) print(paste('mask',sum(few),'years with length <',minlen))
  ## Mask values with few valid data points
  if (!is.na(nmin)) {
    if (verbose) {print(paste('nmin=',nmin)); print(nok)}
    ycd[coredata(nok) < nmin] <-  NA
    if (verbose) print(paste('mask',sum(nok < nmin),'years with nv <',nmin))
  }
  ## If multivariate/matrix: reset dimensions
  if (!is.null(dim(x))) dim(ycd) <- dim(y)
  y <- zoo(ycd,order.by=index(y))
  
  ## Check if the series is gappy:
  if (verbose) print('Fill in gaps')
  it <- seq(min(index(y)),max(index(y)),by=1)
  if (verbose) print(c(range(it),length(it),sum(is.finite(y))))
  if(is.null(dim(y))) {
    z <- zoo(rep(NA,length(it)),order.by=it)
    z[is.element(it,index(y))] <- y
  } else {
    z <- zoo(matrix(rep(NA,ncol(y)*length(it)),
                    ncol=ncol(y), nrow=length(it)),
             order.by=it)
    z[is.element(it,index(y)),] <- y
  }
  y <- z; rm('z')
  
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
    n <- nok  # Not the count above threshold, but number of valid data points
    bad <- n<nmin
    #bad <- coredata(n)==0
    #coredata(n)[bad] <- 1
    if (verbose) cat('n=',n,'\n')
    denominator <- sqrt(coredata(n)-1)
    denominator[denominator==0] <- NA
    if (length(y)==length(denominator)) { 
      std.err <- try(coredata(y)/denominator)
      if (!inherits(std.err,'try-error')) { 
        std.err[bad] <- NA
        attributes(std.err) <- NULL
        dim(std.err) <- dim(y)
        attr(y,'standard.error') <- zoo(std.err,order.by=index(y))
      }
    }
  } else if (FUN=="mean") {
    if (verbose) print("mean")
    sigma <- aggregate(X, year, FUN='sd', ...,
                       regular = regular, frequency = frequency)
    ## KMP 2022-12-09: threshold is not set by default
    ## n should be the numnber of valid data points, not count above threshold
    n <- nok
    bad <- n<nmin
    #    n <- count(x,threshold=threshold)
    #n <- aggregate(X,year,FUN='count', threshold=threshold,...,
    #               regular = regular, frequency = frequency)
    #bad <- coredata(n)==0
    #coredata(n)[bad] <- 1
    std.err <- try(2*coredata(sigma)/sqrt(coredata(n)-1))
    if (!inherits(std.err,'try-error')) { 
      std.err[bad] <- NA
      attributes(std.err) <- NULL
      dim(std.err) <- dim(sigma)
    }
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

#' @exportS3Method
#' @export
annual.station <- function(x,FUN='mean',nmin=NULL,start=NULL,threshold=NULL,verbose=FALSE,...) {
  if (verbose) print(paste('annual.station',FUN))
  attr(x,'names') <- NULL
  y <- annual.default(x,FUN=FUN,nmin=nmin,start=start,threshold=threshold,verbose=verbose,...)
  y[which(is.infinite(y))] <- NA
  invisible(y)
}

#' @exportS3Method
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

#' @exportS3Method
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

#' @exportS3Method
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

#' @exportS3Method
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

#' @export as.annual
as.annual <- function(x, ...) UseMethod("as.annual")

#' @exportS3Method
#' @export
as.annual.default <- function(x, ...) annual(x, ...)

#' @exportS3Method
#' @export
as.annual.numeric <- function(x, ...) annual(x, ...)

#' @exportS3Method
#' @export
as.annual.integer <- function(x, ...) structure(x, class = "annual")

#' @exportS3Method
#' @export
as.annual.yearqtr <- function(x, frac = 0, ...) {
  if (frac == 0) annual(as.numeric(x)) else
    as.annual(as.Date(x, frac = frac), ...)
}

#' @exportS3Method
#' @export
as.annual.station <- function(x, ...) annual.station(x,...)

#' @exportS3Method
#' @export
as.annual.spell <- function(x, ...) annual.spell(x,...)

#' @export
as.monthly <- function(x,...) UseMethod("as.monthly")

yyyymm <- function(x) ym <- as.Date(paste(year(x),month(x),'01',sep='-'))

#' @exportS3Method
#' @export
as.monthly.default <- function(x,...) {
  y <- aggregate(x,by=yyyymm,...)
  return(y)
}

#' @exportS3Method
#' @export
as.monthly.field <- function(x,FUN='mean',...) {
  if (inherits(x,'month')) return(x)
  y <- aggregate(as.zoo(x), yyyymm, #function(tt) as.Date(as.yearmon(tt)),
                 FUN=FUN,...)
  y <- attrcp(x,y)
  attr(y,"dimensions") <- c(attr(x,"dimensions")[1:2],length(index(y)))
  attr(y,'history') <- history.stamp(x)
  class(y) <- class(x)
  class(y)[2] <- "month" 
  return(y)
}


#' @exportS3Method
#' @export
## This is a dublicate of that in as.R
as.monthly.station <- function (x, FUN = "mean", ...) {
  y <- aggregate(zoo(x), yyyymm, #function(tt) as.Date(as.yearmon(tt)), 
                 FUN = FUN, ...)
  y <- attrcp(x, y)
  attr(y, "history") <- history.stamp(x)
  class(y) <- class(x)
  class(y)[2] <- "month"
  return(y)
}

#' @export
as.daily <- function(x,...) UseMethod("as.daily")

yyyymmdd <- function(x) ymd <- as.Date(paste(year(x),month(x),day(x),sep='-'))

#' @exportS3Method
#' @export
as.daily.default <- function(x,...) {
  y <- aggregate(x,by=yyyymmdd,...)
  return(y)
}

#' @exportS3Method
#' @export
as.daily.field <- function(x,FUN='mean',...) {
  if (inherits(x,'month')) return(x)
  y <- aggregate(as.zoo(x), yyyymmdd, #function(tt) as.Date(as.yearmon(tt)),
                 FUN=FUN,...)
  y <- attrcp(x,y)
  attr(y,"dimensions") <- c(attr(x,"dimensions")[1:2],length(index(y)))
  attr(y,'history') <- history.stamp(x)
  class(y) <- class(x)
  class(y)[2] <- "day" 
  return(y)
}

#' @exportS3Method
#' @export
## This is a dublicate of that in as.R
as.daily.station <- function (x, FUN = "mean", ...) {
  y <- aggregate(zoo(x), yyyymmdd, #function(tt) as.Date(as.yearmon(tt)), 
                 FUN = FUN, ...)
  y <- attrcp(x, y)
  attr(y, "history") <- history.stamp(x)
  class(y) <- class(x)
  class(y)[2] <- "day"
  return(y)
}

#' @export
as.4seasons <- function(x,...) UseMethod("as.4seasons")

#' @exportS3Method
#' @export
as.4seasons.default <- function(x,...,FUN='mean',slow=FALSE,verbose=FALSE,nmin=NULL) {
  if(verbose) print('as.4seasons.default')
  if (inherits(x,c('season','seasonal'))) return(x)
  attr(x,'names') <- NULL
  d <- dim(coredata(x))
  #print(d)
  if (is.null(d)) d <- c(length(x),1)
  if (!slow) {
    if (inherits(x,"month")) {
      if ( (is.null(d)) | (d[2]==1) ) {
        X <- c(NA,coredata(x)[1:length(x)-1]) # shift the coredata by 1 to start on December. This works only for monthly data !!!  
      } else {
        # shift the coredata by 1 to start on December. slow...
        if (verbose) print('Shift series')
        X <- rbind(rep(NA,d[2],1),coredata(x)[1:d[1]-1,])
      }
      if (verbose) print(dim(X))
      X <- zoo(X,order.by=index(x))
      ##yrseas <- fourseasons(ix)
      ##print(yrseas)
      if (verbose) print('aggregate')
      #print(names(list(...)))
      yq <- function(t) as.yearqtr(year(t) + 0.25*floor((month(t)-1)/3))
      y <- aggregate(x=as.zoo(X),by=yq,#as.yearqtr,
                     FUN=match.fun(FUN),...)
      # convert yearqtr to yearmon
      ## Remove season values with less than nmin data points
      if(is.null(nmin)) nmin <- 3
      nd <- aggregate(x=as.zoo(X),by=yq,FUN=nv)
      ok <- nd >= nmin  
      coredata(y)[!ok] <- NA
      if (verbose) print('define y')
      y <- zoo(x=y,order.by=as.Date(as.yearmon(index(y))))
    } else y <- as.4seasons.day(x,FUN=FUN,nmin=nmin,verbose=verbose,...)
    #y <- as.4seasons.day(x,FUN=match.fun(FUN),...)
    
    #print(dim(y))
    ok <- length(index(y))
    #print(summary(c(coredata(y))))
  } else {
    yr <- sort(rep(as.integer(rownames(table(year(x)))),4))
    n <- length(yr)
    
    q <- rbind(c(12,1,2),3:5,6:8,9:11)
    X <- matrix(rep(NA,n*d[2]),n,d[2]) 
    t <-rep(NA,n)
    
    if (verbose) print("start loop")
    for (i in 1:n) {
      iq <- (i-1) %% 4 + 1
      if (iq == 1) 
        ii <- (is.element(year(x),yr[i]) &
                 is.element(month(x),c(1,2))) |
          (is.element(year(x),yr[i]-1) &
             is.element(month(x),12))  else
               ii <-  is.element(year(x),yr[i]) &
                 is.element(month(x),q[iq,])
             if (d[2]==1) cline <- paste(FUN,"(coredata(x[ii,]),...)",sep="") else
               cline <- paste("apply(coredata(x[ii,]),2,",FUN,", ...)",sep="")
             #print(cline)
             if ( (inherits(x,'day')) & (sum(ii)>=85) |
                  (inherits(x,'month')) & (sum(ii)>=3) ) {
               X[i,] <- eval(parse(text=cline))
               t[i] <- yr[i] + (iq-1)/4
               print(c(i,yr[i],round(mean(X[i,]),2),iq,sum(ii),round(X[i],2),t[i]))
             }
    } 
    if (verbose) print("end loop")
    ok <- is.finite(rowMeans(X,na.rm=TRUE)) & is.finite(t)
    if (verbose) print(table(as.yearqtr(t[ok])))
    #print(summary(c(X)))
    y <- zoo(X[ok,],order.by=as.Date(as.yearqtr(t[ok])))
    #names(y) <- FUN
  }
  if (!is.null(dim(y))) {
    ok <- is.finite(rowMeans(y,na.rm=TRUE))
    y <- y[ok,]
  }
  if (verbose) print('attributes')
  y <- attrcp(x,y)
  attr(y,'history') <- history.stamp(x)
  attr(y,'season.interval') <- "4seasons"
  if (inherits(x,'field'))
    attr(y,'dimensions') <- c(attr(x,'dimensions')[1:2],sum(ok))
  class(y) <- class(x)
  class(y)[2] <- "season"
  #plot(y)
  return(y) 
}

#' @exportS3Method
#' @export
as.4seasons.day <- function(x,...,FUN='mean',na.rm=TRUE,dateindex=TRUE,nmin=85,verbose=FALSE) {
  if(verbose) print('as.4seasons.day')
  if (inherits(x,'month')) nmin <- 3 # AM 06-07-2015
  attr(x,'names') <- NULL  
  t <- index(x)
  year <- year(t) #as.numeric(format(t,'%Y'))
  month <- month(t) #as.numeric(format(t,'%m'))
  day <- day(t) # as.numeric(format(t,'%d'))
  #shift the time stamps by one month, sneaking December into the subsequent year
  month <- month + 1
  dec <- is.element(month,13)
  year[dec] <- year[dec] + 1
  month[dec] <- 1
  # Change the day to avoid warning that the calendar is wrong (e.g. due to
  # too many days in February). Since the data is aggregated, the exact day
  # in the month doesn't matter here.
  hour <- 12*(day - 2*trunc(day/2))
  day <- trunc(day/2) + 1
  #print(table(year)); print(table(month));
  #print(table(day)); print(table(hour));   
  #tshifted <- as.Date(paste(year,month,day,sep="-"))
  tshifted <-  ISOdate(year=year,month=month,day=day,hour=hour)
  #print(summary(tshifted))
  X <- zoo(coredata(x),order.by=tshifted)
  
  # Test for the presens of 'na.rm' in argument list - this is a crude fix and not a
  # very satisfactory one. Fails for FUN==primitive function.
  if (is.function(FUN)) test.na.rm <- FALSE else
    test.na.rm <- (sum(is.element(FUN,c('mean','min','max','sum','quantile')))>0)
  if ( (sum(is.element(names(formals(FUN)),'na.rm')==1)) | (test.na.rm) )
    y <- aggregate(X,as.yearqtr,FUN=match.fun(FUN),...,na.rm=na.rm) else
      y <- aggregate(X,as.yearqtr,FUN=match.fun(FUN),...)
  
  # Set to missing for seasons with small data samples:
  nd <- aggregate(X,as.yearqtr,FUN=nv)
  ok <- nd >= nmin
  coredata(y)[!ok] <- NA
  # dateindex: convert "1775 Q1" to "1775-01-01"
  if (dateindex) {
    y <- zoo(coredata(y),order.by=as.Date(index(y)))
  }
  unit <- attr(y,'unit')
  y <- attrcp(x,y,ignore=c("unit","names"))
  unit -> attr(y,'unit')
  #str(y); print(unit)
  attr(y,'history') <- history.stamp(x)
  attr(y,'season.interval') <- "4seasons"
  class(y) <- class(x)
  class(y)[2] <- "season"
  invisible(y)
}

#' @exportS3Method
#' @export
as.4seasons.station <- function(x,...,FUN='mean') {
  #print('as.4seasons.station')
  y <- as.4seasons.default(x,FUN=FUN,...)
  #  y <- attrcp(x,y)
  #  attr(y,'history') <- history.stamp(x)
  #  class(y) <- class(x)
  #  class(y) <- gsub("month","season",class(x))
  return(y) 
}

#' @exportS3Method
#' @export
as.4seasons.spell <- function(x,...,FUN='mean') {
  y <- as.4seasons.default(as.station(x),FUN=FUN,...)
  #  y <- attrcp(x,y)
  #  attr(y,'history') <- history.stamp(x)
  #  class(y) <- class(x)
  #  class(y) <- gsub("month","season",class(x))
  return(y) 
}

#' @exportS3Method
#' @export
as.4seasons.field <- function(x,...,FUN='mean',verbose=FALSE) {
  if(verbose) print("as.4seasons.field")
  d <- attr(x,"dimensions")
  y <- as.4seasons.default(x,...,FUN=FUN,verbose=verbose)
  y <- attrcp(x,y)
  attr(y,'history') <- history.stamp(x)
  ## Update dimensions
  attr(y, "dimensions") <- c(d[1],d[2],dim(y)[1])
  class(y) <- class(x)
  class(y) <- gsub("month","season",class(x))
  return(y)
}

#' @exportS3Method
#' @export
as.4seasons.dsensemble <- function(x,...,FUN='mean') {
  cls <- class(x)
  class(x) <- c("station",cls[2],"zoo") ## AM 06-07-2015 Quick fix here, time step added into the class of x
  attrx <- attributes(x)
  y <- as.4seasons.station(x,FUN=FUN,...)
  ##attributes(y) <- attrx
  
  y <- attrcp(x,y)
  attr(y,"station") <- as.4seasons.station(attr(x,"station"))
  
  attr(y,'history') <- history.stamp(x)
  class(y) <- c("dsensemble","season","zoo")
  return(y)
}

# do not export - local function only used in as.seasons
leapdate <- function(years="2000", dates="02-29") {
  yeardate <- paste(years, dates, sep="-")
  nok <- !leapyear(years) & dates=="02-29"
  yeardate[nok] <- paste(years[nok], "02-28", sep="-")
  return(as.Date(yeardate))
}

# Not to confuse with season
# This function extracts a given seasonal interval and aggregates a given statistic
#' @export
as.seasons <- function(x,start='01-01',end='12-31',nmin=NULL,FUN='mean',
                       na.rm=TRUE,verbose=FALSE,...) {
  if(verbose) print("as.seasons")
  yrs <- year(x); d <- dim(x)
  # ns = number of stations
  if (is.null(d)) ns <- 1 else ns <- d[2]
  years <- as.numeric(rownames(table(yrs))); n <- length(years)
  y <- matrix(rep(NA,n*ns),n,ns); k <- y
  if(is.numeric(start) & is.numeric(end)) {
    if(verbose) print("input 'start' and 'end' are likely months")
    if(start>=10) start <- paste0(start,"-01") else start <- paste0("0",start,"-01")
    if(end==2) {
      end <- "02-29"
    } else {
      if(end<10) end <- paste0("0",end)
      ym <- as.yearmon(paste0("2020-",end,"-01"))
      days <- as.Date(ym, frac = 1) - as.Date(ym) + 1
      end <- paste0(end, "-", days)
    }
  } else {
    if(verbose) print("input 'start' and 'end' are likely dates (mm-dd)")
  }
  
  start.1 <- as.numeric(leapdate(years[1], start))
  end.1 <- as.numeric(leapdate(years[1], end))
  if (start.1 > end.1) twoyears <- 1 else twoyears <- 0
  if(verbose) print(paste(start, end))
  
  for (i in 1:n) {
    if(verbose) print(paste("Aggregate for year",years[i]))
    z <- coredata(window(x, start=leapdate(years[i], start),
                         end=leapdate(years[i]+twoyears, end)))
    k[i,] <- apply(matrix(z,ceiling(length(z)/ns),ns),2,nv)
    y[i,] <- apply(matrix(z,ceiling(length(z)/ns),ns),2,FUN, na.rm=na.rm, ...)
  }
  if(!is.null(nmin)) y[ k < nmin ] <- NA
  y <- zoo(y, order.by = as.Date(paste(years, start, sep='-')))
  y <- attrcp(x, y)
  attr(y,'history') <- history.stamp(x)
  if (twoyears==0) {
    attr(y,'season.interval') <- paste(start,'to',end)
  } else {
    attr(y,'season.interval') <- paste(start,'to',end,'the following year')
  }
  attr(y, 'n.valid') <- k
  class(y) <- class(x)
  class(y)[2] <- "season"
  #class(y)[2] <- "annual"
  return(y)
}

## @RasmusBenestad, 2021-04-12
## Summary statistics for the rainy season in Africa: October to March.
## Function that processes the data - it assumes daily precipitation
#' @export
as.OctMar <- function(x,FUN='sum',nmin=90,plot=FALSE,verbose=FALSE) {
  ## x is a station object from esd
  ## Step 1: Oct-Dec from one year
  if (verbose) {print('as.OctMar'); print(class(x))}
  if (!is.precip(x)) warning("as.OctMar was designed for daily rainfall data, but that's OK")
  OctDec <- subset(x,it=month.abb[10:12])
  ## Step 2: Jan - Mar from another year
  JanMar <- subset(x,it=month.abb[1:3])
  ## Now we need to take the sum over the months using 'annual'
  OctDec <- annual(OctDec,FUN=FUN,nmin=nmin)
  JanMar <- annual(JanMar,FUN=FUN,nmin=nmin)
  ## R will add the data with corresponding years
  ## index() controls the time information
  index(OctDec) <- year(OctDec)
  index(JanMar) <- year(JanMar)-1
  ## Check: 
  #print(range(index(OctDec)))
  #print(range(index(JanMar)))
  ## Add the Jan-Mar sum with the previous years Oct-Dec sum
  if (FUN=='sum') OctMar <- JanMar + OctDec else
    OctMar <- 0.5*(JanMar + OctDec)
  #print(range(index(OctMar)))
  OctMar <- attrcp(JanMar,OctMar)
  class(OctMar) <- class(JanMar)
  attr(OctMar,'season.interval') <- "10-01 to 03-31 the following year"
  
  if (plot==TRUE) {
    ## if the plot argument == TRUE, then do this:
    y <- as.matrix(OctMar)
    ## Here we set the y-axix label and wexpress y in % if FUN=='wetfreq'
    if (FUN=='wetfreq') {
      y <- y * 100
      ylab='%'
    } else if (is.precip(x)) ylab <- 'mm'
    ## Assign row names: the year and following year
    rownames(y) <- paste(year(OctMar),year(OctMar)+1,sep='-')
    ## Make a graph with bars
    #plot(y[,1],type='l') ## If you want to plot lines
    ## Loop over all the stations:
    ns <- length((loc(x)))
    ## Loop over stations
    for (i in 1:ns) {
      par(las=2,cex.axis=0.5)
      ## This line is to make a more sensible title for most people:
      fun <- switch(FUN,'mean'='mean','sum'='Total rainfall',
                    'wetfreq'='wet day frequency',
                    'wetmean'='mean intensity')
      ## Here we make a title text which gives additional information
      ## about the plot and the data record.
      main=paste('Rainy season for',loc(x)[i],':',fun,
                 'trend=',round(100*trend.coef(y[,i])/mean(y[,i],na.rm=TRUE),1),'%/decade')
      ## Show a bar plot for the data, but only when it's not missing
      barplot(y[is.finite(y[,i]),i],main=main,col='blue',
              ylab=ylab)
      ## Add a line showing the median value:
      lines(c(1,length(y)),rep(median(y[,i],na.rm=TRUE),2),lty=2)
      lines(c(1,length(y)),rep(mean(y[,i],na.rm=TRUE),2),lty=3)
      lines(trend(y[,i]),col='red',lwd=2)
      grid()
    }
  }
  return(OctMar)
}

## @RasmusBenestad, 2021-04-12
## A function to assist more flexible definition of a 'year' that may start any day between
## January 1st and December 31st. The argument start may be a calendar name in {month.abb} or
## defined 
#' @export
shiftyear <- function(x,start,verbose=FALSE) {
  x0 <- x
  if (verbose) print(paste('shiftyear: The year start is',start))
  if (inherits(x,'season')) warning('annual: does not support start for seasonal aggregates')
  imon <- grep(start,month.abb,ignore.case = TRUE)
  ## If start is a month ('Jan',...'Dec)
  if (length(imon)>0) {
    ## If daily data, need to find he number of days into the year
    if (inherits(x,'day')) {
      yyyy1 <- year(x)[1]
      t1 <- as.Date(paste(yyyy1,'01-01',sep='-'))
      t2 <- as.Date(paste(yyyy1,imon,'01',sep='-'))
      imon <- t2 - t1
    } else imon <- imon - 1
    if (verbose) print(imon)
    x <- lag(x,imon) 
  } else if (is.character(start)) { 
    if (nchar(start)==5) {
      ## if start is in format 'MM-DD'
      yyyy1 <- year(x)[1]
      t1 <- as.Date(paste(yyyy1,'01-01',sep='-'))
      t2 <- as.Date(paste(yyyy1,imon,sep='-'))
      imon <- t2 - t1
      if (verbose) print(imon)
      x <- lag(x,imon)
    }
  } else if (is.numeric(start)) x <- lag(x,start)
  x <- attrcp(x0,x)
  class(x) <- class(x0)
  return(x)
}
