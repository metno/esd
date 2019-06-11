#' Conversion to esd objects.
#' 
#' \code{as.annual} and \code{annual} aggregates time series into annual values (e.g. means).
#'
#' \code{as.monthly} aggregates time series into monthly values (e.g. means).
#'
#' \code{as.4seasons} aggregates to four seasons ('djf': December-February, 'mam': March-May, 'jja': June-August, 'son': September-November)
#' and \code{as.seasons} aggregates to a user defined season (see arguments 'start' and 'end'). 
#'
#' @param x input object, e.g., a 'station' or 'field' object
#' @param FUN a function
#' @param verbose a boolean; if TRUE print information about progress
#' @param nmin minimum number of data points in a season
#' @param threshold threshold used if FUN is, e.g., 'exceedance' or 
#' @param slow a boolean; if FALSE run a fast version that might not work on all data types?
#' @param dateindex a boolean; if TRUE transform index into date format
#' @param start first month and day of the season, e.g., '01-01' for January 1 (argument in \code{season})
#' @param end last month and day of the season, e.g., '12-31' for December 31 (argument in \code{season})
#' @param \dots additional arguments

#' @aliases as.annual as.annual.default as.annual.numeric as.annual.integer as.annual.yearqtr as.annual.station as.annual.spell
#' annual annual.zoo annual.default annual.dsensemble annual.station annual.spell annual.field
#' as.monthly as.monthly.default as.monthly.station as.monthly.field
#' as.4seasons as.4seasons.default as.4seasons.day as.4seasons.station as.4seasons.spell as.4seasons.field as.4seasons.dsensemble as.seasons
#'
#' @seealso aggregate
#' 
#' @param x an input object of, e.g., class 'station', 'field'
#' @param FUN a function, see \code{\link{aggregate.zoo}}
#' @param nmin Minimum number of data points (e.g. days or months) with valid
#' data accepted for annual estimate. NULL demands complete years.
#' @param format 'numeric' or 'character'
#' @param na.rm a boolean; if TRUE, ignore NA - see see \code{\link{mean}}
#' @param verbose a boolean; if TRUE print information about progress
#' @param \dots additional argument
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
#' plot(as.4seasons(bjornholt,threshold=1,FUN="exceedance"))
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

#' @export
as.annual <- function(x, ...) UseMethod("as.annual")

#' @export
as.annual.default <- function(x, ...) annual(x, ...)

#' @export
as.annual.numeric <- function(x, ...) annual(x, ...)

#' @export
as.annual.integer <- function(x, ...) structure(x, class = "annual")

#' @export
as.annual.yearqtr <- function(x, frac = 0, ...) {
    if (frac == 0) annual(as.numeric(x)) else
    as.annual(as.Date(x, frac = frac), ...)
}

#' @export
as.annual.station <- function(x, ...) annual.station(x,...)

#' @export
as.annual.spell <- function(x, ...) annual.spell(x,...)

#' @export
as.monthly <- function(x,...) UseMethod("as.monthly")

yyyymm <- function(x) ym <- as.Date(paste(year(x),month(x),'01',sep='-'))

#' @export
as.monthly.default <- function(x,...) {
  y <- aggregate(x,by=yyyymm,...)
  return(y)
}

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
as.4seasons <- function(x,...) UseMethod("as.4seasons")

#' @export
as.4seasons.default <- function(x,...,FUN='mean',slow=FALSE,verbose=FALSE,nmin=NULL) {
  if(verbose) print('as.4seasons.default')
  if (inherits(x,'season')) return(x)
  attr(x,'names') <- NULL
  d <- dim(coredata(x))
  #print(d)
  if (is.null(d)) d <- c(length(x),1)
  if (!slow) {
    if (inherits(x,"month")) {
      if ( (is.null(d)) | (d[2]==1) ) {
        X <- c(NA,coredata(x)[1:length(x)-1]) # shift the coredata by 1 to start on December. This works only for monthly data !!!  
      } else {
        X <- rbind(rep(NA,d[2],1),coredata(x)[1:d[1]-1,])
      }
      #print(dim(X))
      X <- zoo(X,order.by=index(x))
      ##yrseas <- fourseasons(ix)
      ##print(yrseas)
      #print('aggregate')
      #print(names(list(...)))
      yq <- function(t) as.yearqtr(year(t) + 0.25*floor((month(t)-1)/3))
      y <- aggregate(x=as.zoo(X),by=yq,#as.yearqtr,
                     FUN=match.fun(FUN),...)
      # convert yearqtr to yearmon
      y <- zoo(x=y,order.by=as.Date(as.yearmon(index(y))))
    } else y <- as.4seasons.day(x,FUN=FUN,nmin=nmin,...)
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

    #print("start loop")
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
    #print("end loop")
    ok <- is.finite(rowMeans(X,na.rm=TRUE)) & is.finite(t)
    #print(table(as.yearqtr(t[ok])))
    #print(summary(c(X)))
    y <- zoo(X[ok,],order.by=as.Date(as.yearqtr(t[ok])))
    #names(y) <- FUN
  }
  if (!is.null(dim(y))) {
    ok <- is.finite(rowMeans(y,na.rm=TRUE))
    y <- y[ok,]
  } 
  y <- attrcp(x,y)
  attr(y,'history') <- history.stamp(x)
  if (inherits(x,'field'))
    attr(y,'dimensions') <- c(attr(x,'dimensions')[1:2],sum(ok))
  class(y) <- class(x)
  class(y)[2] <- "season"
  #plot(y)
  return(y) 
}

#' @export as.4seasons.day
as.4seasons.day <- function(x,...,FUN='mean',na.rm=TRUE,dateindex=TRUE,nmin=85,verbose=FALSE) {
  if(verbose) print('as.4seasons.day')
  IV <- function(x) sum(is.finite(x))
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
  nd <- aggregate(X,as.yearqtr,FUN=IV)
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
  class(y) <- class(x)
  class(y)[2] <- "season"
  invisible(y)
}

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

#' @export
as.4seasons.spell <- function(x,...,FUN='mean') {
  y <- as.4seasons.default(as.station(x),FUN=FUN,...)
#  y <- attrcp(x,y)
#  attr(y,'history') <- history.stamp(x)
#  class(y) <- class(x)
#  class(y) <- gsub("month","season",class(x))
  return(y) 
}

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

# Not to confuse with season
# This function extracts a given seasonal interval and aggregates a given statistic
#' @export
as.seasons <- function(x,start='01-01',end='12-31',FUN='mean',verbose=FALSE,...) {
  if(verbose) print("as.seasons")
  IV <- function(x) sum(is.finite(x))
  yrs <- year(x); d <- dim(x)
  # ns = number of stations
  if (is.null(d)) ns <- 1 else ns <- d[2]
  years <- as.numeric(rownames(table(yrs))); n <- length(years)
  y <- matrix(rep(NA,n*ns),n,ns); k <- y
  start.1 <- as.numeric(as.Date(paste(years[1],start,sep='-')))
  end.1 <- as.numeric(as.Date(paste(years[1],end,sep='-')))
  if (start.1 > end.1) twoyears <- 1 else twoyears <- 0
  
  for (i in 1:n) {
    z <- coredata(window(x,start=as.Date(paste(years[i],start,sep='-')),
                           end=as.Date(paste(years[i]+twoyears,end,sep='-'))))
    k[i,] <- apply(matrix(z,length(z),ns),2,IV)
    y[i,] <- apply(matrix(z,length(z),ns),2,FUN, ...)
  }
  y <- zoo(y,order.by=as.Date(paste(years,start,sep='-')))
  y <- attrcp(x,y)
  attr(y,'history') <- history.stamp(x)
  if (twoyears==0) attr(y,'season.interval') <- paste(start,'to',end) else
                   attr(y,'season.interval') <- paste(start,'to',end,'the following year')
  attr(y,'n.valid') <- k
  class(y) <- class(x)
  class(y)[2] <- "annual"
  return(y)
}
