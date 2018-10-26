
# Used to estimate annual statistics
## class creation
#annual <- function(x) structure(floor(x + .0001), class = "annual")

annual <- function(x, ...) UseMethod("annual")


  
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
  }  else
  if (inherits(x,'month')) {
    if (is.null(nmin)) nmin <- 12
  } else if (inherits(x,'season')) {
    if (is.null(nmin)) nmin <-  length(levels(factor(month(x))))
  } else nmin <- NA
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


annual.station <- function(x,FUN='mean',nmin=NULL,threshold=NULL,verbose=FALSE,...) {
  if (verbose) print('annual.station')
  attr(x,'names') <- NULL
#  if (inherits(x,'day')) {
#    y <- annual.station.day(x,FUN=match.fun(FUN),...)
#  } else {
#  ns <- length(x[1,])
#  for (i in 1:ns) {
#    y <- annual.default(x,FUN=match.fun(FUN),nmin=nmin,...)
    y <- annual.default(x,FUN=FUN,nmin=nmin,threshold=threshold,verbose=verbose,...) ## threshold=threshold,
#    if (i==1) y <- z else y <- c(y,z)
#  }
   y[which(is.infinite(y))] <- NA
  invisible(y)
}

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

# AM need to enhance this function
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

year <- function(x) {
  #str(x); print(class(x)); print(index(x))
  if (inherits(x,'integer')) x <- as.numeric(x)
  
  if ( (inherits(x,'numeric')) & (min(x,na.rm=TRUE) > 0) &
      (max(x,na.rm=TRUE) < 3000) )
    return(x)
  
  if (inherits(x,c('station','field','zoo'))) {
    y <- year(index(x))
    return(y)
  }
  if (inherits(x,'trajectory')) {
    y <- strptime(x[,colnames(x)=='start'],format="%Y%m%d%H")$year + 1900
    return(y)
  }
  if (inherits(x,'events')) {
    y <- strptime(x$date,format="%Y%m%d")$year + 1900
    return(y)
  }
  if (inherits(x,c("POSIXt","PCICt"))) {
    y <- as.numeric(format(x, '%Y'))
    return(y)
  }
  if ( (class(x)[1]=="character") & (nchar(x[1])==10)) {
    y <- year(as.Date(x))
    return(y)
  }
  if (class(index(x))=="numeric") {
    y <- trunc(index(x))
    return(y)
  }
  #print("here"); print(index(x))
  if (class(x)[1]=="Date")
    y <- as.numeric(format(x, '%Y')) else
  if (class(x)[1]=="yearmon") y <- trunc(as.numeric(x)) else
  if (class(x)[1]=="yearqtr") y <- trunc(as.numeric(x)) else
  if (class(x)[1]=="character") y <- trunc(as.numeric(x)) else
  if (class(x)[1]=="numeric") y <- trunc(x) else
  if (class(x)[1]=="season") {
    # If season, then the first month is really the December month of the previous year
     month <- round(12*(as.numeric(index(x)) - trunc(as.numeric(index(x)))) + 1)
     y[is.element(month,1)] <- y[is.element(month,1)] - 1
   } else print(paste(class(x)[1],' confused...'))
  return(y)
}


month <- function(x) {
  if ( (inherits(x,c('numeric','integer'))) & (min(x,na.rm=TRUE) > 0) & (max(x,na.rm=TRUE) < 13) )
    return(x) 
  if ( (inherits(x,c('numeric','integer'))) & (min(x,na.rm=TRUE) > 0) ) y <- rep(1,length(x))
  if (inherits(x,c('station','field','zoo'))) {
    #
    y <- month(index(x))
    return(y)
  }
  if (inherits(x,'trajectory')) {
    y <- strptime(x[,colnames(x)=='start'],format="%Y%m%d%H")$mon + 1
    return(y)
  }
  if (inherits(x,'events')) {
    y <- strptime(x$date,format="%Y%m%d")$mon + 1
    return(y)
  }
  if (inherits(x,c("POSIXt","PCICt"))) {
    y <- as.numeric(format(x, '%m'))
    return(y)
  }
  if ( (class(x)[1]=="character") & (nchar(x[1])==10)) {
    y <- month(as.Date(x))
    return(y)
  }
  #print(class(index(x)))
  if (class(x)[1]=="Date")
    y <- as.numeric(format(x, '%m')) else
  if (class(x)[1]=="yearmon")
    y <- round(12*(as.numeric(x) - trunc(as.numeric(x))) + 1) else
  if (class(x)[1]=="yearqtr") y <- round(12*(as.numeric(x) - trunc(as.numeric(x))) + 1)
  if (class(x)[1]=="season") {
    # If season, then the first month is really the months are DJF, MAM, JJA, and OND:
    y <- y - 1
    y[y==-1] <- 12
  }
  return(y)
}

day <- function(x) {
  if (inherits(x,c('station','field','zoo'))) {
    y <- day(index(x))
    return(y)
  }
  if (inherits(x,'trajectory')) {
    y <- strptime(x[,colnames(x)=='start'],format="%Y%m%d%H")$mday
    return(y)
  }
  if (inherits(x,'events')) {
    y <- strptime(x$date,format="%Y%m%d")$mday
    return(y)
  }
  if (inherits(x,c('numeric','integer'))) x <- as.numeric(x)
  if ( (inherits(x,c('numeric','integer'))) & (min(x,na.rm=TRUE) > 0) & (max(x,na.rm=TRUE) < 32) )
    return(x)
  if ( (inherits(x,c('numeric','integer'))) & (min(x,na.rm=TRUE) > 0) ) y <- rep(1,length(x))  
  if (inherits(x,c("POSIXt","PCICt"))) {
    y <- as.numeric(format(x, '%d'))
    return(y)
  }
  if ( (class(x)[1]=="character") & (nchar(x[1])==10) ) {
    y <- day(as.Date(x))
    return(y)
  }
  if ( (class(x)[1]=="character") & (nchar(x[1])==4) ) {
    y <- rep(1,length(x)) 
    return(y)
  }
  if (class(x)[1]=="Date") y <- as.numeric(format(x, '%d'))
  return(y)
}

# Used to estimate Dec-Feb, Mar-May, Jun-Aug, and Sep-Nov statistics
# Manipulate the zoo-object by shifting the year/chonology so that
# zoo thinks the year defined as December-November is January-December.

#season <- function(x,format="numeric", ...) {
#  if (inherits(x,'integer')) x <- as.numeric(x)
#  if ( (inherits(x,'numeric')) & (min(x,na.rm=TRUE) > 0) & (max(x,na.rm=TRUE) < 5) )
#    return(x)
#  
#  if (inherits(x,c('station','field','zoo'))) {
#    y <- season(index(x))
#    return(y)
#  }
#  if ( (class(x)[1]=="character") & (nchar(x[1])==10) ) {
#    y <- season(as.Date(x))
#    return(y)
# }  
#  if (class(x)[1]=="Date") {
#    xl <- x
#    n <- length(x)
#    xl[2:n] <- x[1:n-1]
#    xl[1] <- NA
#    y <- yearqtr(xl)    
#    y <- as.numeric(format(x, '%q'))
#  }
#  if (format=="character") y <- season.name[y+1]
#  return(y)
#}

season <- function(x, ...) UseMethod("season")

season.default <- function(x,format="character") {
  nt <- length(index(x))
  season <- rep('',nt)
  m <- month(x)
  if ( (inherits(x,'zoo')) & (format=="character") ) {
    for (i in 1:nt)  season[i] <- switch(m[i],
                                        '1'='djf','2'='djf','12'='djf',
                                         '3'='mam','4'='mam','5'='mam',
                                         '6'='jja','7'='jja','8'='jja',
                                         '9'='son','10'='son','11'='son')
  } else if ( (inherits(x,'zoo')) & (format=="numeric") ){
    for (i in 1:nt)  season[i] <- switch(m[i],'1'=1,'2'=1,'12'=1,
                                         '3'=2,'4'=2,'5'=2,
                                         '6'=3,'7'=3,'8'=3,
                                         '9'=4,'10'=4,'11'=4)
    season <- as.numeric(season)
  } else {
    season <- paste(substr(month.abb[as.numeric(rownames(table(month(x))))],1,1),sep='')
  }
#
  season
}

seasonal.yearmon <- function(x) {

  attr(season,'history') <- history.stamp(x)
  season
}

season.abb <- function() {
  season.abb <- c('annual','djf','jfm','fma','mam','amj',
                  'mjj','jja','jas','aso',
                  'son','ond','ndj','ondjfm','amjjas',
                  'ndjf','jjas','mjjas',
                  'djfm','djfma','ndjfma',
                  'ndjfmam','ondjfma')
  season<-list(1:12,c(12,1,2),1:3,2:4,3:5,4:6,5:7,6:8,7:9,8:10,9:11,10:12,
               c(11,12,1),c(10:12,1:3),4:9,c(11,12,1,2),6:9,5:9,
               c(12,1,2,3),c(12,1,2,3,4),c(11,12,1,2,3,4),c(11,12,1,2,3,4,5),
               c(10,11,12,1,2,3,4))
  season.abb <- c(season.abb,toupper(season.abb))
  season <- rep(season,2)
  names(season) <- season.abb
  season
}



pentad <- function(x,l=5,it0=NULL,...) {
  if (!is.null(it0)) yr <- year(x) - it0 else yr <- year(x)
  yrl <- l*trunc(yr/l)
  if (!is.null(it0)) yrl <- yrl + it0  
  index(x) <- yrl
  
  xl <- aggregate(x,yrl,...)
  attr(xl,'dimnames') <- NULL
  xl
}

