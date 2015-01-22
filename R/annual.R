
# Used to estimate annual statistics
## class creation
#annual <- function(x) structure(floor(x + .0001), class = "annual")

annual <- function(x, ...) UseMethod("annual")


  
annual.zoo <- function(x,FUN='mean',na.rm=TRUE,nmin=NULL, ...) {
  #print("annual.zoo")
  if (inherits(x,'annual')) return(x)
  attr(x,'names') <- NULL
#  yr <- year(x)  REB: 08.09.2014
  class(x) <- 'zoo'
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
  invisible(y)
}

#annual.station.day <- function(x,FUN=mean,na.rm=TRUE,nmin=350, ...) {
#  attr(x,'names') <- NULL
  #print(class(index(x)))
#  yr <- year(x)
#  years <- as.numeric(rownames(table(yr)))
  #print(table(years))
#  ndyr <- as.numeric(table(yr))
  #print(table(ndyr))
#  ok <- is.element(yr,years[ndyr > nmin]) 
  #print(sum(ok))
#  X <- zoo(x[ok,],order.by=index(x)[ok])
#  yr <- yr[ok]
#  attr(X,'unit') <- attr(x,'unit')
  #print(ndyr)
#  cls <- class(x)
#  class(x) <- "zoo"
#  y <- aggregate.station(X,yr,match.fun(FUN),...,na.rm=na.rm)
#  unit <- attr(y,'unit')
#  y <- attrcp(x,y,ignore="unit")
#  unit -> attr(y,'unit')
#  #str(y); print(unit)
#  attr(y,'history') <- history.stamp(x)
#  class(y) <- cls
#  class(y)[2] <- "annual"
#  invisible(y)
#}



annual.default <- function(x,FUN='mean',na.rm=TRUE, nmin=NULL,...,
                           regular=NULL,frequency=NULL) {

  #print('annual.default')
  if (inherits(x,'annual')) return(x)
  nv <- function(x) sum(is.finite(x))

  #browser()
  if (inherits(FUN,'function')) FUN <- deparse(substitute(FUN)) # REB110314
  attr(x,'names') <- NULL
  yr <- year(x)
  nmo <- length(levels(factor(month(x))))
  d <- dim(x)
  if(is.null(d)) d <- c(length(x),1)
  #print(table(yr))
  YR <- as.numeric(rownames(table(yr)))
  nyr <- as.numeric(table(yr))
  nval <- aggregate(zoo(x,order.by=index(x)),year,FUN=nv)
  #print(c(length(nval),length(nyr),length(yr),length(YR))); plot(nval)
  #print(YR)
  #print(class(x))
  # Need to accomodate for the possibility of more than one station series.
  if (inherits(x,'day')) {
    if (is.null(nmin)) nmin <- 30*nmo
    fewd <- coredata(nval) < nmin
    nyr[fewd] <- coredata(nval)[fewd]
    #print(fewd)
    #x[fewd] <- NA
    ok <- is.element(yr,YR[nyr >= nmin])
    #browser()
  } else   if ( (inherits(x,'mon')) & is.null(nmin) ) {
      iy <- year(x)
      nmy <- as.numeric(table(iy))
      full <- nmy[nmy==12]
      x[is.element(iy,!full),] <- NA
      na.rm=FALSE
      N <- 12 
  } else
  if (inherits(x,'month')) {
    OK <- nyr == 12
    ok <- is.element(yr,YR[OK])
  } else if (inherits(x,'season')) {
    OK <- nyr == 4
    ok <- is.element(yr,YR[OK])
  } else ok <- is.finite(yr)
  #print(c(sum(ok),length(ok),nmin)); print(YR[is.element(YR,yr[ok])])

  # Make a new zoo-object without incomplete years
#  if (length(d)==2) X <- zoo(coredata(x[ok,]),order.by=index(x)[ok]) else
#                    X <- zoo(coredata(x[ok]),order.by=index(x)[ok])
  # REB 2015-01-16: the two commented-out lines produced errors in some cases; lines below are more robust.
  X <- zoo(coredata(x),order.by=index(x))
  if (sum(ok)>0) coredata(X)[!ok] <- NA 
  #print(summary(X))
  if (FUN == 'sum') na.rm <- FALSE ## AM
  #y <- aggregate(X,yr[ok],FUN=FUN,...,na.rm=na.rm) ## AM
  #browser()
  #y <- aggregate(X,year,FUN==FUN,...,na.rm=na.rm) ## AM
  #print(names(list(...))); print(names(formals(FUN)))
  #browser()
  #print(FUN); print(sum(is.element(names(formals(FUN)),'na.rm')))
  if ( (sum(is.element(names(formals(FUN)),'na.rm')==1)) |
       (sum(is.element(FUN,c('mean','min','max','sum','quantile')))>0) )
    y <- aggregate(X,year,FUN=FUN,...,na.rm=na.rm) else
    y <- aggregate(X,year,FUN=FUN,...) # REB.
  y[!is.finite(y)] <- NA ## AM
  y <- attrcp(x,y,ignore="names")

  args <- list(...)
  #print(names(args))
  ix0 <- grep('threshold',names(args))
  if (length(ix0)>0) threshold <- args[[ix0]] else threshold <- 1
  if (FUN=="counts")  {
    #print("Count")
    attr(y,'unit') <-
      rep(paste("counts | X >",threshold," * ",attr(x,'unit')),d[2])
  } else if (FUN=="freq") {
    #print("Frequency")
    attr(y,'variable') <- rep('f',d[2])
    attr(y,'unit') <- 'fraction'
#    attr(y,'unit') <- rep(paste("frequency | X >",threshold," * ",attr(x,'unit')),d[2])
  } else if (FUN=="wetfreq") {
    #print("Wet-day frequency")
    attr(y,'variable') <- rep('f[w]',d[2])
    attr(y,'unit') <- 'fraction'
#    attr(y,'unit') <- rep(paste("frequency | X >",threshold," * ",attr(x,'unit')),d[2])
  } else if (FUN=="wetmean") {
    #print("Wet-day mean")
    attr(y,'variable') <- rep('mu',d[2])
    attr(y,'unit') <- rep('mm/day',d[2])
#    n <- count(X,threshold=threshold) # REB
    n <- aggregate(X,year,FUN='count',threshold=threshold, ...,
                   regular = regular, frequency = frequency)
    bad <- coredata(n)==0
    coredata(n)[bad] <- 1
    std.err <- 2*coredata(y)/sqrt(coredata(n)-1)
    std.err[bad] <- NA
    attributes(std.err) <- NULL
    dim(std.err) <- dim(y)
    attr(y,'standard.error') <- zoo(std.err,order.by=index(y))
  } else if (FUN=="mean") {
    #print("mean")
    sigma <- aggregate(X, year, FUN='sd', ...,
                       regular = regular, frequency = frequency)
#    n <- count(x,threshold=threshold)
    n <- aggregate(X,year,FUN='count',threshold=threshold, ...,
                   regular = regular, frequency = frequency)
    #browser()
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
  return(y)
}


annual.station <- function(x,FUN='mean',nmin=NULL,...) {
  #print('annual.station')
  attr(x,'names') <- NULL
#  if (inherits(x,'day')) {
#    y <- annual.station.day(x,FUN=match.fun(FUN),...)
#  } else {
#  ns <- length(x[1,])
#  for (i in 1:ns) {
#    y <- annual.default(x,FUN=match.fun(FUN),nmin=nmin,...)
    y <- annual.default(x,FUN=FUN,nmin=nmin,...)
#    if (i==1) y <- z else y <- c(y,z)
#  }
   y[which(is.infinite(y))] <- NA
  return(y)
}

annual.spell <- function(x,FUN='mean',nmin=0,...) {
  attr(x,'names') <- NULL
  if ( (inherits(x,'mon'))  & is.null(nmin) ) {
    iy <- year(x)
    nmy <- as.numeric(table(iy))
    full <- nmy[nmy==12]
    x[is.element(iy,!full),] <- NA
    na.rm=FALSE
  }
  #y <- annual.default(x,FUN=match.fun(FUN),...)
  y <- annual.default(x,FUN=FUN,nmin=nmin,...)
  return(y)
}

annual.dsensemble <- function(x,FUN='mean') {
  #print("annual.dsensemble")
  y <- subset(x,it=0)
  return(y)
}



annual.field <- function(x,FUN='mean',na.rm=TRUE,nmin=NULL, ...) {
  attr(x,'names') <- NULL
  if (inherits(FUN,'function')) FUN <- deparse(substitute(FUN)) # REB110314
  yr <- year(x)
  cls <- class(x)
  class(x) <- "zoo"
  if ( (inherits(x,'mon'))  & is.null(nmin) ) {
    iy <- year(x)
    nmy <- as.numeric(table(iy))
    full <- nmy[nmy==12]
    x[is.element(iy,!full),] <- NA
    na.rm=FALSE
  }
#  y <- aggregate(x,yr,FUN=match.fun(FUN),...,na.rm=na.rm)
  y <- aggregate(x,yr,FUN=FUN,...,na.rm=na.rm)
  y <- attrcp(x,y)
  attr(y,'history') <- history.stamp(x)
  attr(y,'dimensions') <- c(attr(x,'dimensions')[1:2],length(index(y)))
  class(y) <- cls
  class(y)[2] <- "annual"
  y[which(is.infinite(y))] <- NA
  invisible(y)
}

year <- function(x) {
  #str(x); print(class(x)); print(index(x))
  if (inherits(x,'integer')) x <- as.numeric(x)
  
  if ( (inherits(x,'numeric')) & (min(x,na.rm=TRUE) > 0) & (max(x,na.rm=TRUE) < 3000) )
    return(x)
  
  if (inherits(x,c('station','field','zoo'))) {
    y <- year(index(x))
    return(y)
  }
  if ( (class(x)[1]=="character") & (nchar(x[1])==10) ) {
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
    #browser()
    y <- month(index(x))
    return(y)
  }
  if ( (class(x)[1]=="character") & (nchar(x[1])==10) ) {
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
  if (inherits(x,c('numeric','integer'))) x <- as.numeric(x)
  if ( (inherits(x,c('numeric','integer'))) & (min(x,na.rm=TRUE) > 0) & (max(x,na.rm=TRUE) < 32) )
    return(x)
  if ( (inherits(x,c('numeric','integer'))) & (min(x,na.rm=TRUE) > 0) ) y <- rep(1,length(x))  
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
#  browser()
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
                  'dfjm','dfjma','ndjma')
  season<-list(1:12,c(12,1,2),1:3,2:4,3:5,4:6,5:7,6:8,7:9,8:10,9:11,10:12,
                c(11,12,1),c(10:12,1:3),4:9,c(11,12,1,2),6:9,5:9,
                c(12,1,2,3),c(12,1,2,3,4),c(11,12,1,2,3,4))
  names(season) <- season.abb
  season
}



pentad <- function(x,l=5,...) {
  yrl <- l*trunc(year(x)/l)
  index(x) <- yrl
  
  xl <- aggregate(x,yrl,...)
  xl
}

