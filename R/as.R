
# KMP 2019-05-14: seems obsolete and does not return a ds object
#as.ds <- function(x,verbose=FALSE) {
#  if(verbose) print(as.ds)
#  y <- zoo(x,order.by=index(x))
#  attr(y,'location') <- attr(x,'location')
#  attr(y,'variable') <- attr(x,'variable')
#  attr(y,'unit') <- attr(x,'unit')
#  attr(y,'longitude') <- attr(x,'longitude')
#  attr(y,'latitude') <- attr(x,'latitude')
#  attr(y,'altitude') <- attr(x,'altitude')
#  attr(y,'country') <- attr(x,'country')
#  attr(y,'longname') <- attr(x,'longname')
#  attr(y,'station_id') <- attr(x,'station_id')
#  attr(y,'quality') <- attr(x,'quality') 
#  attr(y,'calendar') <- attr(x,'calendar')
#  attr(y,'source') <- attr(x,'source')
#  attr(y,'URL') <- attr(x,'URL')
#  #attr(y,'history') <- attr(x,'history')
#  #attr(y,'date-stamp') <- attr(x,'date-stamp')
#  attr(y,'type') <- attr(x,'type')
#  attr(y,'aspect') <- attr(x,'aspect')
#  attr(y,'reference') <- attr(x,'reference')
#  attr(y,'info') <- attr(x,'info')
#  attr(y,'history') <- history.stamp(x)
#  class(y) <- c("station","month","zoo")
#  return(y)
#}



as.annual <- function(x, ...) UseMethod("as.annual")

as.annual.default <- function(x, ...) annual(x, ...)

as.annual.numeric <- function(x, ...) annual(x, ...)

as.annual.integer <- function(x, ...) structure(x, class = "annual")

as.annual.yearqtr <- function(x, frac = 0, ...) {
    if (frac == 0) annual(as.numeric(x)) else
    as.annual(as.Date(x, frac = frac), ...)
}

as.annual.station <- function(x, ...) annual.station(x,...)

as.annual.spell <- function(x, ...) annual.spell(x,...)

as.monthly <- function(x,...) UseMethod("as.monthly")

yyyymm <- function(x) ym <- as.Date(paste(year(x),month(x),'01',sep='-'))

as.monthly.default <- function(x,...) {
  y <- aggregate(x,by=yyyymm,...)
  return(y)
}

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

# Not to confuse with season
# This function extracts a given seasonal interval and aggregates a given statistic
as.seasons <- function(x,start='01-01',end='12-31',FUN='mean',...) {
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

as.4seasons <- function(x,...) UseMethod("as.4seasons")



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

as.4seasons.day <- function(x,...,FUN='mean',na.rm=TRUE,dateindex=TRUE,nmin=85) {
  IV <- function(x) sum(is.finite(x))
    if (inherits(x,'month')) nmin <- 3 # AM 06-07-2015
  #print('as.4seasons.day')
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

as.4seasons.station <- function(x,...,FUN='mean') {
  #print('as.4seasons.station')
  y <- as.4seasons.default(x,FUN=FUN,...)
#  y <- attrcp(x,y)
#  attr(y,'history') <- history.stamp(x)
#  class(y) <- class(x)
#  class(y) <- gsub("month","season",class(x))
  return(y) 
}

as.4seasons.spell <- function(x,...,FUN='mean') {
  y <- as.4seasons.default(as.station(x),FUN=FUN,...)
#  y <- attrcp(x,y)
#  attr(y,'history') <- history.stamp(x)
#  class(y) <- class(x)
#  class(y) <- gsub("month","season",class(x))
  return(y) 
}


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

as.anomaly <- function(x,...) UseMethod("as.anomaly")

# REB 2015-03-23 Tidy up - use anomaly

as.anomaly.default <- function(x,...,ref=NULL,na.rm=TRUE) anomaly.default(x,...)

#
#as.anomaly.default <- function(x,ref=NULL,na.rm=TRUE) {
## The argument monthly can be used to force the method to be
## julian-day regression-based or based on monthly mean
 
##  print('as.anomaly.default')
##  yr <- as.integer(format(index(x),'%Y'))
##  mon <- as.integer(format(index(x),'%m'))
##  dy <- as.integer(format(index(x),'%d'))
##  str(x)
##  browser()
#  if ((is.numeric(x)) & (!is.null(attr(x, "names"))))
#  if ((!is.zoo(x)) & (!is.null(attr(x, "names"))))
#    x <- zoo(x,order.by=as.Date(attr(x, "names")))
#  yr <- year(x);  mon <- month(x);  dy <- day(x); seas <- season(x,format='numeric')
#  if (is.null(ref))
#    ref <- seq(min(yr,na.rm=TRUE),max(yr,na.rm=TRUE),by=1)
#  # Check whether the object is a time series of monthly data.
#  ndd <- length(table(dy)); nmm <- max(diff(mon))
#  #print(c(ndd,nmm))
#  if ((ndd==0) & (nmm==0)) {
#    # Only one month/season
#    return(x - mean(x,na.rm=TRUE))
#  }
#  
#  #print(ndd); print(class(x))
#  #if ( (ndd==1) & is.null(monthly) ) monthly <- TRUE else
#  #                                   monthly <- FALSE
#
#  if ( (inherits(x,'month')) | ((ndd==1) & (nmm==1)) ) {
#    # Monthly data
#    print('monthly')
#    monthly <- TRUE
#    nm <- 12
#  } else monthly <- FALSE
#  if ( (inherits(x,'seasonal')) | ((ndd==1) & (nmm==3)) ) {
#    # seasonal data
#    #print('seasonal')
#    seasonal <- TRUE
#    mon <- seas
#    nm <- 4
#  } else seasonal <- FALSE
#  
#  # Check if x contains one or more stations
#  d <- dim(x)
#
#  if (is.null(d)) {
#    # single records
#    #cat('.')
#    y <- coredata(x)
#    
#    if ( (monthly | seasonal) ) {
#    # If monthly, subtract the monthly mean
##    print(mon)
##    if (is.null(dim(x))) {
#      print("1D")
#      clim <- rep(0,nm)
#      for (i in 1:nm) {
#        im <- is.element(mon,i)
#        z <- mean(y[im & is.element(yr,ref)],na.rm=na.rm)
#        clim[i] <- z
#        y[im] <- y[im] - clim[i]
#      }
#      
#      y <- zoo(y,order.by=index(x))
##    } else {
##      #print("2D")
##      y <- t(y)
##      clim <- rep(0,length(y[,1])*12); dim(clim) <- c(length(y[,1]),12)
##      for (i in 1:12) {
##        im <- is.element(mon,i)
##        z <- rowMeans(y[,im & is.element(yr,ref)],na.rm=na.rm) 
##        #print(c(length(z),length(clim[,1]))); plot(z)
##        clim[,i] <- z
##        y[,im] <- y[,im] - clim[,i]
##        #print(c(i,sum(im),mean(clim[,i])))
##      }
##      y <- zoo(t(y),order.by=index(x))
#
#    } else {
##      print("daily"); print(class(x))
#      t0 <- julian(index(x)[is.element(yr,ref)]) -
#            julian(as.Date(paste(yr[is.element(yr,ref)],"-01-01",sep="")))
#      t <- julian(index(x)) -
#         julian(as.Date(paste(yr,"-01-01",sep="")))
#      c1 <- cos(pi*t0/365.25); s1 <- sin(pi*t0/365.25)
#      c2 <- cos(2*pi*t0/365.25); s2 <- sin(2*pi*t0/365.25)
#      c3 <- cos(3*pi*t0/365.25); s3 <- sin(3*pi*t0/365.25)
#      c4 <- cos(4*pi*t0/365.25); s4 <- sin(4*pi*t0/365.25)
#      C1 <- cos(pi*t/365.25);   S1 <- sin(pi*t/365.25)
#      C2 <- cos(2*pi*t/365.25); S2 <- sin(2*pi*t/365.25)
#      C3 <- cos(3*pi*t/365.25); S3 <- sin(3*pi*t/365.25)
#      C4 <- cos(4*pi*t/365.25); S4 <- sin(4*pi*t/365.25)
#      cal <- data.frame(y=coredata(x),c1=c1,c2=c2,c3=c3,c4=c4,
#                        s1=s1,s2=s2,s3=s3,s4=s4)
#      pre <- data.frame(c1=C1,c2=C2,c3=C3,c4=C4,
#                        s1=S1,s2=S2,s3=S3,s4=S4)
#      i1 <- is.element(year(x),year(x)[1])
#      pre1 <- data.frame(c1=C1[i1],c2=C2[i1],c3=C3[i1],c4=C4[i1],
#                         s1=S1[i1],s2=S2[i1],s3=S3[i1],s4=S4[i1])
#      acfit <- lm(y ~ c1 + s1 + c2 + s2 + c3 + s3 + c4 + s4,data=cal)
#      clim <- predict(acfit,newdata=pre)
#      y <- zoo(coredata(x) - clim,order.by=index(x))
#      clim <-  predict(acfit,newdata=pre1)
#    }
#  } else {
#    #print("many stations")
#    rownames(x) <- as.character(index(x))
#    y <- apply(x,2,FUN='as.anomaly.default',ref=ref)
#    y <- zoo(y,order.by=index(x))
#    clim <- x - y
#  }
#  #print("attributes")
#  y <- attrcp(x,y)
#  attr(y,"aspect") <- 'anomaly'
#  attr(y,"anomaly_method") <- monthly
#  attr(y,"climatology") <- clim
#  #attr(y,"date") <- date()
#  #attr(y,"call") <- match.call()
#  attr(y,'history') <- history.stamp(x)
#  class(y) <- class(x)
#  invisible(y)
#}

as.anomaly.zoo <- function(x,ref=NULL,na.rm=TRUE,...) {
  y <- as.anomaly.station(x,ref=ref,na.rm=na.rm,...)
  attr(y,'history') <- history.stamp(x)
  invisible(y)
}

as.anomaly.list <- function(x,ref=NULL,na.rm=TRUE,...) {
  y <- lapply(x,anomaly(x))
  attr(y,'history') <- history.stamp(x)
  invisible(y)
}

as.anomaly.station <- function(x,ref=NULL,na.rm=TRUE,...) {
  y <- as.anomaly.default(x,ref=ref,na.rm=na.rm,...)
  attr(y,'history') <- history.stamp(x)
  invisible(y)
}

as.anomaly.field<- function(x,ref=NULL,na.rm=TRUE,...) {
   y <- anomaly.default(x,ref=ref,na.rm=na.rm,...)
   attr(y,'history') <- history.stamp(x)
   attr(y,'dimensions') <- attr(x,'dimensions')
   invisible(y)
}

# Handy conversion algorithms:
as.climatology <- function(x,...) {
    ya <- as.anomaly(x,...)
    clim <- coredata(attr(ya,'climatology'))
    if (!is.null(dim(clim)))
        len.clim <- dim(clim)[1]
    else
        len.clim <- length(clim)
    y <- zoo(clim,order.by=1:len.clim)      
    y <- attrcp(x,y)
    attr(y,'aspect') <- 'climatology'
    attr(y,'history') <- history.stamp(x)
    class(y) <- class(x)
    invisible(y)
}

as.residual <- function(x,...) UseMethod("as.residual")

as.residual.ds <- function(x,...,verbose=FALSE){
  if (verbose) print('as.residual.ds')
  if (is.ds(x)) {
    ## If the predictand was originally an EOF or PCA product, then
    ## the residual needs to inherits their attributes
      if (verbose) print('Re-construct and re-compute')
      ## Need to reconstruct the data matrix and re-calculate the EOFs/PCAs 
      if (is.eof(x)) {
        if (verbose) print('eof/field')
        z0 <- as.field(attr(x,'original_data'))
        z1 <- as.field(x)
        y <- z1 - z0
        y <- attrcp(z0,y); class(y) <- class(z0)
      } else
      if (is.pca(x)) {
        if (verbose) print('pca/station')
        z0 <- as.station(attr(x,'original_data'))
        z1 <- as.station(x)
        y <- z1 - z0
        y <- attrcp(z0,y); class(y) <- class(z0)
      } else
      if (is.station(x)) {
        if (verbose) print('station')
        z0 <- attr(x,'original_data')
        z1 <- x
        y <- z1 - z0
        y <- attrcp(z0,y); class(y) <- class(z0)
      }      
   } else
  ## If the results are a field object, then the residuals are stored as EOFs.
  if (is.field(x)) {
    if (verbose) {print('x is a field object'); print(class(x))}
    y <- as.field(attr(x,'original_data')) -
         as.field(attr(x,'fitted_values'))
    y <- attrcp(attr(x,'original_data'),y)
    class(y) <- class(attr(x,'original_data'))
    y <- as.field(y)
    attr(y,'aspect') <- 'residual'
  }
  attr(y,'aspect') <- 'residual'
  attr(y,'history') <- history.stamp(x)
  invisible(y)
}

as.residual.station <- function(x,...) {
  if (!is.null(attr(x,'calibration_data')))
    y <- as.residual.ds(x) else y <- NULL
  invisible(y)
}

as.residual.field <- function(x,...) {
  if (!is.null(attr(x,'calibration_data')))
    y <- as.residual.ds(x) else y <- NULL
  invisible(y)
}


as.calibrationdata <- function(x) UseMethod("as.calibrationdata")

as.calibrationdata.ds <- function(x) {
  y <- attr(x,'calibration_data')
  y <- attrcp(x,y)
  attr(y,'history') <- history.stamp(x)
  class(y) <- class(attr(x,'calibration_data'))
  invisible(y)
}

as.calibrationdata.station <- function(x) {
   if (!is.null(attr(x,'calibration_data')))
    y <- as.calibrationdata.ds(x) else y <- NULL
  invisible(y)
}

as.fitted.values <- function(x) UseMethod("as.fitted.values")

as.fitted.values.ds <- function(x) {
  y <- attr(x,'fitted_values')
  y <- attrcp(x,y)
  attr(y,'history') <- history.stamp(x)
  class(y) <- class(attr(x,'fitted_values'))
  invisible(y)
}

as.fitted.values.station <- function(x) {
   if (!is.null(attr(x,'fitted.values')))
    y <- as.fitted.values.ds(x) else y <- NULL
  invisible(y)
}



as.original.data <- function(x) UseMethod("as.original.data")

as.original.data.ds <- function(x) {
  y <- attr(x,'original_data')
  y <- attrcp(x,y)
  attr(y,'history') <- history.stamp(x)
  class(y) <- class(attr(x,'original_data'))
  invisible(y)
}

as.original.data.station <- function(x) {
  y <- as.original.data.ds(x)
  invisible(y)
}


as.pattern <- function(x,...) UseMethod("as.pattern")

as.pattern.ds <- function(x,...,verbose=FALSE) {
  if(verbose) print("as.pattern.ds")
  y <- attr(x,'pattern')
  attr(y,'longitude') <- lon(x)
  attr(y,'latitude') <- lat(x)
  attr(y,'variable') <- paste('weights(',varid(x),')',sep='')
  attr(y,'unit') <- 'dimensionless'
  attr(y,'dimensions') <- dim(y)
  attr(y,'history') <- history.stamp(x)
  class(y) <- c('pattern',class(y))
  invisible(y)
}

as.pattern.eof <- function(x,...,verbose=FALSE) {
  if(verbose) print("as.pattern.eof")
  y <- attr(x,'pattern')
  attr(y,'longitude') <- lon(x)
  attr(y,'latitude') <- lat(x)
  attr(y,'variable') <- paste('weights(',varid(x),')',sep='')
  attr(y,'unit') <- 'dimensionless'
  attr(y,'dimensions') <- dim(y)
  attr(y,'history') <- history.stamp(x)
  class(y) <- c('pattern',class(y))
  invisible(y)  
}

as.pattern.mvr <- function(x,...,verbose=FALSE) {
  if(verbose) print("as.pattern.mvr")
  y <- attr(x,'pattern')
  attr(y,'longitude') <- lon(x)
  attr(y,'latitude') <- lat(x)
  attr(y,'variable') <- paste('weights(',varid(x),')',sep='')
  attr(y,'unit') <- 'dimensionless'
  attr(y,'dimensions') <- dim(y)
  attr(y,'history') <- history.stamp(x)
  class(y) <- c('pattern',class(y))
  invisible(y)  
}

as.pattern.cca <- function(x,...,verbose=FALSE) {
  if(verbose) print("as.pattern.cca")
  y <- attr(x,'pattern')
  attr(y,'longitude') <- lon(x)
  attr(y,'latitude') <- lat(x)
  attr(y,'variable') <- paste('weights(',varid(x),')',sep='')
  attr(y,'unit') <- 'dimensionless'
  attr(y,'dimensions') <- dim(y)
  attr(y,'history') <- history.stamp(x)
  class(y) <- c('pattern',class(y))
  invisible(y)  
}

as.pattern.trend <- function(x,...,verbose=FALSE) {
  if(verbose) print("as.pattern.trend")
  y <- attr(x,'pattern')
  attr(y,'longitude') <- lon(x)
  attr(y,'latitude') <- lat(x)
  attr(y,'variable') <- paste('weights(',varid(x),')',sep='')
  attr(y,'unit') <- 'dimensionless'
  attr(y,'dimensions') <- dim(y)
  attr(y,'history') <- history.stamp(x)
  class(y) <- c('pattern',class(y))
  invisible(y)  
}

as.pattern.matrix <- function(x,...,verbose=FALSE) {
  if(verbose) print("as.pattern.matrix")
  if(verbose) print("Unfinished function - returning input object")
  return(x)
}

as.pattern.array <- function(x,...,verbose=FALSE) {
  if(verbose) print("as.pattern.array")
  if(verbose) print("Unfinished function - returning input object")
  return(x)
}

as.pattern.field <- function(x,...,FUN=NULL,verbose=FALSE) {
  if(verbose) print("as.pattern.field")
  if (!is.null(FUN)) {
    y <- apply(x,2,FUN,...)
    dim(y) <- attr(x,'dimension')[1:2]
  } else {
    y <- t(coredata(x))
    dim(y) <- attr(x,'dimension')
  }
  attr(y,'longitude') <- lon(x)
  attr(y,'latitude') <- lat(x)
  attr(y,'variable') <- varid(x)
  attr(y,'unit') <- unit(x)
  class(y) <- c('pattern',class(y))
  attr(y,'time') <- index(x) 
  attr(y,'history') <- history.stamp(x)
  invisible(y)
}

as.pattern.corfield <- function(x,...,verbose=FALSE) {
  if(verbose) print("as.pattern.corfield")
  y <- coredata(x)
  dim(y) <- attr(x,'dimension')[1:2]
  y <- attrcp(x,y)
  attr(y,'longitude') <- lon(x)
  attr(y,'latitude') <- lat(x)
  attr(y,'variable') <- varid(x)
  attr(y,'unit') <- unit(x)
  class(y) <- c('pattern',class(y))
  attr(y,'time') <- index(x)
  invisible(y)
}

as.eof <- function(x,...) UseMethod("as.eof")

as.eof.zoo <- function(x,...) {
  class(x) <- c('eof','zoo')
  return(x)
}

as.eof.ds <- function(x,...,iapp=NULL) {
  y <- as.eof(attr(x,'eof'),iapp=iapp) 
  return(y)
}

as.eof.eof <-function(x,...,iapp=NULL) {
  if (inherits(x,'comb')) {
    x <- as.eof.comb(x,iapp=iapp) 
  } 
  return(x)
}
  
as.eof.comb <- function(x,...,iapp=NULL) {
  #print("as.eof.comb")
  stopifnot(inherits(x,'comb'))

  # if x is a 'field'
  if (!inherits(x,'eof')) x <- EOF(x)

  # assume x from now on is an 'eof'
  if (!is.null(iapp)) {
    y <- as.eof.appendix(x,iapp=iapp)
    return(y)
  }
  class(x) <- class(x)[-grep('comb',class(x))]
  napps <- attr(x,'n.apps')
  for (i in seq(napps)) {
    eval(parse(text=paste("attr(x,'appendix.",i,"') <- NULL",sep="")))
  }
  attr(x,'n.apps') <- NULL
  attr(x,'history') <- history.stamp(x)
  return(x)
}

as.eof.field <- function(x,...,iapp=NULL) {
  y <- EOF(x,...)
  if (!is.null(iapp)) y <- as.eof.appendix(y,iapp=iapp)
  return(y)
}

as.eof.appendix <- function(x,...,iapp=1,verbose=FALSE) {
  if (verbose) print("as.eof.appendix")
  clim <- eval(parse(text=paste("attr(attr(x,'appendix.",iapp,"'),'climatology')",sep="")))
  aveg <- eval(parse(text=paste("attr(attr(x,'appendix.",iapp,"'),'mean')",sep="")))
  stopifnot(inherits(x,'comb'))
  y <- eval(parse(text=paste("attr(x,'appendix.",iapp,"')",sep="")))
  x <- as.eof.comb(x)
  y <- attrcp(x,y)
  if (!is.null(clim)) attr(y,'climatology') <- clim 
  if (!is.null(aveg)) attr(y,'mean') <- aveg
  attr(y,'history') <- history.stamp(x)
  class(y) <- class(x)
  return(y)
}

as.eof.list <- function(x,...,verbose=FALSE) {
  stopifnot(inherits(x,'list'),inherits(x[[1]],'eof'))
  if (verbose) print('as.eof.list')
  
  wPC <- function(z,iapp=NULL) {
    eigv <- attr(z,'eigenvalues')
    w <- eigv/sum(eigv)
    if (is.null(iapp)) Z <- z %*% diag(w) else
                       Z <- attr(z,paste('appendix.',iapp,sep='')) %*% diag(w)
    Z <- zoo(Z,order.by=index(z))
    return(Z)
  }

  if (verbose) try(print(summary(x)))
  if (inherits(x[[1]],'character')) x[[1]] <- NULL
  if (inherits(x[[1]],'eof')) {eof <- x[[1]]; x[[1]] <- NULL}
  X.list <- lapply(x,wPC)
  X <- do.call("merge", X.list)
  if (verbose) print(summary(X))
  t <- index(X)
  udv <- svd(coredata(X))
  eof <- zoo(udv$u[,1:20],order.by=t)
  attr(eof,'eigenvalues') <- udv$d
  pattern <- rep(1,dim(udv$v)[1])
  names(pattern) <- names(X)
  attr(eof,'pattern') <- pattern
  if (inherits(x[[1]],'comb')) {
    if (verbose) print('Combined field: appendix.1')
    for (i in 1:attr( attr(x[[1]],'n.apps'))) {
      z.list <- lapply(x,wPC,iapp=i)
      udv1 <- svd(coredata(do.call("merge", z.list)))
      attr(eof,paste('appendix.',i,sep='')) <- zoo(udv1$u[,1:20],
               order.by=index(attr(x,paste('appendix.',i,sep=''))))
      names(attr(eof,paste('appendix.',i,sep=''))) <- paste("X.",1:20,sep="")
    }
  }
  attr(eof,'original.list.of.eofs') <- x
  attr(eof,'udv') <- udv
  id <- c()
  for (i in 1:length(x)) id <- c(id,rep(i,length(attr(x[[i]],'eigenvalues'))))
  attr(eof,'id') <- id
  names(eof) <- paste("X.",1:20,sep="")
  class(eof) <- class(x[[1]])
  return(eof)
}

as.eof.dsensemble <- function(x,...,FUN='mean',verbose=FALSE) {
  ## R.E. Benestad, 2017-05-19
  ## Convert the dsensemble object to an EOF of the multi-model mean
  stopifnot(inherits(x,'dsensemble'),inherits(x[[2]],'eof')|inherits(x[[2]],'pca'))
  if (verbose) print('as.eof.dsensemble')
  eof0 <- x[[2]]; x[[2]] <- NULL
  x[[1]] -> info; x[[1]] <- NULL
  d <- c(dim(x[[1]]),length(x))
  y <- unlist(x)
  dim(y) <- c(d[1]*d[2],d[3])
  Y <- apply(y,1,FUN)
  dim(Y) <- c(d[1],d[2])
  eof <- zoo(Y,order.by=index(x[[1]]))
  eof <- attrcp(eof0,eof)
  class(eof) <- class(eof0)
  attr(eof,'info') <- info
  attr(eof,'history') <- history.stamp()
  return(eof)
}



as.appended <- function(x,...) UseMethod("as.appended")

as.appended.ds.comb <- function(x,...,iapp=1,verbose=FALSE) {
  if(verbose) print("as.appended.ds.comb")
  eval(parse(text=paste("X <- attr(x,'appendix.",iapp,"')",sep="")))
  X <- attrcp(x,X,ignore='appendix')
  attr(X,'history') <- history.stamp(x)
  invisible(X)
}

as.appended.eof.comb <- function(x,...,iapp=1) {
  X <- as.appended.ds.comb(x,iapp=iapp)
  invisible(X)
}

as.appended.field.comb <- function(x,...,iapp=1) {
  X <- as.appended.ds.comb(x,iapp=iapp)
  invisible(X)
}

as.stand <- function(x,...) UseMethod("as.stand")

as.stand.station <- function(x,...,verbose=FALSE,na.rm=TRUE) {
  if(verbose) print("as.stand.station")
  if (is.precip(x)) {
    mu <- apply(x,2,mean,na.rm=na.rm)
    X <- 100*x/mu
    attr(X,'clim') <- mu
    attr(X,'aspect') <- 'proportional'
    attr(X,'unit') <- '%'
    attr(X,'oldunit') <- attr(x,'unit')
  } else if (is.T(x)) {
    mu <- apply(x,2,mean,na.rm=na.rm)
    sigma <- apply(x,2,sd,na.rm=na.rm)
    X <- (x - mu)/sigma
    attr(X,'mean') <- mu
    attr(X,'sigma') <- sigma
    attr(X,'aspect') <- 'standardised'    
  }
  attr(X,'history') <- history.stamp(x)
  return(X)
}



as.original <- function(x) UseMethod("as.original")

as.original.station <- function(x) {
  if (attr(x,'aspect')=='proportional') {
    X <- attr(x,'clim')*x/100
    attr(X,'clim') <- NULL
    attr(X,'unit') <- attr(x,'oldunit')
    attr(X,'oldunit') <- NULL
    attr(X,'aspect') <- 'original'
  } else if (attr(x,'aspect')=='standardised') {
    X <- x * attr(x,'sigma') + attr(x,'mean')
    attr(X,'mean') <- NULL
    attr(X,'sigma') <- NULL
    attr(X,'aspect') <- 'original'     
  } else X <- x
  attr(X,'history') <- history.stamp(x)
  return(X)
}

as.events <- function(x,...) UseMethod("as.events")

as.events.default <- function(x,...,label=NULL,dx=NULL,dy=NULL,
                      units=NULL,longname=NULL,variable=NULL,calendar=NULL,
                      qflabel=NULL,method=NULL,src=NULL,reference=NULL,
                      file=NULL,version=NULL,url=NULL,verbose=FALSE) {
  if (verbose) print("as.events")
  X <- data.frame(x)
  n <- names(X)
  if (!all(c("date","time","lon","lat") %in% names(X))) {
    print(paste("Missing input:",
     names(X)[!c("date","time","lon","lat")%in%names(X)]))
  }
  attr(X,"label") <- label
  attr(X,"dx") <- dx
  attr(X,"dy") <- dy
  attr(X,"longname") <- longname
  attr(X,"calendar") <- calendar
  attr(X,"variable") <- variable
  attr(X,"quality") <- qflabel
  attr(X,"units") <- units
  attr(X,"source") <- src
  attr(X,"file") <- file
  attr(X,"version") <- version
  attr(X,"method") <- method
  attr(X,"URL") <- url
  attr(X,"reference") <- reference
  class(X) <- c("events",class(X))
  attr(X,"history") <- history.stamp(X)
  invisible(X)
}

as.events.trajectory <- function(x,...,verbose=FALSE) {
  if(verbose) print("as.events.trajectory")
  stopifnot(inherits(x,"trajectory"))
  invisible(x)
}

as.field.events <- function(x,...) {
  y <- events2field(x,...)
  return(y)
}

as.field.trajectory <- function(x,...) {
  y <- trajectory2field(x,...)
  return(y)
}

as.station.trajectory <- function(x,...) {
  y <- trajectory2station(x,...)
  return(y)
}
