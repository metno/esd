#' Aggregate to seasonal values
#'
#' @aliases as.4seasons as.4seasons.default as.4seasons.station as.4seasons.day as.4seasons.spell as.4seasons.field as.4seasons.dsensemble
#'
#' @seealso aggregate as.annual as.monthly
#'
#' @param x input object, e.g., a 'station' or 'field' object
#' @param FUN a function
#' @param verbose a boolean; if TRUE print information about progress
#' @param nmin minimum number of data points in a season
#' @param slow a boolean; if FALSE run a fast version that might not work on all data types?
#' @param dateindex a boolean; if TRUE transform index into date format
#' @param \dots additional arguments
#'
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

#' @export
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
