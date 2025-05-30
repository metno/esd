#' Anomaly and Climatology
#' 
#' S3-method that computes anomalies and/or climatology for time series and
#' fields.
#'
#' In 'anomaly.dsensemble', the default value of the reference period is taken
#' as the available time period from observations, i.e., same time period as in
#' attribute `station' is used as base period to compute anomalies of GCM
#' downscaled results.
#'
#' @aliases anomaly anomaly.default anomaly.comb anomaly.field anomaly.station 
#' anomaly.annual anomaly.month anomaly.season anomaly.day anomaly.dsensemble
#' as.anomaly as.anomaly.default as.anomaly.zoo as.anomaly.list as.anomaly.station as.anomaly.field
#' climatology as.climatology
#' @seealso as.stand
#' 
#' @param x A station or field object
#' @param ref vector defining the reference interval
#' @param na.rm a boolean; if TRUE remove NA values
#' @param verbose a boolean; if TRUE print information about progress
#' @param \dots additional arguments
#'
#' @return a similar object as x containing anomalies and climatology
#'
#' @keywords utilities
#'
#' @examples 
#' data(ferder)
#' plot(anomaly(ferder))
#' 
#' 
#' @export
anomaly <-function(x,...) UseMethod("anomaly")

#' @exportS3Method
#' @export anomaly.default
anomaly.default <- function(x,...) {

  ### REB 2023-07-17
  arguments <<- c(as.list(environment()), list(...))
  #if (!is.null(arguments$na.rm)) ref <- arguments$ref else ref <- rep(TRUE,length(index(x)))
  
  if (!is.null(arguments$verbose)) verbose <- arguments$verbose else verbose <- FALSE
  if (verbose) print('anomaly.default')
  if (!is.null(arguments$na.rm)) na.rm <- arguments$na.rm else na.rm <- TRUE
  if (!is.null(arguments$verbose)) verbose <- arguments$verbose else verbose <- FALSE
  
  if(verbose) {
    print('anomaly.default')
    print(match.call())
    print(class(x))
  }
      
  if (!is.null(arguments$ref)) {
    ref <- arguments$ref
    it <- is.element(year(x),ref) 
  } else {
    ref <- unique(year(x))
    it <- rep(TRUE,length(index(x)))
  }
  if (verbose) print(paste('ref=',year(x)[it][1],'-',year(x)[it][sum(it)],
                           'n=',sum(it)))
  if (inherits(x,'annual')) {
    y <- anomaly.annual(x,...) 
  } else if (inherits(x,'month')) {
    y <- anomaly.month(x,...) 
  } else if (inherits(x,'day')) {
    y <- anomaly.day(x,...) 
  } else if (inherits(x,'season')) {
    y <- anomaly.season(x,...) 
  } else if (is.null(dim(x))) {
    y <- x - mean(x[it],na.rm=TRUE)
    attr(y, "mean") <- mean(x[it],na.rm=TRUE)
  } else {
    y <- zoo(apply(x,2,function(x,it) x - mean(x[it],na.rm=TRUE),it), order.by=index(x))
  }
  y <- attrcp(x,y)
  class(y) <- class(x)
  # KMP 2019-02-12: attr(y,'station') sometimes gets the station_id instead of station
  # so !is.null(attr...) is not enough to check if the attribute is available
  if ('station' %in% names(attributes(y)) & !is.null(attr(y,'station'))) {
    attr(y,'station') <- anomaly(attr(y,'station'), ref=ref)
  }
  attr(y,'aspect') <- 'anomaly' 
  attr(y,'history') <- history.stamp(x)
  return(y)
}

#' @exportS3Method
#' @export anomaly.dsensemble
anomaly.dsensemble <- function(x,...) {
  ### REB 2023-07-17
  arguments <<- c(as.list(environment()), list(...))
  ref <- arguments$ref
  if (!is.null(arguments$na.rm)) na.rm <- arguments$na.rm else na.rm <- TRUE
  if (!is.null(arguments$verbose)) verbose <- arguments$verbose else verbose <- FALSE
    if(verbose) print("anomaly.dsensemble")
    yr.obs <- year(attr(x,'station'))
    ref <- range(yr.obs[!is.na(attr(x,'station'))],na.rm=TRUE)
    stopifnot(inherits(x,"dsensemble"))
    x <- anomaly.default(x,...)
    attr(x,'station') <- anomaly(attr(x,'station'),ref=ref,...)
    return(x)
}

#' @exportS3Method
#' @export anomaly.field
anomaly.field <- function(x,...,verbose=FALSE) {
  ### REB 2023-07-17
  arguments <<- c(as.list(environment()), list(...))
  ref <- arguments$ref
  if (!is.null(arguments$na.rm)) na.rm <- arguments$na.rm else na.rm <- TRUE
  if (!is.null(arguments$verbose)) verbose <- arguments$verbose else verbose <- FALSE
  if(verbose) print("anomaly.field")
  stopifnot(inherits(x,"field"))
 
  x <- as.anomaly(x,...)
  return(x)
}

#' @exportS3Method
#' @export anomaly.comb
anomaly.comb <- function(x,...) {
  arguments <<- c(as.list(environment()), list(...))
  ref <- arguments$ref
  if (!is.null(arguments$na.rm)) na.rm <- arguments$na.rm else na.rm <- TRUE
  if (!is.null(arguments$verbose)) verbose <- arguments$verbose else verbose <- FALSE
  if(verbose) print("anomaly.comb")
  stopifnot(inherits(x,"field"),inherits(x,"comb"))
  arguments <<- c(as.list(environment()), list(...))
  ref <- arguments$ref
  if (is.null(arguments$na.rm)) na.rm <- arguments$na.rm else na.rm <- TRUE
  y <- anomaly(x)
  n.apps <- attr(x,'n.apps')
  for (i in 1:n.apps) {
    z <- NULL
    eval(parse(text=paste("z <- attr(x,'appendix.",i,"')",sep="")))
    Z <- anomaly(z)
    eval(parse(text=paste("Z -> attr(x,'appendix.",i,"')",sep="")))
  }
  y <- attrcp(x,y)
  n.apps -> attr(x,'n.apps')
  attr(x,'history') <- history.stamp(x)
  return(y)
}

#' @exportS3Method
#' @export anomaly.station
anomaly.station <- function(x, ...) {
  arguments <<- c(as.list(environment()), list(...))
  ref <- arguments$ref
  if (!is.null(arguments$na.rm)) na.rm <- arguments$na.rm else na.rm <- TRUE
  if (!is.null(arguments$verbose)) verbose <- arguments$verbose else verbose <- FALSE
  if(is.null(ref)) ref <- seq(min(year(x)),max(year(x)),by=1)
  if(verbose) { 
    print("anomaly.station")
    print(match.call())
  }
  x <- anomaly.default(x,verbose=verbose,...)
  return(x)
}

#' @exportS3Method
#' @export anomaly.annual
anomaly.annual <- function(x,...) {
  arguments <<- c(as.list(environment()), list(...))
  # print('anomaly.annual'); print(names(arguments))
  # str(arguments)
  ref <- arguments$ref
  if (!is.null(arguments$na.rm)) na.rm <- arguments$na.rm else na.rm <- TRUE
  if (!is.null(arguments$verbose)) verbose <- arguments$verbose else verbose <- FALSE
  if(is.null(ref)) ref <- seq(min(year(x)),max(year(x)),by=1)
  if (verbose) {
    print('anomaly.annual')
    print(match.call())
    print(range(ref)); print(dim(x))
  }
  X <- x;  x <- coredata(X)
  t <- index(X)
  d <- dim(x)
  if (!is.null(d)) ns <- d[2] else ns <- 1
  datetype <- class(t)
  if (verbose) print(datetype)
  if (datetype=="Date") years <- year(X) else
  if (datetype=="numeric") years <- t
  if (ns==1) { 
    if (verbose) print('single series')
    clim <- mean(x[is.element(years,ref)],na.rm=TRUE) 
  } else {
    if (verbose) print('field')
    clim <- apply(x[is.element(years,ref),],2,FUN=mean,na.rm=TRUE) 
  }
  if (verbose) print(clim)
  x <- t(t(x) - clim)
  x <- zoo(x,order.by=index(X))
  x <- attrcp(X,x)
  #nattr <- softattr(X)
  #for (i in 1:length(nattr))
  #  attr(x,nattr[i]) <- attr(X,nattr[i])
  attr(x,'climatology') <- clim
  attr(x,'aspect') <- 'anomaly'
  attr(x,'history') <- history.stamp(X)
  class(x) <- class(X)
  return(x)
}

#' @exportS3Method
#' @export anomaly.month
anomaly.month <- function(x,...) {
  arguments <<- c(as.list(environment()), list(...))
  ref <- arguments$ref
  if (!is.null(arguments$na.rm)) na.rm <- arguments$na.rm else na.rm <- TRUE
  if (!is.null(arguments$verbose)) verbose <- arguments$verbose else verbose <- FALSE
  if(verbose) print("anomaly.month")
  clim.month <- function(x,months,years,ref=NULL,na.rm=TRUE,verbose=FALSE) {
    ## This function calculated the mean for each calendar month.
    if (verbose) cat('.')
    if (!is.null(ref)) {
      it <- is.element(years,c(min(ref):max(ref))) 
    } else if (na.rm) {
      it=is.finite(x) 
    } else {
      it <- rep(TRUE,length(x))
    }
    clim <- rep(NA,12)
    for (i in 1:12) {
      ii <- is.element(months,i) & it
      clim[i] <- mean(x[ii],na.rm=TRUE)
    }
    names(clim) <- month.abb
    return(clim)
  }

  X <- x
  if (verbose) print(paste('anomaly.month: ref=',ref))
  if (is.null(dim(x))) {
    if (verbose) print('Estimate clim for a single series')
    Yc <- clim.month(x,months=month(x),years=year(x),ref=ref,na.rm=na.rm,verbose=verbose)
    Y <- X-Yc[month(x)] 
  } else {
    if (verbose) print('Estimate clim for multiple seies')
    Yc <- apply(coredata(x),2,FUN='clim.month',months=month(x),years=year(x),
                ref=ref,na.rm=na.rm,verbose=verbose)
    Y <- X-Yc[month(x),] 
  }
    
  y <- Y
  x <- zoo(y,order.by=index(X))
  x <- attrcp(X,x)
  attr(x,'climatology') <- Yc
  attr(x,'aspect') <- 'anomaly'
  class(x) <- class(X)
  return(x)
}

#' @exportS3Method
#' @export anomaly.season
anomaly.season <- function(x,...) {
  arguments <<- c(as.list(environment()), list(...))
  ref <- arguments$ref
  if (!is.null(arguments$na.rm)) na.rm <- arguments$na.rm else na.rm <- TRUE
  if (!is.null(arguments$verbose)) verbose <- arguments$verbose else verbose <- FALSE
  if (verbose) {print('anomaly.season'); print(match.call())}
  
  anomaly.season1 <- function(x,yr=NULL,ref=NULL,verbose=FALSE,what='anomaly') {
    l <- length(x); n <- ceiling(l/4)
    pad <- 4*n - l
    if (verbose) print(paste('anomaly.season1: l=',l,' n=',n,' pad=',pad))
    
    ## base-line period
    if (!is.null(yr) & !is.null(ref)) iref <- is.element(yr,ref) else
                                      iref <- rep(TRUE,n)
    ## If the record is not full years, pad the extra months of the last year
    if (pad>0) x <- c(rep(NA,pad),x)
    ##Fast way to compute the climatology: clim
    dim(x) <- c(4,n)
    clim <- rowMeans(x,na.rm=TRUE)
    x <- c(x - clim)
    if (pad>0) x <- x[-(1:pad)]
    if (substr(what,1,4)=='clim') x <- clim
    return(x)
  }
  anomaly.season1.v2 <- function(x,ref=NULL,t=NULL,verbose=FALSE,what='anomaly') {
    x0 <- x
    yr <- year(t)
    mn <- month(t)
    l <- length(x)
    m <- length(unique(mn))
    n <- length(unique(yr))
    pad <- m*n - l
    if (pad>0) x <- c(rep(NA,pad),x)
    dim(x) <- c(m, n)
    clim <- rowMeans(x, na.rm=TRUE)
    x <- c(x - clim)
    if (pad>0) x <- x[-(1:pad)]
    if (substr(what,1,4)=='clim') x <- clim
    return(x)
  }
  X <- x
  t <- index(x); yr <- year(x)
  if (is.null(ref)) ref=seq(min(yr),max(yr),by=1)
  if (verbose) print(paste('anomaly.season: ref=',min(ref),'-',max(ref)))
  if (is.null(dim(x))) {
    #y <- anomaly.season1(coredata(x),yr=yr,ref=ref,verbose=verbose)
    #clim <- anomaly.season1(coredata(x),yr=yr,ref=ref,verbose=verbose,what='clim')
    y <- anomaly.season1.v2(coredata(x),t=index(x),verbose=verbose,ref=ref)
    clim <- anomaly.season1.v2(coredata(x),t=index(x),verbose=verbose,what='clim',ref=ref)
  } else {
    #y <- apply(coredata(x),2,FUN='anomaly.season1',yr=yr,ref=ref,verbose=verbose)
    #clim <- apply(coredata(x),2,FUN='anomaly.season1',yr=yr,ref=ref,verbose=verbose,what='clim')
    y <- apply(x,2,FUN='anomaly.season1.v2',t=index(x),verbose=verbose,ref=ref)
    clim <- apply(x,2,FUN='anomaly.season1.v2',t=index(x),verbose=verbose,what='clim',ref=ref)
    if(is.null(dim(clim))) dim(clim) <- c(1, length(clim))
  }
  x <- zoo(y,order.by=t)
  x <- attrcp(X,x)
  clim <- zoo(clim, order.by=unique(month(index(x))))
  attr(x,'climatology') <- clim
  attr(x,'aspect') <- 'anomaly'
  class(x) <- class(X)
  return(x)
}

#' @export anomaly.day
anomaly.day <- function(x,...) {
  arguments <<- c(as.list(environment()), list(...))
  ref <- arguments$ref
  if (!is.null(arguments$verbose)) verbose <- arguments$verbose else verbose <- FALSE
  if (verbose) {
    print('anomaly.day')
  }
  if (!is.null(arguments$na.rm)) na.rm <- arguments$na.rm else na.rm <- TRUE
  if (!is.null(arguments$verbose)) verbose <- arguments$verbose else verbose <- FALSE
  ## Internal function
  anomaly.day.1 <- function(x,t0,t,ref=NULL,verbose=FALSE) {
    if (verbose) print(' anomaly.day.1')
    ## One station 
    c1 <- cos(pi*t0/365.25); s1 <- sin(pi*t0/365.25)
    c2 <- cos(2*pi*t0/365.25); s2 <- sin(2*pi*t0/365.25)
    c3 <- cos(3*pi*t0/365.25); s3 <- sin(3*pi*t0/365.25)
    c4 <- cos(4*pi*t0/365.25); s4 <- sin(4*pi*t0/365.25)
    C1 <- cos(pi*t/365.25);   S1 <- sin(pi*t/365.25)
    C2 <- cos(2*pi*t/365.25); S2 <- sin(2*pi*t/365.25)
    C3 <- cos(3*pi*t/365.25); S3 <- sin(3*pi*t/365.25)
    C4 <- cos(4*pi*t/365.25); S4 <- sin(4*pi*t/365.25)
    if (is.null(dim(x))) x0 <- x[ref] else x0 <- x[ref,]
    cal <- data.frame(y=x0,c1=c1,c2=c2,c3=c3,c4=c4,
                      s1=s1,s2=s2,s3=s3,s4=s4)
    if (verbose) str(cal)
    pre <- data.frame(c1=C1,c2=C2,c3=C3,c4=C4,
                      s1=S1,s2=S2,s3=S3,s4=S4)
    if (verbose) str(cal)
    ## The following lines don't seem to make much difference REB 2023-12-19
    # i1 <- is.element(year(x),year(x)[1])
    # pre1 <- data.frame(c1=C1[i1],c2=C2[i1],c3=C3[i1],c4=C4[i1],
    #                    s1=S1[i1],s2=S2[i1],s3=S3[i1],s4=S4[i1])
    acfit <- lm(y ~ c1 + s1 + c2 + s2 + c3 + s3 + c4 + s4,data=cal)
    clim <- predict(acfit,newdata=pre)
    y <- zoo(coredata(x) - clim,order.by=index(x))
    if (verbose) print('exit  anomaly.day.1')
    return(y)
  }
  
  if (verbose) {print('anomaly.day');print(class(x))}
  yr <- year(x)
  if (is.null(ref)) ref <- seq(min(yr,na.rm=TRUE),max(yr,na.rm=TRUE),by=1)
  if (verbose) {print('reference:'); print(range(ref))}
  ## Time indices with julian time repeated for each year during reference
  t0 <- julian(index(x)[is.element(yr,ref)]) -
    julian(as.Date(paste(yr[is.element(yr,ref)],"-01-01",sep="")))
  ## Time indices with julian time repeated for each year for entire data record
  t <- julian(index(x)) - julian(as.Date(paste(yr,"-01-01",sep="")))
  if (verbose) {print(c(length(t),length(t0)))}
  ## ref is TRUE or FALSE
  ref <- is.element(yr,ref)
  if (verbose) print(paste('Length of reference period=',sum(ref)))
  if (is.null(dim(x))) 
    y <- anomaly.day.1(x=coredata(x),t0=t0,t=t,ref=ref,verbose=verbose) else 
    y <- apply(coredata(x),2,FUN='anomaly.day.1',t0=t0,t=t,ref=ref,verbose=verbose)
  y <- zoo(y,order.by=index(x))
  y <- attrcp(x,y)
  class(y) <- class(x)
  ## find the first year with complete data (365 valid points)
  z <- y[is.finite(y)]
  fullyear <- table(year(z)) == 365
  completeyear <- as.numeric(rownames(table(year(z)))[fullyear])
  clim <- (x - y)[is.element(yr,completeyear[1])] 
  if (verbose) print(paste('Climatology has',sum(is.finite(clim)),'valid data'))
  attr(y,'climatology') <- clim
  attr(y,'aspect') <- 'anomaly'
  return(y)
}


#' @export as.anomaly
as.anomaly <- function(x,...) UseMethod("as.anomaly")

#' @exportS3Method
#' @export
as.anomaly.default <- function(x,...) anomaly.default(x,...)

#' @exportS3Method
#' @export as.anomaly.zoo
as.anomaly.zoo <- function(x,...) {
  arguments <<- c(as.list(environment()), list(...))
  ref <- arguments$ref
  if (!is.null(arguments$na.rm)) na.rm <- arguments$na.rm else na.rm <- TRUE
  if (!is.null(arguments$verbose)) verbose <- arguments$verbose else verbose <- FALSE
  y <- as.anomaly.station(x,...)
  attr(y,'history') <- history.stamp(x)
  invisible(y)
}

#' @exportS3Method
#' @export as.anomaly.list
as.anomaly.list <- function(x,...) {
  arguments <<- c(as.list(environment()), list(...))
  ref <- arguments$ref
  if (!is.null(arguments$na.rm)) na.rm <- arguments$na.rm else na.rm <- TRUE
  if (!is.null(arguments$verbose)) verbose <- arguments$verbose else verbose <- FALSE
  y <- lapply(x,anomaly(x),ref=ref,na.rm=na.rm,verbose=verbose)
  attr(y,'history') <- history.stamp(x)
  invisible(y)
}

#' @exportS3Method
#' @export as.anomaly.station
as.anomaly.station <- function(x,...) {
  y <- as.anomaly.default(x,...)
  attr(y,'history') <- history.stamp(x)
  invisible(y)
}

#' @exportS3Method
#' @export as.anomaly.field
as.anomaly.field <- function(x,...) {
   y <- anomaly.default(x,...)
   attr(y,'history') <- history.stamp(x)
   attr(y,'dimensions') <- attr(x,'dimensions')
   invisible(y)
}

#' @export
climatology <- function(x,...) {
  x <- as.climatology(x,...)
  return(x)
}

# Handy conversion algorithms:
#' @export as.climatology
as.climatology <- function(x,...) {
  if (is.null(attr(x, "aspect"))) {
    warning('esd::as.climatolody - aspect not set, but assumed to be "measured"')
    attr(y,'aspect') <- 'measured'
  }
  if(attr(x, "aspect")=="anomaly") ya <- x else ya <- as.anomaly(x,...)
  clim <- attr(ya,'climatology')
  if(inherits(clim,"zoo")) {
    y <- clim
  } else {
    if (!is.null(dim(clim))) {
      len.clim <- dim(clim)[1]
    } else {
      len.clim <- length(clim)
    }
    y <- zoo(clim,order.by=1:len.clim)
  }
  y <- attrcp(x,y)
  attr(y,'aspect') <- 'climatology'
  attr(y,'history') <- history.stamp(x)
  class(y) <- class(x)
  invisible(y)
}

