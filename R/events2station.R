#' Transform an 'event' object into a 'station' object
#'
#' Aggregate some aspects of an 'events' object in time and space and transform the time series into a 'station' object.
#'
#' @param x input object of class 'events'
#' @param param parameter to aggregate, e.g., 'count' or some characteristic such as 'pcent' or 'radius' (for options, see names(x))
#' @param FUN a function, e.g., 'mean' or 'sum'
#' @param verbose a boolean; If TRUE print information about progress
#' @param longname long name of variable
#' @param unit name of unit
#' @param \dots additional arguments
#'
#' @return a 'station' object
#'
#' @export
events2station <- function(x,param="count",FUN="mean",verbose=FALSE,
                         longname=NULL,unit=NULL,...) {
  if (verbose) print("events2station")
  if(is.null(attr(x,"calendar"))) calendar <- "gregorian" else calendar <- attr(x,"calendar")
  if (requireNamespace("PCICt", quietly = TRUE)) {
    dates <- PCICt::as.PCICt(paste(x$date,x$time),format="%Y%m%d %H",cal=calendar)
    fn <- function(x) PCICt::as.PCICt(paste(format(x,"%Y-%m"),"01",sep="-"),cal=calendar)
  } else {
    dates <- as.POSIXct(paste(x$date,x$time),format="%Y%m%d %H")
    fn <- function(x) as.Date(as.yearmon(x))
  }
  if("location" %in% c(FUN,param)) {
    dt <- unique(dates)
    y <- matrix(rep(NA, nrow(x)*length(dt)), 
                ncol=nrow(x), nrow=length(dt))
    if(param %in% names(x)) j <- which(names(x)==param) else j <- 5
    for(i in seq(1,nrow(x))) y[dt==dates[i],i] <- x[i,j]
    N <- as.station(zoo(y, order.by=dt), 
                    lon=x$lon, lat=x$lat)#,
                    #variable=colnames(x)[j], unit=attr(x,"unit")[j])
    if("trajectory" %in% names(x)) attr(N, "trajectory") <- x$trajectory
    param <- colnames(x)[j]
    longname <- paste(FUN," of ",param,sep="")
    unit <- attr(x,"unit")[j]
    class(N) <- c("station", "hourly", "zoo")
  } else if (param=="count") {
    N <- count.events(x,...,verbose=verbose)
    longname <- paste(attr(x,"variable"),param)
    unit <- "events/months"
  } else if (param %in% names(x)) {
    y <- zoo(x[,param],order.by=dates)
    N <- aggregate(y,by=fn,FUN=FUN)
    if (requireNamespace("PCICt", quietly = TRUE)) {
      nrt <- PCICt::as.PCICt(as.character(range(year(dates))*1E4+range(month(dates))*1E2+1),format="%Y%m%d",cal=calendar)
    } else {
      nrt <- as.Date(strptime(range(year(dates))*1E4+range(month(dates))*1E2+1,format="%Y%m%d"))
    }
    N0 <- zoo(,seq(from = nrt[1], to = nrt[2], by = "month"))
    N <- merge(N, N0)
    N <- attrcp(x,N)
    N <- as.station(N)
    attr(N,"lat") <- attr(x,"lat")
    attr(N,"lon") <- attr(x,"lon")
    N <- subset(N, it=paste(range(format(dates,format="%Y-%m")),"01",sep="-"))
  } else {
    print(paste("input error: param =",param))
  }
  if (inherits(x,c("season","month"))) {
    mn <- unique(month(N[!is.na(N) & N>0]))
    N <- subset(N,it=month.abb[mn])
    class(N) <- c("station",class(x)[2],"zoo")
  }
  if (is.null(longname)) {
    longname <- paste(FUN," ",param," (",attr(x,"variable"),")",sep="")
  }
  if (is.null(unit)) {
    if (length(names(x))==length(attr(x,"unit"))) {
      unit <- attr(x,"unit")[names(x)==param]
    } else {
      unit <- NA
    }
  }
  # KMP 2024-09-24: I'm changing the time index class to Date which shouldn't be necessary, 
  #  but the POSIXlt and PCICt formats suddenly started causing problems in a plot function
  #  and this was the easiest solution I could think of.
  index(N) <- as.Date(strptime(index(N), format="%Y-%m-%d"))
  attr(N,"variable") <- param
  attr(N,"longname") <- longname
  attr(N,"calendar") <- calendar
  attr(N,"unit") <- unit
  #N <- subset(N, it=paste(range(strftime(dates,format="%Y-%m")),"01",sep="-"))
  invisible(N)
}