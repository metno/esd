#' Transform an 'event' object into a 'station' object
#'
#' Aggregate some aspects of an 'events' object in time and space and transform the time series into a 'station' object.
#'
#' If the aggregated aspect is the number of events (param='count', the default option), the events can be aggregated in several different ways.
#' The frequency of aggregation can be set by the input parameter 'by' (default: "month").
#' If by.trajectory=FALSE (default), the function counts the average number of events per timestep at the selected interval.
#' If by.trajectory=TRUE, the function instead counts the number of individual trajectories that are identified within that time period.
#' The latter option requires 'track' to have been applied to the events.
#'
#' @param x input object of class 'events'
#' @param param parameter to aggregate, e.g., 'count' or some characteristic such as 'pcent' or 'radius' (for options, see names(x))
#' @param FUN a function, e.g., 'mean' or 'sum'
#' @param by.trajectory a boolean; This input is only used when param is 'count'. See description above.
#' @param by the time interval of aggregation, only used when param is 'count'. The default is "month", other options are "year", "day", and "timestep".
#'        When param is not 'count', by cannot be changed and the aggregation is always on a monthly time scale. 
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
    y <- matrix(rep(NA, as.numeric(nrow(x))*as.numeric(length(dt))), 
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
    N <- count.events(x,FUN=FUN,...,verbose=verbose)
    longname <- paste(attr(x,"variable"),param)
    unit <- attr(N, "unit")
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
    index(N) <- as.Date(strptime(index(N), format="%Y-%m-%d"))
    N <- subset(N, it=paste(range(format(dates, format="%Y-%m")),"01",sep="-"))
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

#' Count the number of events per month
#'
#' @param x input object of class 'events'
#' @param by.trajectory if TRUE count every trajectory once, otherwise count every time step separately
#' @param verbose if TRUE print progress
#' @param \dots additional arguments
#'
#' @export count.events
count.events <- function(x,by.trajectory=TRUE,FUN=NULL,by="month",dhr=NULL,verbose=FALSE,...) {
  if (verbose) print("count.events")
  if(is.null(attr(x,"calendar"))) calendar <- "gregorian" else calendar <- attr(x,"calendar")
  if (requireNamespace("PCICt", quietly = TRUE)) {
    dates <- PCICt::as.PCICt(paste(x$date,x$time),format="%Y%m%d %H",cal=calendar)
    #dates <- PCICt::as.PCICt(x$date,format="%Y%m%d",cal=calendar)
    if(by=="month") {
      fn <- function(x) PCICt::as.PCICt(paste(format(x,"%Y-%m"),"01",sep="-"),cal=calendar)
    } else if(by=="day") {
      fn <- function(x) PCICt::as.PCICt(paste(format(x,"%Y-%m-%d"),sep="-"),cal=calendar)
    } else if(by=="year") {
      fn <- function(x) PCICt::as.PCICt(paste(format(x,"%Y"),sep="-"),cal=calendar)
    } else if(by=="timestep") {
      fn <- function(x) PCICt::as.PCICt(x, cal=calendar)
    }
  } else {
    dates <- as.Date(strptime(x$date,format="%Y%m%d"))
    if(by=="month") {
      fn <- function(x) as.Date(as.yearmon(x))
    } else if(by=="year") {
      fn <- function(x) as.Date(as.year(x))
    } else if(by=="day") {
      fn <- function(x) x
    } else if(by=="timestep") {
      dates <- strptime(paste(x$date, x$time), format="%Y%m%d %H")
      fn <- function(x) x
    }
  }
  if (by.trajectory & !(by!="timestep")) {
    if (!"trajectory" %in% names(x)) x <- track(x)
    z <- zoo(x$trajectory,order.by=dates)
    N <- aggregate(z, by=fn, FUN=function(x) length(unique(x,na.rm=TRUE)))
    z2 <- zoo(paste(x$date,x$time), order.by=dates)
    #N_timesteps <- aggregate(z2, by=fn, FUN=function(x) length(unique(x)))
  } else {
    ## KMP 2024-10-09: Changing the function for individual untracked cyclones
    ## so that you get the mean count per time step rather than the sum (which will depend on the length of each month)
    #z <- zoo(x$date,order.by=dates)
    z <- zoo(paste(x$date,x$time), order.by=dates)
    N <- aggregate(z,by=fn,FUN=length)
    #N_timesteps <- aggregate(z, by=fn, FUN=function(x) length(unique(x)))
    #N_days <- aggregate(z2, by=fn, FUN=function(x) length(unique(x)))
    #N <- N/N2_days
  }
  # fill in missing months by merging with an empty time series
  if (requireNamespace("PCICt", quietly = TRUE)) {
    nrt <- PCICt::as.PCICt(as.character(range(year(dates))*1E4+range(month(dates))*1E2+1),format="%Y%m%d",cal=calendar)
  } else {
    nrt <- as.Date(strptime(range(year(dates))*1E4+range(month(dates))*1E2+1,format="%Y%m%d"))
  }
  if(by=="timestep") {
    N0 <- zoo(,seq(from = nrt[1], to = nrt[2], by = difftime(index(N)[2], index(N)[1]) ))
  } else {
    N0 <- zoo(,seq(from = nrt[1], to = nrt[2], by = by))
  }
  N <- merge(N, N0)
  N[is.na(N)] <- 0
  if(by=="timesteps") N_timesteps <- rep(1, length(N0)) else {
    if(is.null(dhr)) { dhr <- diff(x$time); dhr <- min(dhr[dhr>0]) }
    nrt <- as.Date(strptime(c(min(year(dates)), max(year(dates))+1)*1E4 + 101, format="%Y%m%d"))
    N0 <- zoo(,seq(from = nrt[1], to = nrt[2], by = by))
    dday <- difftime(index(N0)[2:length(index(N0))], 
                     index(N0)[1:(length(index(N0))-1)])
    N_timesteps <- as.numeric(dday*(24/dhr))
  }
  if(!is.null(FUN)) if(FUN=="mean") {
    cb <- coredata(N)
    cb <- cb/N_timesteps
    coredata(N) <- cb
  }
  ## ============================================================================
  ## KMP 2026-01-14: The package PCICt, which deals with non-standard time formats 
  ## (e.g 360-day calendars), is deprecated and should be removed from esd.
  ## It is causing trouble here and there down the line. However, we don't have a good 
  ## replacement if we want esd to work for unphysical calendar data. 
  ## For now, I am adding some line here to avoid returning a timeseries with 
  ## a PCICt index when possible (i.e. monthly/annual resolution and when the calendar is normal)
  if(inherits(index(N), "PCICt")) {
    if(by %in% c("month","year")) {
      index(N) <- as.Date(strptime(index(N), format="%Y-%m-%d"))
    } else if(!attr(x, "calendar")=="360_day") {
      index(N) <- as.POSIXct(index(N))
    }
  }
  ## ============================================================================
  N <- attrcp(x, N)
  N <- as.station(N)
  attr(N, "timesteps") <- N_timesteps
  unit <- paste0("events/", by)
  if(!is.null(FUN)) if(FUN=="mean" & !by.trajectory) unit <- "events/timestep"
  attr(N, "unit") <- unit
  invisible(N)
}