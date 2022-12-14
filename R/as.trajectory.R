#' Transform an input object into a \code{trajectory} object
#'
#' @param x input object
#' @param \dots additional arguments
#'
#' @aliases as.trajectory as.trajectory.default as.trajectory.events events2trajectory trajectory
#'
#' @export as.trajectory
as.trajectory <- function(x,...) UseMethod("as.trajectory")

#' @exportS3Method
#' @export
as.trajectory.default <- function(x,verbose=FALSE,...) {
  if(verbose) print("as.trajectory.default")
  X <- trajectory(x,...)
  invisible(X)
}

#' @exportS3Method
#' @export
as.trajectory.events <- function(x,verbose=FALSE,...) {
  if (verbose) print("as.trajectory.events")
  stopifnot(inherits(x,"events"))
  if (!("trajectory" %in% names(x))) x <- track(x,verbose=verbose,...)
  if (!("tracklength" %in% names(x))) x <- trackstats(x,verbose=verbose)
  y <- trajectory(x,verbose=verbose,...)
  invisible(y)
}

#' @export
events2trajectory <- function(x,verbose=FALSE,...) {
  if(verbose) print("events2trajectory")
  y <- as.trajectory.events(x,verbose=verbose,...)
  invisible(y)
}

#' @export
trajectory <- function(x,verbose=FALSE,loc=NA,param=NA,longname=NA,
                          quality=NA,src=NA,url=NA,reference=NA,info=NA,
                          method=NA,unit=NA,nmin=5,n=15) {
  if (verbose) print("trajectory")
  if (verbose) print(paste('dim: ',paste(dim(x),collapse=" x ")))
  if (verbose) print(paste('names: ',paste(names(x),collapse=", ")))
  names(x) <- tolower(names(x))
  names(x)[grep("latitude",names(x))] <- "lat"
  names(x)[grep("longitude",names(x))] <- "lon"
  names(x)[grep("step",names(x))] <- "timestep"
  if(is.na(loc) & !is.null(attr(x,"loc"))) loc <- attr(x,"loc")
  if(is.na(param) & !is.null(attr(x,"variable"))) param <- attr(x,"variable")
  if(is.na(longname) & !is.null(attr(x,"longname"))) longname <- attr(x,"longname")
  if(is.null(attr(x,"calendar"))) calendar <- "gregorian" else calendar <- attr(x,"calendar")
  if(is.na(quality) & !is.null(attr(x,"quality"))) quality <- attr(x,"quality")
  if(is.na(unit) & !is.null(attr(x,"unit"))) unit <- attr(x,"unit")
  if(is.na(src) & !is.null(attr(x,"source"))) src <- attr(x,"source")
  if(is.na(url) & !is.null(attr(x,"URL"))) url <- attr(x,"URL")
  if(is.na(reference) & !is.null(attr(x,"reference"))) reference <- attr(x,"refe
rence")
  if(is.na(info) & !is.null(attr(x,"info"))) info <- attr(x,"info")
  if(is.na(method) & !is.null(attr(x,"method"))) method <- attr(x,"method")
  
  # transform data in data.frame x to numeric values
  if(verbose) print('data.frame to numeric values')
  x <- data.frame(x)
  fn <- function(x) {
    if(!is.null(levels(x))) {suppressWarnings(as.numeric(levels(x))[x])
    } else as.numeric(as.character(x))
  }
  for (name in names(x)) {
    eval(parse(text=paste("x$",name,"<-fn(x$",name,")",sep="")))
  }
  
  if(verbose) print("remove short trajectories")
  #if (!("trackcount" %in% names(x))) {
  x <- trackstats(x,verbose=verbose)
  #}
  nlist1 <- c("trajectory","code99","date","time","trackcount",
              "start","end","n","d","distance","tracklength","timestep","lon","lat")
  nlist1 <- nlist1[nlist1 %in% names(x)]
  nlist2 <- names(x)[!(names(x) %in% nlist1)]
  
  x <- subset.events(x,it=x$trackcount>=nmin)
  X <- matrix(,nrow=length(unique(x$trajectory)),ncol=5+2*n+length(nlist2)*n)
  cnames <- c("trajectory","start","end","n","d",rep("lon",n),rep("lat",n))
  if(!is.null(nlist2)) cnames <- c(cnames,as.vector(sapply(nlist2,function(x) rep(x,n))))
  colnames(X) <- cnames
  if(dim(X)[1]==0) {
    if(verbose) print(paste("No trajectories >",nmin-1,"time steps"))
  } else {
    x["trajectory"] <- enumerate(x)
    # interpolate all trajectories to same length n
    if(verbose) print('interpolate trajectories')
    aggregate(x$date, list(x$trajectory), function(x) x[1])$x -> d1
    aggregate(x$time, list(x$trajectory), function(x) x[1])$x -> t1
    aggregate(x$date, list(x$trajectory), function(x) x[length(x)])$x -> d2
    aggregate(x$time, list(x$trajectory), function(x) x[length(x)])$x -> t2
    ## KMP 2017-10-04: Solution to problem with 06 time step when using format HHMM:
    if(max(c(t1,t2))>24) {
      t1 <- t1*1E-2
      t2 <- t2*1E-2
    }
    dt1 <- as.numeric(paste(d1,sapply(t1,function(h) if(h<10) paste("0",h,sep="") else as.character(h)),sep=""))
    dt2 <- as.numeric(paste(d2,sapply(t2,function(h) if(h<10) paste("0",h,sep="") else as.character(h)),sep=""))
    #aggregate(x$date, list(x$trajectory), length)$x -> len
    aggregate(x$trackcount, list(x$trajectory), function(x) x[1])$x -> nlen
    aggregate(x$distance, list(x$trajectory), function(x) x[1])$x -> dlen
    #aggregate(x$tracklength, list(x$trajectory), function(x) x[1])$x -> dlen
    a <- 6.378e06
    xx <- a * cos( x$lat*pi/180 ) * cos( x$lon*pi/180 )
    yy <- a * cos( x$lat*pi/180 ) * sin( x$lon*pi/180 )
    zz <- a * sin( x$lat*pi/180 )
    xa <- aggregate(xx, list(x$trajectory), function(x) approx(x,n=n)$y)$x
    ya <- aggregate(yy, list(x$trajectory), function(x) approx(x,n=n)$y)$x
    za <- aggregate(zz, list(x$trajectory), function(x) approx(x,n=n)$y)$x
    lon <- atan2( ya, xa )*180/pi 
    lat <- asin( za/sqrt( xa^2 + ya^2 + za^2 ))*180/pi
    colnames(lon) <- rep('lon',n)
    colnames(lat) <- rep('lat',n)
    X[,"trajectory"] <- unique(x$trajectory)
    X[,"start"] <- dt1
    X[,"end"] <- dt2
    X[,"n"] <- nlen
    X[,"d"] <- dlen
    X[,colnames(X)=="lon"] <- lon
    X[,colnames(X)=="lat"] <- lat
    for(name in nlist2) {
      if(verbose) print(name)
      eval(parse(text=paste("aggregate(x$",name,
        ",list(x$trajectory),function(x) approx(x,n=",
        n,")$y)$x -> y",sep="")))
      eval(parse(text=paste("colnames(y)<-rep('",name,"',",n,")",sep="")))
      eval(parse(text=paste("X[,colnames(X)=='",name,"'] <- y",sep="")))
    }
  }
  # add attributes to trajectory matrix X
  attr(X, "location")= loc
  attr(X, "variable")= param
  attr(X, "longname")= longname
  attr(X, "quality")= quality
  attr(X, "calendar")= calendar
  attr(X, "source")= src
  attr(X, "URL")= url
  attr(X, "unit")= unit
  attr(X, "type")= "analysis"
  attr(X, "aspect")= "interpolated"
  attr(X, "reference")= reference
  attr(X, "info")= info
  attr(X, "method")= method
  attr(X,"lon") <- NA
  attr(X,"lat") <- NA
  attr(X,"alt") <- NA
  attr(X,"cntr") <- NA
  attr(X,"stid") <- NA
  attr(X, "history")= history.stamp()
  class(X) <- c('trajectory','matrix')
  invisible(X)
}

#' Transform an input object into a \code{station} object  
#' 
#' Transform a \code{trajectory} object into a \code{station} object by aggregating it in time and space.
#'
#' @param x a \code{trajectory} object
#' @param it a time index, e.g., a range of years: c(1984,2019)
#' @param is a spatial index, e.g., a list with longitude and latitude ranges: list(lon=c(0,45), lat=c(45,70))
#' @param param a characteristic of the trajectories (to see options: \code{colnames(x)})
#' @param FUN a function. If 'count' return number of trajectories, otherwise apply \code{FUN} to \code{param}
#' @param longname variable name
#' @param unit name of unit
#' @param loc name of location/region
#'
#' @seealso as.station as.station.trajectory
#'
#' @export
trajectory2station <- function(x,it=NULL,is=NULL,param=NULL,FUN='count',
                               longname=NULL,unit=NULL,loc=NULL) {
  y <- subset(x,it=it,is=is)
  if (FUN=='count') {
    z <- count.trajectory(y,by='month')
    shortname <- 'trajectories'
    unit <- 'events/month'
  } else {
    z <- param.trajectory(y,param=param,FUN=FUN)
    t <- strptime(y[,colnames(y)=='start'],"%Y%m%d%H")
    z <- zoo(z,order.by=t)
    z <- aggregate(z,by=as.Date(t),FUN='mean')
    shortname <- paste(param,FUN,sep="")
  }
  z <- attrcp(y,z)
  attr(z,'variable') <- shortname
  attr(z,'unit') <- unit
  attr(z,'long name') <- longname
  attr(z,'location') <- loc
  attr(z,'longitude') <- lon(y)
  attr(z,'latitude') <- lat(y)
  z <- as.station.zoo(z)
  invisible(z)
}

#' Transform an input object into a \code{field} object  
#' 
#' Transform a \code{trajectory} object into a \code{field} object by aggregating it in time and space.
#'
#' @param x a \code{trajectory} object
#' @param dt frequency of output: 'month', 'season', 'quarter' (same as 'season') or 'year'
#' @param dx resolution in longitude direction (unit: degrees)
#' @param dy resolution in latitude direction (unit: degrees)
#' @param radius radius within which to look for trajectories for each grid point (unit: m)
#' @param it a time index, e.g., a range of years: c(1984,2019)
#' @param is a spatial index, e.g., a list with longitude and latitude ranges: list(lon=c(0,45), lat=c(45,70))
#' @param verbose a boolean; if TRUE print information about progress
#'
#' @return a \code{field} object
#'
#' @aliases as.field.trajectory
#' @seealso as.field as.trajectory CCI track.events
#'
#' @export
trajectory2field <- function(x,dt='month',dx=2,dy=2,radius=5E5,
                             it=NULL,is=NULL,verbose=FALSE) {
  if(verbose) print("trajectory2field")
  stopifnot(is.trajectory(x))
  x <- subset(x,it=it,is=is)
  lons <- round(x[,colnames(x)=='lon']/dx)*dx
  lons <- seq(min(lons),max(lons),dx)
  lats <- round(x[,colnames(x)=='lat']/dy)*dy
  lats <- seq(min(lats),max(lats),dy)
  dates <- strptime(x[,'start'],format='%Y%m%d%H')
  if(is.null(attr(x,"calendar"))) calendar <- "gregorian" else calendar <- attr(x,"calendar")
  if (requireNamespace("PCICt", quietly = TRUE)) {
    dates <- PCICt::as.PCICt(as.character(x[,'start']),format='%Y%m%d%H',cal=calendar)
  } else {
    dates <- as.POSIXct(as.character(x[,'start']),format='%Y%m%d%H')
  }
  fn <- function(x,it=NULL) {
    x <- subset(x,it=it)
    X.out <- array(rep(0,),dim=c(length(lons),length(lats)))
    if(length(x)>0) {
      if(verbose) {
        if(!is.null(dim(x))) {
          print(paste(range(x[,'start']),collapse="-"))
          print(dim(x))
        } else {
          print(x['start'])
        }
      }
      if(!is.null(dim(x))) {
        A <- trajectory2density(x,dx=dx,dy=dy,radius=radius)
        lat <- A$lat
        lon <- A$lon
        den <- A$density
        for(j in 1:length(lat)) {
          X.out[lons==lon[j],lats==lat[j]] <- den[j]
        }
        dim(X.out) <- c(length(lons),length(lats))#length(X)
      } else {
        if(verbose) print("no storms")
      }
    }
    invisible(X.out)
  }
  if(verbose) print("calculate trajectory density")
  if (grepl('month',dt)) {   
    if(verbose) print("monthly")
    dall <- as.Date(paste(format(dates,format="%Y-%m"),"01",sep="-"))
    dvec <- seq(min(dall),max(dall),by='month')
    unit <- 'events/month/area'
    fn2 <- function(di) fn(x,it=(dall==di))
    #X <- t(sapply(dvec,function(di) fn2))
  } else if (grepl('season',dt) | grepl('quarter',dt)) {
    if(verbose) print("seasonal")
    dall <- as.Date(as.yearqtr(as.Date(paste(format(dates,format="%Y-%m"),"01",sep="-"))))
    dvec <- seq(min(dall),max(dall),by='quarter')
    fn2 <- function(di) fn(x,it=(dall==di))
    #X <- t(sapply(dvec,fn2))
    unit <- 'events/quarter/area'
  } else if (grepl('year',dt) | grepl('annual',dt)) {
    if(verbose) print("annual")
    dall <- as.Date(paste(format(dates,format="%Y"),"01-01",sep="-"))
    dvec <- seq(min(dall),max(dall),by='year')
    unit <- 'events/year/area'
    fn2 <- function(di) fn(x,it=c(year(di),year(di)))
    #X <- t(sapply(dvec,fn2)
  } else {
    print(paste("WARNING! invalid time resolution dt",dt))
  }
  X <- array(rep(0,),dim=c(length(dvec),length(lons),length(lats)))
  if (verbose) print("looping...")
  t1 <- Sys.time()
  pb <- txtProgressBar(style=3)
  for (i in 1:length(dvec)) {
    setTxtProgressBar(pb,i/length(dvec))
    di <- try(fn2(dvec[i]))
    X[i,,] <- di
  }
  t2 <- Sys.time()
  if (verbose) print(paste('Calculating the density took',
                round(as.numeric(t2-t1,units="secs")),'s'))
  d <- dim(X)
  dim(X) <- c(d[1],d[2]*d[3])
  longname <- paste("trajectory density",attr(x,'longname'),sep=', ')
  param <- 'density'
  if(verbose) print("transform to field")
  Y <- as.field(X,index=dvec,lon=lons,lat=lats,
          unit=unit,longname=longname,param=param,
          quality=attr(x,'quality'),src=attr(x,'source'),
          url=attr(x,'URL'),reference=attr(x,'reference'),
          info=attr(x,'info'),calendar=attr(x,'calendar'),
          method=attr(x,'method'),aspect=attr(x,'aspect'))
  if(verbose) print("done")
  invisible(Y)
}

