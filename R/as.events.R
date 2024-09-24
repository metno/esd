#' Transform an input object into an \code{events} object
#'
#' @param x input object
#' @param verbose if TRUE print progress
#' @param \dots additional arguments
#'
#' @aliases as.events as.events.default as.events.trajectory events
#' trajectory2events
#'
#' @export as.events
as.events <- function(x,verbose=FALSE,...) UseMethod("as.events")

#' @exportS3Method
#' @export
as.events.default <- function(x,verbose=FALSE,...) {
  if(verbose) print("as.events.default")
  X <- events(x,verbose=verbose,...)
  invisible(X)
}

#' @exportS3Method
#' @export
as.events.trajectory <- function(x,verbose=FALSE,...) {
  if (verbose) print("as.events.trajectory")
  stopifnot(inherits(x,"trajectory"))
  y <- trajectory2events(x,verbose=verbose,...)
  invisible(y)
}

#' @export
events <- function(x,verbose=FALSE,loc=NULL,param=NULL,longname=NULL,calendar=NULL,
                   quality=NULL,src=NULL,url=NULL,reference=NULL,greenwich=NULL,
                   info=NULL,method=NULL,unit=NULL,file=NULL,version=NULL) {
  if (verbose) print("events")
  if (verbose) print(paste('dim: ',paste(dim(x),collapse=" x ")))
  if (verbose) print(paste('names: ',paste(names(x),collapse=", ")))
  if(inherits(x,"matrix")) {
    cnames <- colnames(x)
    x <- data.frame(x)
    names(x) <- cnames
  }
  if(is.null(loc) & !is.null(attr(x,"loc"))) loc <- attr(x,"loc")
  if(is.null(param) & !is.null(attr(x,"variable"))) param <- attr(x,"variable")
  if(all(is.null(unit)) & !is.null(attr(x,"unit"))) unit <- attr(x,"unit")
  if(is.null(version) & !is.null(attr(x,"version"))) version <- attr(x,"version")
  if(is.null(longname) & !is.null(attr(x,"longname"))) longname <- attr(x,"longname")
  if(is.null(calendar) & !is.null(attr(x,"calendar"))) calendar <- attr(x,"calendar")
  if(is.null(quality) & !is.null(attr(x,"quality"))) quality <- attr(x,"quality")
  if(is.null(src) & !is.null(attr(x,"source"))) src <- attr(x,"source")
  if(is.null(file) & !is.null(attr(x,"file"))) src <- attr(x,"file")
  if(is.null(url) & !is.null(attr(x,"URL"))) url <- attr(x,"URL")
  if(is.null(reference) & !is.null(ref(x))) reference <- ref(x)
  if(is.null(info) & !is.null(attr(x,"info"))) info <- attr(x,"info")
  if(is.null(method) & !is.null(attr(x,"method"))) method <- attr(x,"method")
  names(x) <- tolower(names(x))
  names(x)[grep("latitude",names(x))] <- "lat"
  names(x)[grep("longitude",names(x))] <- "lon"
  names(x)[grep("step",names(x))] <- "timestep"
  names(x)[grep("^p$",names(x)) | grep("^slp$",names(x))] <- "pcent"
  if(is.null(greenwich)) {
    if(!is.null(attr(x,"greenwich"))) {
      greenwich <- attr(x,"greenwich")
    } else {
      lon <- as.numeric(x$lon)
      greenwich <- !(min(lon,na.rm=TRUE)<0 | max(lon,na.rm=TRUE)<=180)
    }
  }
  if (inherits(x,'trajectory')) {
    y <- trajectory2events(x,verbose=verbose)
  } else {
    y <- x
  }
  attr(y, "location") <- loc
  attr(y, "variable") <- param
  attr(y, "longname") <- longname
  attr(y, "calendar") <- calendar
  attr(y, "unit") <- unit
  attr(y, "quality") <- quality
  attr(y, "source") <- src
  attr(y, "URL") <- url
  attr(y, "type") <- "analysis"
  attr(y, "aspect") <- "interpolated"
  attr(y, "reference") <- reference
  attr(y, "info") <-info
  attr(y, "method") <- method
  attr(y,"lon") <- NA
  attr(y,"lat") <- NA
  attr(y,"alt") <- NA
  attr(y,"cntr") <- NA
  attr(y,"stid") <- NA
  attr(y, "history") <- history.stamp()
  class(y) <- c("events","data.frame")
  y <- g2dl(y, greenwich=greenwich, verbose=verbose)
  invisible(y) 
}

#' @export
trajectory2events <- function(x,minlen=3,verbose=FALSE) {
  stopifnot(inherits(x,"trajectory"))
  
  if (verbose) print("trajectory2events")
  if (verbose) print(paste('dim: ',paste(dim(x),collapse=" x ")))
  if (verbose) print(paste('names: ',paste(names(x),collapse=", ")))
  cnames <- unique(colnames(x))
  params <- cnames[!cnames %in% c("trajectory","lat","lon","start","end","n")]
  cnames <- c("trajectory","lon","lat","time","n",params)
  if(verbose) print(paste('remove trajectories shorter than',minlen))
  if(verbose) print('interpolate trajectories to original length')
  y <- data.frame(matrix(rep(NA,sum(x[,"n"])*length(cnames)),
                         nrow=sum(x[,"n"]),ncol=length(cnames)))
  names(y) <- cnames
  i.n <- colnames(x)=="n"
  i.lat <- colnames(x)=="lat"
  i.lon <- colnames(x)=="lon"
  if(!"trajectory" %in% colnames(x)) {
    tr <- seq(nrow(x))
  } else {
    tr <- x[,"trajectory"]
  }
  y$trajectory <- unlist(mapply(rep, tr, x[,i.n]))
  y$n <- unlist(apply(x,1,function(z) rep(z[i.n],z[i.n])))
  if(verbose) print('interpolate time')
  i.start <- colnames(x)=="start"
  i.end <- colnames(x)=="end"
  datetime <- unlist(apply( x, 1, function(z) format(
        seq(strptime(z[i.start],format="%Y%m%d%H"),
        strptime(z[i.end],format="%Y%m%d%H"),
        length.out=z[i.n]),format="%Y%m%d%H")))
  y$date <- as.numeric(format(strptime(datetime,"%Y%m%d%H"),"%Y%m%d"))
  y$time <- as.numeric(format(strptime(datetime,"%Y%m%d%H"),"%H"))
  if(verbose) print('interpolate lon and lat')
  a <- 6.378e06
  xx <- a * cos( x[,i.lat]*pi/180 ) * cos( x[,i.lon]*pi/180 )
  yy <- a * cos( x[,i.lat]*pi/180 ) * sin( x[,i.lon]*pi/180 )
  zz <- a * sin( x[,i.lat]*pi/180 )
  xa <- unlist(apply(cbind(xx,x[,i.n]),1,
               function(z) approx(z[1:(length(z)-1)],n=z[length(z)])$y))
  ya <- unlist(apply(cbind(yy,x[,i.n]),1,
               function(z) approx(z[1:(length(z)-1)],n=z[length(z)])$y))
  za <- unlist(apply(cbind(zz,x[,i.n]),1,
               function(z) approx(z[1:(length(z)-1)],n=z[length(z)])$y))
  y$lon <- atan2( ya, xa )*180/pi
  y$lat <- asin( za/sqrt( xa^2 + ya^2 + za^2 ))*180/pi
  if (length(params)>0) {
    for (p in params) {
      if(verbose) print(paste('interpolate',p))
      i.p <- colnames(x)==p
      if(sum(i.p)==1) {
        y[p] <- unlist(apply(x,1,function(z) rep(z[i.p],z[i.n])))
      } else {
        y[p] <- unlist(apply(x,1,function(z) approx(z[i.p],n=z[i.n])$y))
      }
    }
  }
  y <- attrcp(x,y)
  class(y) <- c("events","data.frame")
  return(y)
}

#' Count the number of events per month
#'
#' @param x input object of class 'events'
#' @param by.trajectory if TRUE count every trajectory once, otherwise count every time step separately
#' @param verbose if TRUE print progress
#' @param \dots additional arguments
#'
#' @export count.events
count.events <- function(x,by.trajectory=TRUE,verbose=FALSE,...) {
  if (verbose) print("count.events")
  if(is.null(attr(x,"calendar"))) calendar <- "gregorian" else calendar <- attr(x,"calendar")
  if (requireNamespace("PCICt", quietly = TRUE)) {
    dates <- PCICt::as.PCICt(paste(x$date,x$time),format="%Y%m%d %H",cal=calendar)
    #dates <- PCICt::as.PCICt(x$date,format="%Y%m%d",cal=calendar)
    fn <- function(x) PCICt::as.PCICt(paste(format(x,"%Y-%m"),"01",sep="-"),cal=calendar)
  } else {
    dates <- as.Date(strptime(x$date,format="%Y%m%d"))
    fn <- function(x) as.Date(as.yearmon(x))
  }
  if (by.trajectory) {
    if (!"trajectory" %in% names(x)) x <- track(x)
    z <- zoo(x$trajectory,order.by=dates)
    N <- aggregate(z,by=fn,FUN=function(x) length(unique(x,na.rm=TRUE)))
  } else {
    z <- zoo(x$date,order.by=dates)
    N <- aggregate(z,by=fn,FUN=length)
  }
  # fill in missing months by merging with an empty time series
  if (requireNamespace("PCICt", quietly = TRUE)) {
    nrt <- PCICt::as.PCICt(as.character(range(year(dates))*1E4+range(month(dates))*1E2+1),format="%Y%m%d",cal=calendar)
  } else {
    nrt <- as.Date(strptime(range(year(dates))*1E4+range(month(dates))*1E2+1,format="%Y%m%d"))
  }
  N0 <- zoo(,seq(from = nrt[1], to = nrt[2], by = "month"))
  N <- merge(N, N0)
  N[is.na(N)] <- 0
  N <- attrcp(x,N)
  N <- as.station(N)
  invisible(N)
}


