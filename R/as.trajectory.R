as.trajectory <- function(x,...) UseMethod("as.trajectory")

as.trajectory.default <- function(x,verbose=FALSE,...) {
  if(verbose) print("as.trajectory.default")
  X <- trajectory(x,...)
  invisible(X)
}

as.trajectory.events <- function(x,verbose=FALSE,...) {
  if (verbose) print("as.trajectory.events")
  stopifnot(inherits(x,"events"))
  if (!("trajectory" %in% names(x))) x <- track(x,verbose=verbose,...)
  if (!("tracklength" %in% names(x))) x <- Trackstats(x,verbose=verbose)
  y <- trajectory(x,verbose=verbose,...)
  invisible(y)
}

events2trajectory <- function(x,verbose=FALSE,...) {
  if(verbose) print("events2trajectory")
  y <- as.trajectory.events(x,verbose=verbose,...)
  invisible(y)
}

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
  x <- Trackstats(x,verbose=verbose)
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
    x["trajectory"] <- Enumerate(x)
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
    dt1 <- as.numeric(strftime(strptime(paste(d1,t1),format="%Y%m%d %H"),format="%Y%m%d%H"))
    dt2 <- as.numeric(strftime(strptime(paste(d2,t2),format="%Y%m%d %H"),format="%Y%m%d%H"))
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
    #aggregate(x$lat, list(x$trajectory),function(x) approx(x,n=n)$y)$x -> lat
    #aggregate(x$lon,list(x$trajectory),function(x) approxlon(x,n=n)$y)$x -> lon
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
  attr(X, "calendar")= "gregorian"
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

matrix2tdf <- function(x) {
  i.t <- colnames(x) %in% c('start','end')
  i.n <- colnames(x)=='n'  
  pnames <- unique(colnames(x))
  pnames <- pnames[!(pnames %in% c('start','end','n'))]
  trajectory <- unlist(mapply(function(i,n) rep(i,n),seq(1,dim(x)[1]),x[,i.n]))
  for(p in pnames) {
    i.p <- colnames(x)==p
    if (p=='lon') {
      lon <- unlist(apply(x,1,function(x) approxlon(x[i.p],n=x[i.n])$y))
    } else {
      eval(parse(text=paste(p,
       " <- unlist(apply(x,1,function(x) approx(x[i.p],n=x[i.n])$y))",sep="")))
    }
  }
  fn <- function(x) {
    d <- strptime(x[i.t],format="%Y%m%d%H")
    d <- seq(d[1],d[2],length=x[i.n])
    d <- as.numeric(strftime(d,"%Y%m%d%H"))
    invisible(d)
  }
  date <- unlist(apply(x,1,fn))
  year <- as.integer(date*1E-6)
  code99 <- rep(as.numeric(substring(attr(x,"method"),2)),length(date))
  eval(parse(text=paste('X <- data.frame(trajectory,date,year,',
               paste(pnames,collapse=","),',code99)',sep="")))
  invisible(X)
}


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
