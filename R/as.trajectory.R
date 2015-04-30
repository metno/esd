
as.trajectory <- function(x,verbose=FALSE,loc=NA,param=NA,unit=NA,lon=NA,
                          lat=NA,alt=NA,cntr=NA,longname=NA,stid=NA,
                          quality=NA,src=NA,url=NA,reference=NA,info=NA,
                          method= NA,n=10) {
  if (verbose) print("as.trajectory")
  if (verbose) print(paste('dim: ',paste(dim(x),collapse=" x ")))
  if (verbose) print(paste('names: ',paste(names(x),collapse=", ")))
  X <- tdf2matrix(x,n=n,verbose=verbose)
  attr(X, "location")= loc
  attr(X, "variable")= param
  attr(X, "unit")= unit
  attr(X, "longitude")= lon
  attr(X, "latitude")= lat
  attr(X, "altitude")= alt
  attr(X, "country")= cntr
  attr(X, "longname")= longname
  attr(X, "station_id")= stid
  attr(X, "quality")= quality
  attr(X, "calendar")= "gregorian"
  attr(X, "source")= src
  attr(X, "URL")= url
  attr(X, "type")= "analysis"
  attr(X, "aspect")= "interpolated"
  attr(X, "reference")= reference
  attr(X, "info")= info
  attr(X, "method")= method
  attr(X, "history")= history.stamp()
  class(X) <- 'trajectory'
  invisible(X)
}

tdf2matrix <- function(x,n=10,verbose=FALSE) {
  x <- tdf2numeric(x)
  aggregate(x$date, list(x$trajectory), function(x) x[1])$x -> t1
  aggregate(x$date, list(x$trajectory), function(x) x[length(x)])$x -> t2
  aggregate(x$date, list(x$trajectory), length)$x -> len
  aggregate(x$lat, list(x$trajectory),function(x) approx(x,n=n)$y)$x -> lat
  aggregate(x$lon,list(x$trajectory),function(x) approxlon(x,n=n)$y)$x -> lon
  colnames(lon) <- rep('lon',n)
  colnames(lat) <- rep('lat',n)
  if(verbose) print(n[1:n])
  nlist1 <- c('trajectory','lat','lon','year','month','day','time',
             'date','code99','timestep')
  nlist2 <- names(x)[!(names(x) %in% nlist1)]
  if(is.null(nlist2)) {
    X <- cbind(lon=lon,lat=lat,start=t1,end=t2,n=len)
  } else {
    for(name in nlist2) {
      if(verbose) print(name)
      eval(parse(text=paste("aggregate(x$",name,
       ",list(x$trajectory),function(x) approx(x,n=",
                   n,")$y)$x ->",name,sep="")))
      eval(parse(text=paste("colnames(",name,
                   ")<-rep('",name,"',",n,")",sep="")))
    }
    X <- eval(parse(text=paste("cbind(lon=lon,lat=lat",
        paste(nlist2,"=",nlist2,sep="",collapse=","),
        "start=t1,end=t2,n=len)",sep=",")))
  }
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

tdf2numeric <- function(x,verbose=verbose) {
  x <- data.frame(x)
  names(x) <- tolower(names(x))
  if ('code99' %in% names(x)) {
    x$code99 <- as.numeric(as.character(x$code99))
    x <- x[x$code99<90,]
  }
  nlist1 <- c('trajectory','lat','lon','year','month','day','time')
  if (sum(nlist1 %in% names(x))==length(nlist1)) {
    x$trajectory <- factor2numeric(x$trajectory)
    x$lat <- factor2numeric(x$lat)
    x$lon <- factor2numeric(x$lon)
    x$year <- factor2numeric(x$year)
    x$month <- factor2numeric(x$month)
    x$day <- factor2numeric(x$day)
    x$time <- factor2numeric(x$time)
    x$date <- x$year*1E6+x$month*1E4+x$day*1E2+x$time
  } else {
    print(paste(paste(nlist1[!(nlist1 %in% names(x))],collapse=' '),'missing'))
  }
  nlist2 <- names(x)[!(names(x) %in% c('code99','date',nlist1))]
  for (name in nlist2) {
    eval(parse(text=paste("x$",name,"<-factor2numeric(x$",name,")",sep="")))
  }
  invisible(x)
}

factor2numeric <- function(x) {
  suppressWarnings(as.numeric(levels(x))[x])}


trajectory2station <- function(x,it=NULL,is=NULL,param=NULL,FUN='count',
                               longname=NULL,unit=NULL,loc=NULL) {
  y <- subset(x,it=it,is=is)
  if (FUN=='count') {
    z <- count.trajectory(y,by='month')
    shortname <- 'trajectory count'
    unit <- 'events/month'
  } else {
    z <- param.trajectory(y,param=param,FUN=FUN)
    t <- strptime(y[,colnames(y)=='start'],"%Y%m%d%H")
    z <- zoo(z,order.by=t)
    z <- aggregate(z,by=as.Date(t))
    shortname <- paste(param,FUN)
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

trajectory2field <- function(x,it=NULL,is=NULL,param=NULL,FUN='count',
                             unit=NULL) {

  y <- subset(x,it=it,is=is)
  
  lon <- seq(min(lon(y)),max(lon(y)),length=30)
  lat <- seq(min(lat(y)),max(lat(y)),length=30)

  z <- c()
  z <- attrcp(y,z)
  as.field.zoo(z,lon,lat,param,unit)
  
  
}

