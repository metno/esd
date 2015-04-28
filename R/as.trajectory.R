
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
  aggregate(x$Date, list(x$Trajectory), function(x) x[1])$x -> t1
  aggregate(x$Date, list(x$Trajectory), function(x) x[length(x)])$x -> t2
  aggregate(x$Date, list(x$Trajectory), length)$x -> len
  aggregate(x$Lat, list(x$Trajectory),function(x) approx(x,n=n)$y)$x -> lat
  aggregate(x$Lon,list(x$Trajectory),function(x) approxlon(x,n=n)$y)$x -> lon
  colnames(lon) <- rep('lon',n)
  colnames(lat) <- rep('lat',n)
  if(verbose) print(n[1:n])
  nlist1 <- c('Trajectory','Lat','Lon','Year','Month','Day','Time',
             'Date','Code99','timeStep')
  nlist2 <- names(x)[!(names(x) %in% nlist1)]
  if(is.null(nlist2)) {
    X <- cbind(lon=lon,lat=lat,start=t1,end=t2,n=len)
  } else {
    for(name in nlist2) {
      if(verbose) print(name)
      eval(parse(text=paste("aggregate(x$",name,
       ",list(x$Trajectory),function(x) approx(x,n=",
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

tdf2numeric <- function(x,verbose=verbose) {
  x <- data.frame(x)
  if ('Code99' %in% names(x)) {
    x$Code99 <- as.numeric(as.character(x$Code99))
    x <- x[x$Code99<90,]
  }
  nlist1 <- c('Trajectory','Lat','Lon','Year','Month','Day','Time')
  if (sum(nlist1 %in% names(x))==length(nlist1)) {
    x$Trajectory <- factor2numeric(x$Trajectory)
    x$Lat <- factor2numeric(x$Lat)
    x$Lon <- factor2numeric(x$Lon)
    x$Year <- factor2numeric(x$Year)
    x$Month <- factor2numeric(x$Month)
    x$Day <- factor2numeric(x$Day)
    x$Time <- factor2numeric(x$Time)
    x$Date <- x$Year*1E6+x$Month*1E4+x$Day*1E2+x$Time
  } else {
    print(paste(paste(nlist1[!(nlist1 %in% names(x))],collapse=' '),'missing'))
  }
  nlist2 <- names(x)[!(names(x) %in% c('Code99','Date',nlist1))]
  for (name in nlist2) {
    eval(parse(text=paste("x$",name,"<-factor2numeric(x$",name,")",sep="")))
  }
  invisible(x)
}

factor2numeric <- function(x) {
  suppressWarnings(as.numeric(levels(x))[x])}


trajectory2station <- function(x) {

}

trajectory2field <- function(x) {

}


## trajectory2station <- function(x,it=NULL,is=NULL,loc=NULL,FUN='count') {
##   y <- subset(x,it=it,is=is)
##   if (FUN=='count') {
##     n <- count.trajectory(y,by='month') # monthly storm count?
##   } else if (FUN=='lat1') {
##     n <- y[,which(colnames(y)=='lat')[1]]
##     attr(n,'long name') <- 'genesis latitude'
##     attr(n,'unit') <- 'degrees N'
##   } else if (FUN=='lon1') {
##     n <- y[,which(colnames(y)=='lon')[1]]
##     attr(n,'long name') <- 'genesis longitude'
##     attr(n,'unit') <- 'degrees E'
##   }
##   if (!inherits(n,'zoo')) {
##     t <- strptime(y[,colnames(y)=='start'],"%Y%m%d%H")
##     n <- zoo(n,order.by=t)
##     n <- aggregate(n,by=as.Date(t)) # calculate average?
##   }
##   ns <- as.station(n,loc=loc,param='storms',unit='count',
## 		lon=lon(n),lat=lat(n))
##   class(ns) <- c("station",class(n))
##   invisible(ns)
## }
