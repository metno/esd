as.trajectory <- function(x,...) UseMethod("as.trajectory")

as.trajectory.default <- function(x,...) {
  X <- trajectory(x,...)
  invisible(X)
}

as.trajectory.events <- function(x,verbose=FALSE,...) {
  if (verbose) print("as.trajectory.events")
  stopifnot(inherits(x,"events"))
  if (!("trajectory" %in% names(x))) x <- Track.events(x,verbose=verbose,...)
  if (!("trackcount" %in% names(x))) x <- Trackstats(x,verbose=verbose)
  y <- trajectory(x,verbose=verbose,...)
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
  x <- Trackstats(x)
  #}
  x <- subset.events(x,it=x$trackcount>nmin)
  x["trajectory"] <- Enumerate(x)
  
  # interpolate all trajectories to same length n
  if(verbose) print('interpolate trajectories')
  aggregate(x$date, list(x$trajectory), function(x) x[1])$x -> d1
  aggregate(x$time, list(x$trajectory), function(x) x[1])$x -> t1
  aggregate(x$date, list(x$trajectory), function(x) x[length(x)])$x -> d2
  aggregate(x$time, list(x$trajectory), function(x) x[length(x)])$x -> t2
  dt1 <- as.numeric(strftime(strptime(paste(d1,t1),format="%Y%m%d %H"),format="%Y%m%d%H"))
  dt2 <- as.numeric(strftime(strptime(paste(d2,t2),format="%Y%m%d %H"),format="%Y%m%d%H"))
  #aggregate(x$date, list(x$trajectory), length)$x -> len
  aggregate(x$trackcount, list(x$trajectory), function(x) x[1])$x -> len
  aggregate(x$tracklength, list(x$trajectory), function(x) x[1])$x -> distance
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
  #if(verbose) print(n[1:n])
  #nlist1 <- c('trajectory','lat','lon','year','month','day','time',
  #           'date','code99','timestep')
  #nlist2 <- names(x)[!(names(x) %in% nlist1)]
  nlist1 <- c("trajectory","code99","date","time","trackcount",
              "tracklength","timestep","lon","lat")
  nlist1 <- nlist1[nlist1 %in% names(x)]
  nlist2 <- names(x)[!(names(x) %in% nlist1)]
  if(is.null(nlist2)) {
    X <- cbind(trajectory=unique(x$trajectory),lon=lon,lat=lat,
               start=dt1,end=dt2,n=len,d=distance)
  } else {
    for(name in nlist2) {
      if(verbose) print(name)
      eval(parse(text=paste("aggregate(x$",name,
           ",list(x$trajectory),function(x) approx(x,n=",
           n,")$y)$x ->",name,sep="")))
      eval(parse(text=paste("colnames(",name,
           ")<-rep('",name,"',",n,")",sep="")))
    }
    X <- eval(parse(text=paste("cbind(trajectory=unique(x$trajectory),lon=lon,lat=lat",
        paste(nlist2,"=",nlist2,sep="",collapse=","),
        "start=dt1,end=dt2,n=len,d=distance)",sep=",")))
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


read.imilast <- function(fname,path=NULL,verbose=FALSE) {
  fname <- file.path(path,fname)
  # read and rearrange file header
  if(verbose) print(paste("reading file header:"))
  h <- tolower(readLines(fname,1))
  if(verbose) print(h)
  if(verbose) print("rearranging")
  h <- unlist(strsplit(h,","))
  if(length(h)==1) h <- unlist(strsplit(h," "))
  h <- gsub("^\\s+|\\s+$","",h)
  h[grepl("cyclone",h) | grepl("trackn",h)] <- "trajectory"
  h[grepl("lati",h) | grepl("latn",h)] <- "lat"
  h[grep("long",h)] <- "lon"
  h[grepl("99",h) | grepl("code",h)] <- "code99"
  h[grepl("date",h) | grepl("yyyymmddhh",h)] <- "datetime"
  h[grep("yyyy",h)] <- "year"
  h[grep("mm",h)] <- "month"
  h[grep("dd",h)] <- "day"
  h[grepl("hh",h) | grepl("timestep\\.",h) | grepl("timestep\\_",h)] <- "time"
  h[grepl("step",h) | grepl("ptn",h)] <- "timestep"
  h <- unique(h)
  if(verbose) print(paste(h,collapse=", "))
  # check width of columns
  l <- readLines(fname,4)[4]
  l <- gsub("^\\s+|\\s+$","",l)
  blanks <- unlist(gregexpr(" ",l))
  breaks <- blanks[!(blanks %in% as.integer(blanks+1))]
  w <- c(blanks[1]-1,
         breaks[2:length(breaks)]-breaks[1:(length(breaks)-1)],
         nchar(l)-breaks[length(breaks)]+1)
  if(verbose) print(paste("width of columns:",paste(w,collapse=",")))                
  # read data
  if(verbose) print("reading data")
  x <- read.fwf(fname,width=w,col.names=h,skip=1)
  x <- x[x$code99<90,]
  # rearrange date and time information
  dates <- round(x["datetime"][[1]]*1E-2)
  times <- x["datetime"][[1]] - round(dates)*1E2
  x["date"] <- dates
  x["time"] <- times
  x <- x[,-which(names(x) %in% c("datetime","year","month","day","timestep"))]
  #add attributes
  variable <- "storm tracks"
  longname <- "mid-latitude storm trajectories"
  src <- "IMILAST"
  url <- "http://journals.ametsoc.org/doi/abs/10.1175/BAMS-D-11-00154.1"
  ref <- "Neu, et al. , 2013: IMILAST: A Community Effort to Intercompare Extratropical Cyclone Detection and Tracking Algorithms. Bull. Amer. Meteor. Soc., 94, 529–547."
  if (x$code99[1]>9) {
    method <- paste("M",as.character(x$code99[1]),sep="")
  } else {
    method <- paste("M0",as.character(x$code99[1]),sep="")
  }
  x <- as.events(x,label="imilast",dx=min(diff(sort(unique(x["lon"][[1]])))),
                 dy=min(diff(sort(unique(x["lat"][[1]])))),
                 longname=longname,variable=variable,method=method,src=src,
                 reference=ref,file=file.path(path,fname),url=url)
  attr(x, "history")= history.stamp() 
  invisible(x)
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
