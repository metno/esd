
as.trajectory <- function(x,...) {
  X <- trajectory(x,...)
  invisible(X)
}

trajectory <- function(x,verbose=FALSE,loc=NA,param=NA,longname=NA,
                          quality=NA,src=NA,url=NA,reference=NA,info=NA,
                          method=NA,n=10) {
  if (verbose) print("trajectory")
  if (verbose) print(paste('dim: ',paste(dim(x),collapse=" x ")))
  if (verbose) print(paste('names: ',paste(names(x),collapse=", ")))
  names(x) <- tolower(names(x))
  names(x)[grep("lat",names(x))] <- "lat"
  names(x)[grep("lon",names(x))] <- "lon"
  names(x)[grep("step",names(x))] <- "timestep"

  if(is.na(loc) & !is.null(attr(x,"loc"))) loc <- attr(x,"loc")
  if(is.na(param) & !is.null(attr(x,"variable"))) param <- attr(x,"variable")
  if(is.na(unit) & !is.null(attr(x,"unit"))) unit <- attr(x,"unit")
  if(is.na(longname) & !is.null(attr(x,"longname"))) longname <- attr(x,"longname")
  if(is.na(quality) & !is.null(attr(x,"quality"))) quality <- attr(x,"quality")
  if(is.na(src) & !is.null(attr(x,"source"))) src <- attr(x,"source")
  if(is.na(url) & !is.null(attr(x,"URL"))) url <- attr(x,"URL")
  if(is.na(reference) & !is.null(attr(x,"reference"))) reference <- attr(x,"reference")
  if(is.na(info) & !is.null(attr(x,"info"))) info <- attr(x,"info")
  if(is.na(method) & !is.null(attr(x,"method"))) method <- attr(x,"method")
  
  # transform data in data.frame x to numeric values
  x <- data.frame(x)
  fn <- function(x) {
    if(!is.null(levels(x))) {suppressWarnings(as.numeric(levels(x))[x])
    } else as.numeric(as.character(x))
  }
  nlist1 <- c('trajectory','lat','lon','year','month','day','time')
  if (sum(nlist1 %in% names(x))==length(nlist1)) {
    x$trajectory <- fn(x$trajectory)
    x$lat <- fn(x$lat)
    x$lon <- fn(x$lon)
    x$year <- fn(x$year)
    x$month <- fn(x$month)
    x$day <- fn(x$day)
    x$time <- fn(x$time)
    x$date <- x$year*1E6+x$month*1E4+x$day*1E2+x$time
  } else {
    print(paste(paste(nlist1[!(nlist1 %in% names(x))],collapse=' '),'missing'))
  }
  nlist2 <- names(x)[!(names(x) %in% c('code99','date',nlist1))]
  for (name in nlist2) {
    eval(parse(text=paste("x$",name,"<-fn(x$",name,")",sep="")))
  }

  # interpolate all trajectories to same length n
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

  # add attributes to trajectory matrix X
  attr(X, "location")= loc
  attr(X, "variable")= param
  attr(X, "longname")= longname
  attr(X, "quality")= quality
  attr(X, "calendar")= "gregorian"
  attr(X, "source")= src
  attr(X, "URL")= url
  attr(X, "type")= "analysis"
  attr(X, "aspect")= "interpolated"
  attr(X, "reference")= reference
  attr(X, "info")= info
  attr(X, "method")= method
  attr(x,"lon") <- NA
  attr(x,"lat") <- NA
  attr(x,"alt") <- NA
  attr(x,"cntr") <- NA
  attr(x,"stid") <- NA
  attr(X, "history")= history.stamp()
  class(X) <- 'trajectory'
  invisible(X)
}


read.imilast <- function(fname,path=NULL) {
  fname <- paste(path,fname,sep="")
  # read file header
  h <- strsplit(readLines(fname,1),",")
  h <- tolower(unlist(h))
  h <- gsub("^\\s+|\\s+$","",h)
  h[grep("cyclone",h)] <- "trajectory"
  h[grep("99",h)] <- "code99"
  h[grep("step",h)] <- "timestep"
  h[grep("date",h)] <- "date"
  h[grep("lat",h)] <- "lat"
  h[grep("lon",h)] <- "lon"
  # check width of columns
  l <- readLines(fname,4)[4]
  blanks <- unlist(gregexpr(" ",l))
  breaks <- blanks[!(blanks %in% as.integer(blanks+1))]
  w <- c(blanks[1]-1,
         breaks[2:length(breaks)]-breaks[1:(length(breaks)-1)],
         nchar(l)-breaks[length(breaks)]+1)
  # read data
  x <- read.fwf(fname,width=w,col.names=h,skip=1)
  x <- x[x$code99<90,]
  # add attributes
  attr(x, "variable")= "storm tracks"
  attr(x, "longname")= "mid-latitude storm trajectories"
  attr(x, "source")= "IMILAST"
  attr(x, "URL")= "http://journals.ametsoc.org/doi/abs/10.1175/BAMS-D-11-00154.1"
  attr(x, "reference")= "Neu, et al. , 2013: IMILAST: A Community Effort to Intercompare Extratropical Cyclone Detection and Tracking Algorithms. Bull. Amer. Meteor. Soc., 94, 529–547."
  if (x$code99[1]>9) {
    attr(x, "method") <- paste("M",as.character(x$code99[1]),sep="")
  } else attr(x, "method") <- paste("M0",as.character(x$code99[1]),sep="")
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
    shortname <- paste(param,FUN,sep=".")
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

  ## y <- subset(x,it=it,is=is)

  ## lon <- min(y[,colnames(y)=='lon'])
  ## lon <- seq(min(lon(y)),max(lon(y)),length=30)
  ## lat <- seq(min(lat(y)),max(lat(y)),length=30)

  ## z <- c()
  ## z <- attrcp(y,z)
  ## as.field.zoo(z,lon,lat,param,unit)
  
  
}

