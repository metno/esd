as.events <- function(x,...) UseMethod("as.events")

as.events.default <- function(x,verbose=FALSE,...) {
  if(verbose) print("as.events.default")
  X <- events(x,verbose=verbose,...)
  invisible(X)
}

as.events.trajectory <- function(x,verbose=FALSE,...) {
  if (verbose) print("as.events.trajectory")
  stopifnot(inherits(x,"trajectory"))
  y <- trajectory2events(x,verbose=verbose,...)
  invisible(y)
}

events <- function(x,verbose=FALSE,loc=NULL,param=NULL,longname=NULL,
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
  if(is.null(quality) & !is.null(attr(x,"quality"))) quality <- attr(x,"quality")
  if(is.null(src) & !is.null(attr(x,"source"))) src <- attr(x,"source")
  if(is.null(file) & !is.null(attr(x,"file"))) src <- attr(x,"file")
  if(is.null(url) & !is.null(attr(x,"URL"))) url <- attr(x,"URL")
  if(is.null(reference) & !is.null(ref(x))) reference <- ref(x)
  if(is.null(info) & !is.null(attr(x,"info"))) info <- attr(x,"info")
  if(is.null(method) & !is.null(attr(x,"method"))) method <- attr(x,"method")
  if(is.null(greenwich)) {
    if(!is.null(attr(x,"greenwich"))) {
      greenwich <- attr(x,"greenwich")
    } else {
      greenwich <- !(min(x$lon)<0 | max(x$lon)<=180)
    }
  }
  if (inherits(x,'trajectory')) {
    y <- trajectory2events(x,verbose=verbose)
  } else {
    y <- x
  }
  names(y) <- tolower(names(y))
  names(y)[grep("latitude",names(y))] <- "lat"
  names(y)[grep("longitude",names(y))] <- "lon"
  names(y)[grep("step",names(y))] <- "timestep"
  names(y)[grep("^p$",names(y)) | grep("^slp$",names(y))] <- "pcent"
  attr(y, "location") <- loc
  attr(y, "variable") <- param
  attr(y, "longname") <- longname
  attr(y, "unit") <- unit
  attr(y, "quality") <- quality
  attr(y, "calendar") <- "gregorian"
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
  y <- g2dl(y,greenwich=greenwich)
  invisible(y) 
}

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
  i.t <- colnames(x)=="trajectory"
  i.lat <- colnames(x)=="lat"
  i.lon <- colnames(x)=="lon"
  y$trajectory <- unlist(apply(x,1,function(z) rep(z[i.t],z[i.n])))
  y$n <- unlist(apply(x,1,function(z) rep(z[i.n],z[i.n])))
  if(verbose) print('interpolate time')
  i.start <- colnames(x)=="start"
  i.end <- colnames(x)=="end"
  datetime <- unlist(apply( x, 1, function(z) strftime(
        seq(strptime(z[i.start],format="%Y%m%d%H"),
        strptime(z[i.end],format="%Y%m%d%H"),
        length.out=z[i.n]),format="%Y%m%d%H")))
  y$date <- as.numeric(strftime(strptime(datetime,"%Y%m%d%H"),"%Y%m%d"))
  y$time <- as.numeric(strftime(strptime(datetime,"%Y%m%d%H"),"%H"))
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


#map.events <- function(x,it=NULL,is=NULL,
#               projection="sphere",verbose=TRUE,...) {
#  y <- subset(x,it=it,is=is)
#  Y <- as.field(y)
#  map(Y,...)  
#}

subset.events <- function(x,it=NULL,is=NULL,verbose=FALSE,...) {
  if(verbose) print("subset.events")
  if (length(it)==0) it <- NULL
  if (length(is)==0) is <- NULL
  ## date vector
  d <- strptime(paste(x["date"][[1]],x["time"][[1]]),"%Y%m%d %H")
  d <- as.POSIXct(d)
  yr <- year(d)
  mo <- month(d)
  dy <- day(d)
  hr <- as.numeric(format(d,"%H"))

  ii <- rep(TRUE,length(yr))
  if (is.character(it)) {
    if (verbose) print('it is character')
    nlev <- levels(factor(nchar(it)))
    if (length(nlev)==1 & is.element(nlev[1],c(4,8,10,13))) {
      if (nlev==13) {
        it <- as.POSIXct(strptime(it,format="%Y-%m-%d %H")); t <- d
      } else if (nlev==10) {
        if (any(grep("-",it[1]))) {
          it <- as.Date(strptime(it,format="%Y-%m-%d"))
          t <- as.Date(d)
        } else {
          it <- as.POSIXct(strptime(it,format="%Y%m%d%H"))
          t <- d
        }
      } else if (nlev==8) {
        it <- as.Date(strptime(it,format="%Y%m%d"))
        t <- as.Date(d)      
      } else if (nlev==4) {
        t <- yr
      }
      if (verbose) {print(it); print(class(x))}
      if ( length(it) == 2 ) {
        if (verbose) print('Between two dates')
        ii <- (t >= min(it)) & (t <= max(it))
      } else {
        if (verbose) print('Dates:'); print(paste(it,collapse=", "))
        ii <- is.element(t,it)
      }
    } else {
      if (verbose) print('it is a string')
      if (sum(is.element(tolower(substr(it,1,3)),tolower(month.abb)))>0) {
        if (verbose) print('Monthly selected')
        ii <- is.element(mo,(1:12)[is.element(tolower(month.abb),
                                              tolower(substr(it,1,3)))])
      } else if (sum(is.element(tolower(it),names(season.abb())))>0) {
        if (verbose) print("Seasonally selected")
        if (verbose) print(table(mo))
        if (verbose) print(eval(parse(text=paste('season.abb()$',it,sep=''))))
        ii <- is.element(mo,eval(parse(text=paste('season.abb()$',it,sep=''))))
      } else {
        str(it); print(class(it))
        ii <- rep(FALSE,length(d))
        warning("default.subset: did not recognise the selection citerion for 'it'")
      }
    }
  } else if (inherits(it,c("numeric","integer"))) {
    if (verbose) print('it is numeric or integer')
    nlev <- as.numeric(levels(factor(as.character(it)))) # REB 2015-01-15
    if (verbose) {print(nlev); print(it)}
    if (length(it)==2) {
      if ( (min(it) >= 1800) & (max(it) <= 2500) ) {
        if (verbose) print("it most probably contains a years")
        ii <- is.element(yr,year(it[1]):year(it[2]))
      } else if ( (min(it) <= min(yr)) & (max(it) <= length(yr)) ) {
        if (verbose) print("it most probably contains a indices")
        ii <- is.element(1:length(yr),it[1]:it[2])
      } else  if (min(it) >= min(yr)) {
        if (verbose) print("it most probably contains years")
        ii <- is.element(yr,it[1]:max(yr))
      } else  if (max(it) <= max(yr)) {
        if (verbose) print("it most probably contains years")
        ii <- is.element(yr,min(yr):it[2])
      }
    } else if (length(it)>2) {
      if (min(it) > length(yr)) {
        if (verbose) print("match years")
        ii <- is.element(yr,it)
      } else if (max(it) <= length(yr)) {
        if (verbose) print("pick by indices")
        ii <- is.element(1:length(yr),it)
      }
    } else if ((nlev<=4) & (it <=4)) {
      if (verbose) print("it are most probably seasons")
      if (inherits(x,'season') & (length(it)==1)) {
        if (verbose) print(paste("The 'it' value must be a season index",
           "between 1 and 4. If not please use character strings instead,",
           "e.g., it='djf'"))
        it <- switch(it,'1'=1,'2'=4,'3'=7,'4'=10)
        ii <- is.element(mo,it)
    } else if (max(it) <=12) {
      if (verbose) {
        print("The 'it' value must be a month index.")
        print("If not please use character strings instead")
      }
      ii <- is.element(mo,it)
      } else {
        if (verbose) print("it represents indices")
        ii <- it
      }
    } else if (nlev<=12) {
      if (verbose) {
        print("The 'it' value are most probably a month index. ")
        print("If not please use character strings instead")
      }
      ii <- is.element(mo,it)       
    } else {
      if ( (min(it) >= min(yr)) & (max(it) <= max(yr)) ) {
        if (verbose) print("it most probably contains years")
        ii <- is.element(yr,it)
      } else {
        if (verbose)  print("it most probably holds indices")
        ii <- 1:length(d) %in% it
      }
    }
  } else if (inherits(it,c("Date","yearmon","POSIXt"))) {
    if (verbose) print('it is a date object')
    if (inherits(it,"yearmon")) t <- as.yearmon(d) else
    if (inherits(it,"Date")) t <- as.Date(d) else
    if (inherits(it,"POSIXt")) it <- as.POSIXct(it); t <- d
    if (length(it)==2) { ii <- (t >= min(it)) & (t <= max(it)) 
    } else ii <- is.element(t,it)
  } else if (inherits(it,"logical") & length(it)==length(yr)) {
    ii <- it
  } else if (!is.null(it)) {
    ii <- rep(FALSE,length(yr))
    warning("default.subset: did not reckognise the selection citerion for 'it'")
  }

  jj <- rep(TRUE,length(yr))
  if (inherits(is,'list')) {
    if (verbose) print('is is a list:')
    nm.x <- names(x)
    nm.is <- names(is)
    ok <- sapply(nm.is,function(n) any(grep(n,nm.x)))
    if (verbose) print(nm.is[ok])
    for (n in nm.is[ok]) {
      jj <- jj & x[n][[1]]>=min(is[n][[1]]) &
                 x[n][[1]]<=max(is[n][[1]])
    }
  } else if (is.numeric(is)) {
    jj <- 1:length(yr) %in% it
  } else if (is.logical(is) & length(is)==length(yr)) {
    jj <- is  
  } else if (!is.null(is)){
    jj <- rep(FALSE,length(yr))
    warning("default.subset: did not reckognise the selection citerion for 'is'")
  }

  ij <- ii & jj
  y <- x[ij,]
  attr(y,"aspect") <- "subset"
  invisible(y)
}

as.station.events <- function(x,...) {
  y <- param.events(x,...)
  invisible(y)
}

count.events <- function(x,by.trajectory=TRUE,verbose=FALSE,...) {
  if (verbose) print("count.events")
  dates <- as.Date(strptime(paste(x$date,x$time),format="%Y%m%d %H"))
  fn <- function(x) as.Date(as.yearmon(x))
  if (by.trajectory) {
    if (!"trajectory" %in% names(x)) x <- Track.events(x)
    z <- zoo(x$trajectory,order.by=dates)
    N <- aggregate(z,by=fn,FUN=function(x) length(unique(x,na.rm=TRUE)))
  } else {
    z <- zoo(x$date,order.by=dates)
    N <- aggregate(z,by=fn,FUN=length)
  }
  # fill in missing months by merging with an empty time series
  nrt <- as.Date(strptime(range(year(dates))*1E4+range(month(dates))*1E2+1,
                format="%Y%m%d"))
  N0 <- zoo(,seq(from = nrt[1], to = nrt[2], by = "month"))
  N <- merge(N, N0)
  N[is.na(N)] <- 0
  N <- attrcp(x,N)
  N <- as.station(N)
  invisible(N)
}

param.events <- function(x,param="count",FUN="mean",verbose=TRUE,
                         longname=NULL,unit=NULL,...) {
  if (verbose) print("param.events")
  dates <- as.Date(strptime(paste(x$date,x$time),format="%Y%m%d %H"))
  fn <- function(x) as.Date(as.yearmon(x))
  if (param=="count") {
    N <- count.events(x,...)
    longname <- paste(attr(x,"variable"),param)
    unit <- "events/months"
  } else if (param %in% names(x)) {
    y <- zoo(x[,param],order.by=dates)
    N <- aggregate(y,by=fn,FUN=FUN)
    nrt <- as.Date(strptime(range(year(dates))*1E4+range(month(dates))*1E2+1,
              format="%Y%m%d"))
    N0 <- zoo(,seq(from = nrt[1], to = nrt[2], by = "month"))
    N <- merge(N, N0)
    N <- attrcp(x,N)
    N <- as.station(N)
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
  attr(N,"variable") <- param
  attr(N,"longname") <- longname
  attr(N,"unit") <- unit
  attr(N,"lat") <- attr(x,"lat")
  attr(N,"lon") <- attr(x,"lon")
  invisible(N)
}

## Moved events2trajectory to as.trajectory.R
# events2trajectory <- function(x,verbose=FALSE,loc=NA,param=NA,longname=NA,
#                           quality=NA,src=NA,url=NA,reference=NA,info=NA,
#                           method=NA,minlen=5,n=10) {
#   if (verbose) print("events2trajectory")
#   stopifnot(inherits(x,"events"))
#   if (verbose) print(paste('dim: ',paste(dim(x),collapse=" x ")))
#   if (verbose) print(paste('names: ',paste(names(x),collapse=", ")))
# 
#   if (!any("trajectory" %in% names(x))) {
#     x <- Track.events(x,verbose=verbose)
#   }
#   if (!any("trackcount" %in% names(x))) {
#     x <- Trackstats(x,verbose=verbose)
#   }
#   y <- x[order(x$trajectory),]
#   y <- attrcp(x,y)
#   y <- y[y$trackcount>=minlen,]
#   
#   names(y) <- tolower(names(y))
#   names(y)[grep("latitude",names(y))] <- "lat"
#   names(y)[grep("longitude",names(y))] <- "lon"
#   names(y)[grep("step",names(y))] <- "timestep"
# 
#   if(is.na(loc) & !is.null(attr(y,"loc"))) loc <- attr(y,"loc")
#   if(is.na(param) & !is.null(attr(y,"variable"))) param <- attr(y,"variable")
#   if(is.na(longname) & !is.null(attr(y,"longname"))) longname <- attr(y,"longname")
#   if(is.na(quality) & !is.null(attr(y,"quality"))) quality <- attr(y,"quality")
#   if(is.na(src) & !is.null(attr(y,"source"))) src <- attr(y,"source")
#   if(is.na(url) & !is.null(attr(y,"URL"))) url <- attr(y,"URL")
#   if(is.na(reference) & !is.null(attr(y,"reference"))) reference <- attr(y,"refe
# rence")
#   if(is.na(info) & !is.null(attr(x,"info"))) info <- attr(y,"info")
#   if(is.na(method) & !is.null(attr(x,"method"))) method <- attr(y,"method")
#   
#   if(verbose) print('data.frame to numeric values')
#   fn <- function(x) {
#     if(!is.null(levels(x))) {suppressWarnings(as.numeric(levels(x))[x])
#     } else as.numeric(as.character(x))
#   }
#   nlist1 <- c('trajectory','lat','lon','date','time')
#   if (sum(nlist1 %in% names(y))==length(nlist1)) {
#     y$trajectory <- fn(y$trajectory)
#     y$lat <- fn(y$lat)
#     y$lon <- fn(y$lon)
#     y$date <- fn(y$date)
#     y$time <- fn(y$time)
#     y$datetime <- y$date*1E2 + y$time
#   } else {
#     print(paste(paste(nlist1[!(nlist1 %in% names(y))],collapse=' '),'missing'))
#   }
#   nlist2 <- names(y)[!(names(y) %in% nlist1)]
#   for (name in nlist2) {
#     eval(parse(text=paste("y$",name,"<-fn(y$",name,")",sep="")))
#   }
# 
#   if(verbose) print('interpolate trajectories')
#   len <- aggregate(y$datetime, list(y$trajectory), length)$x
#   t1 <- aggregate(y$datetime, list(y$trajectory), function(x) x[1])$x
#   t2 <- aggregate(y$datetime, list(y$trajectory), function(x) x[length(x)])$x
#   a <- 6.378e06
#   xx <- a * cos( y$lat*pi/180 ) * cos( y$lon*pi/180 )
#   yy <- a * cos( y$lat*pi/180 ) * sin( y$lon*pi/180 )
#   zz <- a * sin( y$lat*pi/180 )
#   xa <- aggregate(xx, list(y$trajectory), function(x) approx(x,n=n)$y)$x
#   ya <- aggregate(yy, list(y$trajectory), function(x) approx(x,n=n)$y)$x
#   za <- aggregate(zz, list(y$trajectory), function(x) approx(x,n=n)$y)$x
#   lon <- atan2( ya, xa )*180/pi
#   lat <- asin( za/sqrt( xa^2 + ya^2 + za^2 ))*180/pi
#   colnames(lon) <- rep('lon',n)
#   colnames(lat) <- rep('lat',n)
#   if(verbose) print(n[1:n])
#   nlist1 <- c('trajectory','lat','lon','year','month','day','time',
#              'date','code99','timestep')
#   nlist2 <- names(x)[!(names(x) %in% nlist1)]
#   if(is.null(nlist2)) {
#     X <- cbind(lon=lon,lat=lat,start=t1,end=t2,n=len)
#   } else {
#     for(name in nlist2) {
#       if(verbose) print(name)
#       eval(parse(text=paste("aggregate(x$",name,
#        ",list(x$trajectory),function(x) approx(x,n=",
#                    n,")$y)$x ->",name,sep="")))
#       eval(parse(text=paste("colnames(",name,
#                    ")<-rep('",name,"',",n,")",sep="")))
#     }
#     X <- eval(parse(text=paste("cbind(lon=lon,lat=lat",
#         paste(nlist2,"=",nlist2,sep="",collapse=","),
#         "start=t1,end=t2,n=len)",sep=",")))
#   }
#   attr(X, "location")= loc
#   attr(X, "variable")= param
#   attr(X, "longname")= longname
#   attr(X, "quality")= quality
#   attr(X, "calendar")= "gregorian"
#   attr(X, "source")= src
#   attr(X, "URL")= url
#   attr(X, "type")= "analysis"
#   attr(X, "aspect")= "interpolated"
#   attr(X, "reference")= reference
#   attr(X, "info")= info
#   attr(X, "method")= method
#   attr(x,"lon") <- NA
#   attr(x,"lat") <- NA
#   attr(x,"alt") <- NA
#   attr(x,"cntr") <- NA
#   attr(x,"stid") <- NA
#   attr(X, "history")= history.stamp()
#   class(X) <- c('trajectory','matrix')
#   invisible(X)
# }

