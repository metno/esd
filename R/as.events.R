
map.events <- function(x,it=NULL,is=NULL,
               projection="sphere",verbose=TRUE,...) {
  y <- subset(x,it=it,is=is)
  Y <- as.field(y)
  map(Y,...)  
}

as.field.events <- function(x,...) {
  y <- density.events(x,...)
  invisible(y)  
}

density.events <- function(x,dt="month",dx=2,dy=2,it=NULL,is=NULL,
                         verbose=FALSE,...) {
  if (verbose) print("density.events")
  y <- subset.events(x,it=it,is=is)
  if(is.null(dx)) dx <- min(diff(sort(unique(y["lon"][[1]]))))
  if(is.null(dy)) dy <- min(diff(sort(unique(y["lat"][[1]]))))
  lons <- round(x["lon"]/dx)*dx
  lons <- seq(min(lons),max(lons),dx)
  lats <- round(x["lat"]/dy)*dy
  lats <- seq(min(lats),max(lats),dy)
  if (grepl('day',dt)) {
    d <- as.Date(strptime(y["date"][[1]],"%Y%m%d"))
    dvec <- seq(min(d),max(d),by="day")
    unit <- 'events/day/km^2'
  } else if (grepl('month',dt)) {
    d <- as.Date(as.yearmon(strptime(y["date"][[1]],"%Y%m%d")))
    dvec <- seq(min(d),max(d),by="month")
    unit <- 'events/month/km^2'
  } else if (grepl('season',dt) | grepl('quarter',dt)) {
    d <- as.Date(as.yearqtr(strptime(y["date"][[1]],"%Y%m%d")))
    dvec <- seq(min(d),max(d),by="quarter")
    unit <- 'events/quarter/km^2'
  } else if (grepl('year',dt) | grepl('annual',dt)) {
    d <- as.Date(paste(year(strptime(y["date"][[1]],"%Y%m%d")),"-01-01",sep=""))
    dvec <- seq(min(d),max(d),by="year")
    unit <- 'events/year/km^2'
  } else if (grepl('hour',dt)) {
    d <- as.POSIXct(strptime(paste(y["date"][[1]],y["time"][[1]]),format="%Y%m%d %H"))
    dh <- min(diff(sort(unique(y["time"][[1]]))))
    dvec <- seq(min(d),max(d),by=dh*60*60)
    unit <- paste('events/',dh,'hours/km^2',sep='')
  } else {
    print(paste("WARNING! invalid time resolution dt",dt))
    break
  }
  X <- array(rep(0,),dim=c(length(dvec),length(lons),length(lats)))
  if (verbose) print("looping...")
  for (i in dvec) {
    if (verbose) print(as.Date(i))
    yi <- subset.events(y,it=(d==i))
    if (!is.null(dim(yi)) & length(yi)>0) {
      print("events")
      di <- lonlatdensity(yi["lon"][[1]],yi["lat"][[1]],dx=dx,dy=dy,verbose=verbose)
      for(j in 1:length(di$lat)) {
        X[dvec==i,lons==di$lon[j],lats==di$lat[j]] <- di$density[j]
      }
    }
  }
  d <- dim(X)
  dim(X) <- c(d[1],d[2]*d[3])
  longname <- paste("event density",attr(x,'longname'),sep=', ')
  param <- 'density'
  Y <- as.field(X,index=dvec,lon=lons,lat=lats,
          unit=unit,longname=longname,param=param,
          quality=attr(x,'quality'),src=attr(x,'source'),
          url=attr(x,'URL'),reference=attr(x,'reference'),
          info=attr(x,'info'),calendar=attr(x,'calendar'),
          method=attr(x,'method'),aspect=attr(x,'aspect'))
  invisible(Y)
}

factor2numeric <- function(f) {
  if(!is.null(levels(f))) {return(as.numeric(levels(f))[f])
  } else return(as.numeric(f))
}

lonlatdensity <- function(lons,lats,dx=NULL,dy=NULL,verbose=FALSE,R=6371,...) {
   if (verbose) print("lonlatdensity")
   lons <- round(lons/dx)*dx
   lats <- round(lats/dy)*dy
   hits <- as.data.frame(table(lons,lats))
   hits <- hits[hits$Freq>0,]
   lons <- factor2numeric(hits$lon)
   lats <- factor2numeric(hits$lat)
   A <- dx*(pi/180)*R**2*abs(sin((lats+dy/2)*pi/180)-
                            sin((lats-dy/2)*pi/180))
   d <- hits$Freq/A
   X <- data.frame(lon=lons,lat=lats,density=d)
   invisible(X)
}

subset.events <- function(x,it=NULL,is=NULL,verbose=FALSE,...) {
  if(verbose) print("subset.events")

  ## Sometimes 'it' = 'integer(0)' - reset to NULL!
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
 

as.events.trajectory <- function(x,...) {
  X <- trajectory2events(x,...)
  invisible(X)
}


events2trajectory <- function(x,verbose=FALSE,loc=NA,param=NA,longname=NA,
                          quality=NA,src=NA,url=NA,reference=NA,info=NA,
                          method=NA,minlen=5,n=10) {
  if (verbose) print("events2trajectory")
  stopifnot(inherits(x,"events"))
  if (verbose) print(paste('dim: ',paste(dim(x),collapse=" x ")))
  if (verbose) print(paste('names: ',paste(names(x),collapse=", ")))

  if (!any("trajectory" %in% names(x))) {
    x <- Track.events(x,verbose=verbose)
  }
  if (!any("trackcount" %in% names(x))) {
    x <- Trackstats(x,verbose=verbose)
  }
  y <- x[order(x$trajectory),]
  y <- attrcp(x,y)
  y <- y[y$trackcount>=minlen,]
  
  names(y) <- tolower(names(y))
  names(y)[grep("latitude",names(y))] <- "lat"
  names(y)[grep("longitude",names(y))] <- "lon"
  names(y)[grep("step",names(y))] <- "timestep"

  if(is.na(loc) & !is.null(attr(y,"loc"))) loc <- attr(y,"loc")
  if(is.na(param) & !is.null(attr(y,"variable"))) param <- attr(y,"variable")
  if(is.na(longname) & !is.null(attr(y,"longname"))) longname <- attr(y,"longname")
  if(is.na(quality) & !is.null(attr(y,"quality"))) quality <- attr(y,"quality")
  if(is.na(src) & !is.null(attr(y,"source"))) src <- attr(y,"source")
  if(is.na(url) & !is.null(attr(y,"URL"))) url <- attr(y,"URL")
  if(is.na(reference) & !is.null(attr(y,"reference"))) reference <- attr(y,"refe
rence")
  if(is.na(info) & !is.null(attr(x,"info"))) info <- attr(y,"info")
  if(is.na(method) & !is.null(attr(x,"method"))) method <- attr(y,"method")
  
  if(verbose) print('data.frame to numeric values')
  fn <- function(x) {
    if(!is.null(levels(x))) {suppressWarnings(as.numeric(levels(x))[x])
    } else as.numeric(as.character(x))
  }
  nlist1 <- c('trajectory','lat','lon','date','time')
  if (sum(nlist1 %in% names(y))==length(nlist1)) {
    y$trajectory <- fn(y$trajectory)
    y$lat <- fn(y$lat)
    y$lon <- fn(y$lon)
    y$date <- fn(y$date)
    y$time <- fn(y$time)
    y$datetime <- y$date*1E2 + y$time
  } else {
    print(paste(paste(nlist1[!(nlist1 %in% names(y))],collapse=' '),'missing'))
  }
  nlist2 <- names(y)[!(names(y) %in% nlist1)]
  for (name in nlist2) {
    eval(parse(text=paste("y$",name,"<-fn(y$",name,")",sep="")))
  }

  # interpolate all trajectories to same length n
  if(verbose) print('interpolate trajectories')
  len <- aggregate(y$datetime, list(y$trajectory), length)$x
  t1 <- aggregate(y$datetime, list(y$trajectory), function(x) x[1])$x
  t2 <- aggregate(y$datetime, list(y$trajectory), function(x) x[length(x)])$x
  a <- 6.378e06
  xx <- a * cos( y$lat*pi/180 ) * cos( y$lon*pi/180 )
  yy <- a * cos( y$lat*pi/180 ) * sin( y$lon*pi/180 )
  zz <- a * sin( y$lat*pi/180 )
  xa <- aggregate(xx, list(y$trajectory), function(x) approx(x,n=n)$y)$x
  ya <- aggregate(yy, list(y$trajectory), function(x) approx(x,n=n)$y)$x
  za <- aggregate(zz, list(y$trajectory), function(x) approx(x,n=n)$y)$x
  lon <- atan2( ya, xa )*180/pi
  lat <- asin( za/sqrt( xa^2 + ya^2 + za^2 ))*180/pi
  #aggregate(x$lat, list(x$trajectory),function(x) approx(x,n=n)$y)$x -> lat
  #aggregate(x$lon,list(x$trajectory),function(x) approxlon(x,n=n)$y)$x -> lon
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
  class(X) <- c('trajectory','matrix')
  invisible(X)
}


as.station.events <- function(x,...) {
  y <- param.events(x,...)
  invisible(y)
}

count.events <- function(x,by.trajectory=TRUE,verbose=TRUE,...) {
  if (verbose) print("count.events")
  dates <- as.Date(strptime(paste(x$date,x$time),format="%Y%m%d %H"))
  fn <- function(x) as.Date(as.yearmon(x))
  if (by.trajectory) {
    if (!"trajectory" %in% names(x)) x <- Track.events(x)
    y <- zoo(x$trajectory,order.by=dates)
    N <- aggregate(y,by=fn,FUN=function(x) length(unique(x,na.rm=TRUE)))
  } else {
    y <- zoo(x$date,order.by=dates)
    N <- aggregate(y,by=fn,FUN=length)
  }
  # fill in missing months by merging with an empty series containing
  # a complete set of 1st of the months
  nrt <- as.Date(strptime(range(year(dates))*1E4+range(month(dates))*1E2+1,
                format="%Y%m%d"))
  N0 <- zoo(,seq(from = nrt[1], to = nrt[2], by = "month"))
  N <- merge(N, N0)
  N[is.na(N)] <- 0
  N <- attrcp(x,N)
  N <- as.station(N)
  invisible(N)
}

param.events <- function(x,param="count",FUN="mean",verbose=TRUE,...) {
  if (verbose) print("param.events")
  dates <- as.Date(strptime(paste(x$date,x$time),format="%Y%m%d %H"))
  fn <- function(x) as.Date(as.yearmon(x))
  if (param=="count") {
    N <- count.events(x,...)
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
  attr(N,"lat") <- attr(x,"lat")
  attr(N,"lon") <- attr(x,"lon")
  invisible(N)
}

## trajectory2events <- function(x,verbose=FALSE) {
##   stopifnot(inherits(x,"trajectory"))
  
##   if (verbose) print("trajectory")
##   if (verbose) print(paste('dim: ',paste(dim(x),collapse=" x ")))
##   if (verbose) print(paste('names: ',paste(names(x),collapse=", ")))
##   names(x) <- tolower(names(x))
##   names(x)[grep("latitude",names(x))] <- "lat"
##   names(x)[grep("longitude",names(x))] <- "lon"
##   names(x)[grep("step",names(x))] <- "timestep"

##   if(is.na(loc) & !is.null(attr(x,"loc"))) loc <- attr(x,"loc")
##   if(is.na(param) & !is.null(attr(x,"variable"))) param <- attr(x,"variable")
##   if(is.na(longname) & !is.null(attr(x,"longname"))) longname <- attr(x,"longname")
##   if(is.na(quality) & !is.null(attr(x,"quality"))) quality <- attr(x,"quality")
##   if(is.na(src) & !is.null(attr(x,"source"))) src <- attr(x,"source")
##   if(is.na(url) & !is.null(attr(x,"URL"))) url <- attr(x,"URL")
##   if(is.na(reference) & !is.null(attr(x,"reference"))) reference <- attr(x,"reference")
##   if(is.na(info) & !is.null(attr(x,"info"))) info <- attr(x,"info")
##   if(is.na(method) & !is.null(attr(x,"method"))) method <- attr(x,"method")
  
##   # transform data in data.frame x to numeric values
##   if(verbose) print('data.frame to numeric values')
##   x <- data.frame(x)
##   fn <- function(x) {
##     if(!is.null(levels(x))) {suppressWarnings(as.numeric(levels(x))[x])
##     } else as.numeric(as.character(x))
##   }
##   nlist1 <- c('trajectory','lat','lon','year','month','day','time')
##   if (sum(nlist1 %in% names(x))==length(nlist1)) {
##     x$trajectory <- fn(x$trajectory)
##     x$lat <- fn(x$lat)
##     x$lon <- fn(x$lon)
##     x$year <- fn(x$year)
##     x$month <- fn(x$month)
##     x$day <- fn(x$day)
##     x$time <- fn(x$time)
##     x$date <- x$year*1E6+x$month*1E4+x$day*1E2+x$time
##   } else {
##     print(paste(paste(nlist1[!(nlist1 %in% names(x))],collapse=' '),'missing'))
##   }
##   nlist2 <- names(x)[!(names(x) %in% c('code99','date',nlist1))]
##   for (name in nlist2) {
##     eval(parse(text=paste("x$",name,"<-fn(x$",name,")",sep="")))
##   }

##   # interpolate all trajectories to same length n
##   if(verbose) print('interpolate trajectories')
##   aggregate(x$date, list(x$trajectory), function(x) x[1])$x -> t1
##   aggregate(x$date, list(x$trajectory), function(x) x[length(x)])$x -> t2
##   aggregate(x$date, list(x$trajectory), length)$x -> len
##   aggregate(x$lat, list(x$trajectory),function(x) approx(x,n=n)$y)$x -> lat
##   #for(i in seq(11000,12000)) loni<-approxlon(x$lon[x$trajectory==i]);print(i)
##   aggregate(x$lon,list(x$trajectory),function(x) approxlon(x,n=n)$y)$x -> lon
##   colnames(lon) <- rep('lon',n)
##   colnames(lat) <- rep('lat',n)
##   if(verbose) print(n[1:n])
##   nlist1 <- c('trajectory','lat','lon','year','month','day','time',
##              'date','code99','timestep')
##   nlist2 <- names(x)[!(names(x) %in% nlist1)]
##   if(is.null(nlist2)) {
##     X <- cbind(lon=lon,lat=lat,start=t1,end=t2,n=len)
##   } else {
##     for(name in nlist2) {
##       if(verbose) print(name)
##       eval(parse(text=paste("aggregate(x$",name,
##        ",list(x$trajectory),function(x) approx(x,n=",
##                    n,")$y)$x ->",name,sep="")))
##       eval(parse(text=paste("colnames(",name,
##                    ")<-rep('",name,"',",n,")",sep="")))
##     }
##     X <- eval(parse(text=paste("cbind(lon=lon,lat=lat",
##         paste(nlist2,"=",nlist2,sep="",collapse=","),
##         "start=t1,end=t2,n=len)",sep=",")))
##   }

##   # add attributes to trajectory matrix X
##   attr(X, "location")= loc
##   attr(X, "variable")= param
##   attr(X, "longname")= longname
##   attr(X, "quality")= quality
##   attr(X, "calendar")= "gregorian"
##   attr(X, "source")= src
##   attr(X, "URL")= url
##   attr(X, "type")= "analysis"
##   attr(X, "aspect")= "interpolated"
##   attr(X, "reference")= reference
##   attr(X, "info")= info
##   attr(X, "method")= method
##   attr(x,"lon") <- NA
##   attr(x,"lat") <- NA
##   attr(x,"alt") <- NA
##   attr(x,"cntr") <- NA
##   attr(x,"stid") <- NA
##   attr(X, "history")= history.stamp()
##   class(X) <- 'trajectory'
##   invisible(X)
## }
