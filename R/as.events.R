
as.events <- function(x,...) UseMethod("as.events")

as.events.default <- function(x,label=NULL,dx=NULL,dy=NULL,
                      units=NULL,longname=NULL,variable=NULL,
                      qflabel=NULL,src=NULL,file=NULL,
                      version=NULL,verbose=FALSE) {
  if (verbose) print("as.events")
  X <- data.frame(x)
  n <- names(X)
  if (!all(c("date","time","lon","lat") %in% names(X))) {
    print(paste("Missing input:",
     names(X)[!c("date","time","lon","lat")%in%names(X)]))
  }
  attr(X,"label") <- label
  attr(X,"dx") <- dx
  attr(X,"dy") <- dy
  attr(X,"longname") <- longname
  attr(X,"variable") <- variable
  attr(X,"QF") <- qflabel
  attr(X,"source") <- src
  attr(X,"file") <- file
  attr(X,"version") <- version
  class(X) <- c("events",class(X))
  attr(X,"history") <- history.stamp(X)
  invisible(X)
}


as.events.trajectory <- function(x,...) {
  X <- trajectory2events(x,...)
  invisible(X)
}

trajectory2events <- function(x,verbose=FALSE,loc=NA,param=NA,longname=NA,
                          quality=NA,src=NA,url=NA,reference=NA,info=NA,
                          method=NA,n=10) {
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
  if(is.na(src) & !is.null(attr(x,"source"))) src <- attr(x,"source")
  if(is.na(url) & !is.null(attr(x,"URL"))) url <- attr(x,"URL")
  if(is.na(reference) & !is.null(attr(x,"reference"))) reference <- attr(x,"reference")
  if(is.na(info) & !is.null(attr(x,"info"))) info <- attr(x,"info")
  if(is.na(method) & !is.null(attr(x,"method"))) method <- attr(x,"method")
  
  # transform data in data.frame x to numeric values
  if(verbose) print('data.frame to numeric values')
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
  if(verbose) print('interpolate trajectories')
  aggregate(x$date, list(x$trajectory), function(x) x[1])$x -> t1
  aggregate(x$date, list(x$trajectory), function(x) x[length(x)])$x -> t2
  aggregate(x$date, list(x$trajectory), length)$x -> len
  aggregate(x$lat, list(x$trajectory),function(x) approx(x,n=n)$y)$x -> lat
  #for(i in seq(11000,12000)) loni<-approxlon(x$lon[x$trajectory==i]);print(i)
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
