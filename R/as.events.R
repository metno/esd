
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
  attr(X,"quality") <- qflabel
  attr(X,"source") <- src
  attr(X,"file") <- file
  attr(X,"version") <- version
  class(X) <- c("events",class(X))
  attr(X,"history") <- history.stamp(X)
  invisible(X)
}

map.events <- function(x,it=NULL,is=NULL,
               projection="sphere",verbose=TRUE,...) {
  y <- subset(x,it=it,is=is)
  Y <- as.field(y)
  map(Y,...)  
}

as.field.events <- function(x,...) {
  y <- events2field(x,...)
  invisible(y)  
}

events2field <- function(x,dt="month",dx=2,dy=2,it=NULL,is=NULL,
                         verbose=FALSE,...) {
  y <- subset(x,it=it,is=is)
  lons <- round(y["lon"]/dx)*dx
  lons <- seq(min(lons),max(lons),dx)
  lats <- round(y["lat"]/dy)*dy
  lats <- seq(min(lats),max(lats),dy)
  #d <- mapply(function(x,y) strptime(paste(x,y),"%Y%m%d %H"),
  #            y["date"],y["time"])
  d <- strptime(y["date"][[1]],"%Y%m%d")
  if (grepl('month',dt)) {
    d <- as.Date(as.yearmon(d))
    dvec <- unique(d)
  } else if (grepl('season',dt) | grepl('quarter',dt)) {
    
  
  dens <- density.events()
  hits <- as.data.frame(table(lon,lat))
  invisible(x)
}

factor2numeric <- function(f) {
  if(!is.null(levels(f))) {return(as.numeric(levels(f))[f])
  } else return(as.numeric(f))
}

density.events <- function(x,it=NULL,is=NULL,dx=NULL,dy=NULL,
                           verbose=FALSE,R=6371,...) {
  y <- subset(x,it=it,is=is)
  lons <- y["lon"][[1]]
  lats <- y["lat"][[1]]
  if(!is.null(dx)) { lons <- round(lons/dx)*dx
  } else dx <- min(diff(sort(unique(lons))))
  if(!is.null(dy)) { lats <- round(lats/dy)*dy
  } else dy <- min(diff(sort(unique(lats))))
  hits <- as.data.frame(table(lons,lats))
  hits <- hits[hits$Freq>0,]
  lons <- factor2numeric(hits$lon)
  lats <- factor2numeric(hits$lat)
  A <- dx*(pi/180)*R**2*abs(sin((lat+dy/2)*pi/180)-
                            sin((lat-dy/2)*pi/180))
  dens <- hits$Freq/A
  X <- data.frame(lon=lon,lat=lat,density=dens)
  invisible(X)
}

subset.events <- function(x,it=NULL,is=NULL,verbose=FASLE,...) {
  if(verbose) print("subset.events")

  ## Sometimes 'it' = 'integer(0)' - reset to NULL!
  if (length(it)==0) it <- NULL
  if (length(is)==0) is <- NULL
  ## date vector
  d <- strptime(paste(x["date"][[1]],x["time"][[1]]),"%Y%m%d %H")
  yr <- year(d)
  mo <- month(d)
  dy <- day(d)
  hr <- as.numeric(format(d,"%H"))
  
  if (is.character(it)) {
    if (verbose) print('it is character')
    if ((levels(factor(nchar(it)))==10)) ##  convert as Date
      it <- strptime(it,"%Y-%m-%d")
      #it <- as.Date(it)
    if ( length(it) == 2 ) {
      if (verbose) print('Between two dates')
      if(length(unique(hr))>1)
        it <- seq(it[1],it[2],by=) else
      if (nchar(it[1])==4) it <- as.Date(c(paste(it[1],'-01-01',sep=''), 
                                         paste(it[2],'-12-31',sep='')))
      if (verbose) {print(it); print(class(x))}
      if (inherits(x,"month"))
        it <- seq(it[1],it[2],by='month') else
      if (inherits(x,"season"))
        it <- seq(it[1],it[2],by='season') else
      if (inherits(x,"annual"))
        it <- seq(it[1],it[2],by='year') else 
      if (inherits(x,"day"))
        it <- seq(it[1],it[2],by='day') else 
        
      ii <- is.element(t,it)
      if (verbose) print(paste('sum(ii)=',sum(ii)))
        } else { ## added AM 10-06-2014
            if (verbose) print('it is a string')
            if (sum(is.element(tolower(substr(it,1,3)),tolower(month.abb)))>0) {
                if (verbose) print('Monthly selected')
                ii <- is.element(month(x),(1:12)[is.element(tolower(month.abb),tolower(substr(it,1,3)))])
                                        #y <- x[ii,is] #  REB Not here
            } else if (sum(is.element(tolower(it),names(season.abb())))>0) {
                if (verbose) print("Seasonally selected")
                if (verbose) print(table(month(x)))
                if (verbose) print(eval(parse(text=paste('season.abb()$',it,sep=''))))
                ii <- is.element(month(x),eval(parse(text=paste('season.abb()$',it,sep=''))))
                                        #y <- x[ii,is] # REB Not here
            }
            else if (inherits(it,"Date")) {
                if (verbose) print('it is a Date object')
                ii <- is.element(t,it)
            } else {
                str(it); print(class(it))
                ii <- rep(FALSE,length(t))
                warning("default.subset: did not recognise the selection citerion for 'it'")
            }
        }
        ##    } else if ((class(it)=="numeric") | (class(it)=="integer")) {
##        if (min(it) > 1500) ## it is a year
##            it <- seq(it[1],it[2],by=1)
##        ##print("HERE"); print(it)
##    }
##    ii <- (t >= it[1]) & (t <= it[length(it)])
##
    ## get the subset indices in ii
    } else if ((class(it)=="numeric") | (class(it)=="integer")) {
        if (verbose) print('it is numeric or integer')
# REB bug        nlev <- as.numeric(levels(factor(nchar(it))))
# nchar returns the string length, but these lines need to find the number of different levels/categories
        nlev <- as.numeric(levels(factor(as.character(it)))) # REB 2015-01-15
        if (verbose) {print(nlev); print(it)}
#         if ((length(nlev)==1)) { REB 2015-01-20: the lines below will never happen with this line:
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
          # if it is years:
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
                    if (verbose)  print(paste("The 'it' value must be a season index between 1 and 4.",
                                              "If not please use character strings instead. e.g. it='djf'"))
                    it <- switch(it,'1'=1,'2'=4,'3'=7,'4'=10)
                    ii <- is.element(mo,it)
                 } else if ( (max(it) <=12) &
                             (inherits(x,'month') | (inherits(x,'day')))) {
                     if (verbose) {
                         print("The 'it' value must be a month index.")
                         print("If not please use character strings instead")
                       }
                     ii <- is.element(mo,it)
                 }  else {
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
            #  length(nlev) > 1
                if ( (min(it) >= min(yr)) & (max(it) <= max(yr)) ) {
                  if (verbose) print("it most probably contains years")
                  ii <- is.element(yr,it)
                } else {
                  if (verbose)  print("it most probably holds indices")
                  ii <- it
                }
              }
      } else if (inherits(it,c("Date","yearmon"))) {       
        ##        ii <- is.element(t,it)
        if (verbose) print('it is a date object')
        ii <- (t >= min(it)) & (t <= max(it))
      } else if (inherits(it,"logical") & length(it)==length(yr)) {
        ii <- it
      } else if (!is.null(it)) {
        ii <- rep(FALSE,length(t))
        warning("default.subset: did not reckognise the selection citerion for 'it'")
      } 
  
  invisible(x)
}
 

as.events.trajectory <- function(x,...) {
  X <- trajectory2events(x,...)
  invisible(X)
}

trajectory2events <- function(x,verbose=FALSE) {
  stopifnot(inherits(x,"trajectory"))
  
  ## if (verbose) print("trajectory")
  ## if (verbose) print(paste('dim: ',paste(dim(x),collapse=" x ")))
  ## if (verbose) print(paste('names: ',paste(names(x),collapse=", ")))
  ## names(x) <- tolower(names(x))
  ## names(x)[grep("latitude",names(x))] <- "lat"
  ## names(x)[grep("longitude",names(x))] <- "lon"
  ## names(x)[grep("step",names(x))] <- "timestep"

  ## if(is.na(loc) & !is.null(attr(x,"loc"))) loc <- attr(x,"loc")
  ## if(is.na(param) & !is.null(attr(x,"variable"))) param <- attr(x,"variable")
  ## if(is.na(longname) & !is.null(attr(x,"longname"))) longname <- attr(x,"longname")
  ## if(is.na(quality) & !is.null(attr(x,"quality"))) quality <- attr(x,"quality")
  ## if(is.na(src) & !is.null(attr(x,"source"))) src <- attr(x,"source")
  ## if(is.na(url) & !is.null(attr(x,"URL"))) url <- attr(x,"URL")
  ## if(is.na(reference) & !is.null(attr(x,"reference"))) reference <- attr(x,"reference")
  ## if(is.na(info) & !is.null(attr(x,"info"))) info <- attr(x,"info")
  ## if(is.na(method) & !is.null(attr(x,"method"))) method <- attr(x,"method")
  
  ## # transform data in data.frame x to numeric values
  ## if(verbose) print('data.frame to numeric values')
  ## x <- data.frame(x)
  ## fn <- function(x) {
  ##   if(!is.null(levels(x))) {suppressWarnings(as.numeric(levels(x))[x])
  ##   } else as.numeric(as.character(x))
  ## }
  ## nlist1 <- c('trajectory','lat','lon','year','month','day','time')
  ## if (sum(nlist1 %in% names(x))==length(nlist1)) {
  ##   x$trajectory <- fn(x$trajectory)
  ##   x$lat <- fn(x$lat)
  ##   x$lon <- fn(x$lon)
  ##   x$year <- fn(x$year)
  ##   x$month <- fn(x$month)
  ##   x$day <- fn(x$day)
  ##   x$time <- fn(x$time)
  ##   x$date <- x$year*1E6+x$month*1E4+x$day*1E2+x$time
  ## } else {
  ##   print(paste(paste(nlist1[!(nlist1 %in% names(x))],collapse=' '),'missing'))
  ## }
  ## nlist2 <- names(x)[!(names(x) %in% c('code99','date',nlist1))]
  ## for (name in nlist2) {
  ##   eval(parse(text=paste("x$",name,"<-fn(x$",name,")",sep="")))
  ## }

  ## # interpolate all trajectories to same length n
  ## if(verbose) print('interpolate trajectories')
  ## aggregate(x$date, list(x$trajectory), function(x) x[1])$x -> t1
  ## aggregate(x$date, list(x$trajectory), function(x) x[length(x)])$x -> t2
  ## aggregate(x$date, list(x$trajectory), length)$x -> len
  ## aggregate(x$lat, list(x$trajectory),function(x) approx(x,n=n)$y)$x -> lat
  ## #for(i in seq(11000,12000)) loni<-approxlon(x$lon[x$trajectory==i]);print(i)
  ## aggregate(x$lon,list(x$trajectory),function(x) approxlon(x,n=n)$y)$x -> lon
  ## colnames(lon) <- rep('lon',n)
  ## colnames(lat) <- rep('lat',n)
  ## if(verbose) print(n[1:n])
  ## nlist1 <- c('trajectory','lat','lon','year','month','day','time',
  ##            'date','code99','timestep')
  ## nlist2 <- names(x)[!(names(x) %in% nlist1)]
  ## if(is.null(nlist2)) {
  ##   X <- cbind(lon=lon,lat=lat,start=t1,end=t2,n=len)
  ## } else {
  ##   for(name in nlist2) {
  ##     if(verbose) print(name)
  ##     eval(parse(text=paste("aggregate(x$",name,
  ##      ",list(x$trajectory),function(x) approx(x,n=",
  ##                  n,")$y)$x ->",name,sep="")))
  ##     eval(parse(text=paste("colnames(",name,
  ##                  ")<-rep('",name,"',",n,")",sep="")))
  ##   }
  ##   X <- eval(parse(text=paste("cbind(lon=lon,lat=lat",
  ##       paste(nlist2,"=",nlist2,sep="",collapse=","),
  ##       "start=t1,end=t2,n=len)",sep=",")))
  ## }

  ## # add attributes to trajectory matrix X
  ## attr(X, "location")= loc
  ## attr(X, "variable")= param
  ## attr(X, "longname")= longname
  ## attr(X, "quality")= quality
  ## attr(X, "calendar")= "gregorian"
  ## attr(X, "source")= src
  ## attr(X, "URL")= url
  ## attr(X, "type")= "analysis"
  ## attr(X, "aspect")= "interpolated"
  ## attr(X, "reference")= reference
  ## attr(X, "info")= info
  ## attr(X, "method")= method
  ## attr(x,"lon") <- NA
  ## attr(x,"lat") <- NA
  ## attr(x,"alt") <- NA
  ## attr(x,"cntr") <- NA
  ## attr(x,"stid") <- NA
  ## attr(X, "history")= history.stamp()
  ## class(X) <- 'trajectory'
  ## invisible(X)
}
