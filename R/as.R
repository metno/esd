
as.station <- function(x,...) UseMethod("as.station")

as.station.zoo <- function(x,loc=NA,param=NA,unit=NA,lon=NA,lat=NA,alt=NA,
                          cntr=NA,longname=NA,stid=NA,quality=NA,src=NA,
                          url=NA,reference=NA,info=NA, method= NA,type=NA,
                           aspect=NA) {
  #print(c(length(X),length(index)))
  y <- zoo(x,order.by=index(x))

  if (is.null(loc)) loc <- NA
  if ((is.na(loc)) & !is.null(attr(x,'location')))
    loc <- attr(x,'location')
  attr(y,'location') <- loc
  if (is.null(param)) param <- NA
  if ((is.na(param)) & !is.null(attr(x,'variable')))
    param <- attr(x,'variable')
  attr(y,'variable') <- param
  if (is.null(unit)) unit <- NA
  if ((is.na(unit)) & !is.null(attr(x,'unit')))
    unit <- attr(x,'unit')
  attr(y,'unit') <- unit
  if (is.null(lon)) lon <- NA
  if ((is.na(lon[1])) & !is.null(attr(x,'longitude')))
    lon <- attr(x,'longitude')  
  attr(y,'longitude') <- lon
  if (is.null(lat)) lat <- NA
  if ((is.na(lat[1])) & !is.null(attr(x,'latitude')))
    lat <- attr(x,'latitude')    
  attr(y,'latitude') <- lat
  if (is.null(alt)) alt <- NA
  if ((is.na(alt)) & !is.null(attr(x,'altitude')))
    alt <- attr(x,'altitude')    
  attr(y,'altitude') <- alt
  if (is.null(cntr)) cntr <- NA
  if ((is.na(cntr)) & !is.null(attr(x,'country')))
    cntr <- attr(x,'country')      
  attr(y,'country') <- cntr
  if (is.null(longname)) longname <- NA
  if ((is.na(longname)) & !is.null(attr(x,'longname')))
    longname <- attr(x,'longname')    
  attr(y,'longname') <- longname
  if (is.null(stid)) stid <- NA
  if ((is.na(stid)) & !is.null(attr(x,'station_id')))
    stid <- attr(x,'station_id')    
  attr(y,'station_id') <- stid
  if (is.null(quality)) quality <- NA
  if ((is.na(quality[1])) & !is.null(attr(x,'quality')))
    quality <- attr(x,'guality')    
  attr(y,'quality') <- quality
  attr(y,'calendar') <- 'gregorian'
  if (is.null(src)) src <- NA 
  if ((is.na(src)) & !is.null(attr(x,'source')))
    src <- attr(x,'source')    
  attr(y,'source') <- src
  if (is.null(url)) url <- NA 
  if ((is.na(url)) & !is.null(attr(x,'URL')))
    url <- attr(x,'URL')    
  attr(y,'URL') <- url
  #attr(y,'history') <- 'as.station.data.frame'
  #attr(y,'date-stamp') <- date()
  if (is.null(type)) type <- NA 
  if ((is.na(type)) & !is.null(attr(x,'type')))
    type <- attr(x,'type')    
  attr(y,'type') <- type
  if (is.null(aspect)) aspect <- NA 
  if ((is.na(aspect)) & !is.null(attr(x,'aspect')))
    aspect <- attr(x,'aspect')    
  attr(y,'aspect') <- aspect
  if (is.null(reference)) reference <- NA 
  if ((is.na(reference)) & !is.null(attr(x,'reference')))
    reference <- attr(x,'reference')    
  attr(y,'reference') <- reference
  if (is.null(info)) info <- NA 
  if ((is.na(info)) & !is.null(attr(x,'info')))
    info <- attr(x,'info')    
  attr(y,'info') <- info
  if (is.null(method)) method <- NA
  if ((is.na(method)) & !is.null(attr(x,'method')))
    method <- attr(x,'method')    
  attr(y,'method') <- method
  #attr(y,'call') <- match.call()
  attr(y,'history') <- history.stamp(x)
  dt <- as.numeric(levels(factor(diff(index(yy)))))
  if (dt==1) tscale <- 'day' else
  if ( ((dt>=28) & (dt <=31)) |
       (dt < 0.1) ) tscale <- 'month' else
  if ( (dt>=89) & (dt <=93) ) tscale <- 'season' else 
  if ( (dt>=360) & (dt <=366) ) tscale <- 'annual' else
                                tscale <- 'annual'
  class(y) <- c("station",tscale,"zoo")
  return(y)
}



as.station.data.frame <-  function (x, loc = NA, param = NA, unit = NA,
                                    lon = NA, lat = NA,
    alt = NA, cntr = NA, longname = NA, stid = NA, quality = NA,
    src = NA, url = NA, reference = NA, info = NA, method= NA,type=NA,aspect=NA)
{
    cnam <- names(x)
    year <- sort(rep(x[[1]], 12))
    month <- rep(1:12, length(x[[1]]))
    index <- as.Date(paste(year, month, 1, sep = "-"))
    X <- c(t(as.matrix(x)[, 2:13]))
    y <- zoo(X, order.by = index)
    attr(y, "location") <- loc
    attr(y, "variable") <- param
    attr(y, "unit") <- unit
    attr(y, "longitude") <- lon
    attr(y, "latitude") <- lat
    attr(y, "altitude") <- alt
    attr(y, "country") <- cntr
    attr(y, "longname") <- longname
    attr(y, "station_id") <- stid
    attr(y, "quality") <- quality
    attr(y, "calendar") <- "gregorian"
    attr(y, "source") <- src
    attr(y, "URL") <- url
    #attr(y, "history") <- "as.station.data.frame"
    #attr(y, "date-stamp") <- date()
    attr(y, "type") <- type
    attr(y, "aspect") <- aspect
    attr(y, "reference") <- reference
    attr(y, "info") <- info
    attr(y, "method") <- method
    #attr(y, "call") <- match.call()
    attr(y,'history') <- history.stamp(x)
    class(y) <- c("station", "month", "zoo")
    return(y)
}

as.station.ds <- function(x) {
  if (inherits(x,'pca')) {
    class(x) <- class(x)[-1]
    y <- as.station.pca(x)
  } else {
    y <- zoo(coredata(x),order.by=index(x))
    if (!is.null(attr(x,"original_data")))
      y <- attrcp(attr(x,"original_data"),y) else
      y <- attrcp(x,y)
    if (!is.null(attr(x,"original_data")))
      class(y) <- class(attr(x,"original_data")) else
      class(y) <-class(x)
  }
  attr(y,'history') <- history.stamp(x)
  attr(y,'method') <- attr(x,'method')
  attr(y,'info') <- attr(x,'info')
  return(y)
}

as.station.pca <- function(x) {
  y <- pca2station(x)
  return(y)
}



as.station.list <- function(x) {
  #print("as.station.ds")
#  Jan <- x$Jan + attr(x$Jan,'mean')
#  Feb <- x$Feb + attr(x$Feb,'mean')
#  Mar <- x$Mar + attr(x$Mar,'mean')
#  Apr <- x$Apr + attr(x$Apr,'mean')
#  May <- x$May + attr(x$May,'mean')
#  Jun <- x$Jun + attr(x$Jun,'mean')
#  Jul <- x$Jul + attr(x$Jul,'mean')
#  Aug <- x$Aug + attr(x$Aug,'mean')
#  Sep <- x$Sep + attr(x$Sep,'mean')
#  Oct <- x$Oct + attr(x$Oct,'mean')
#  Nov <- x$Nov + attr(x$Nov,'mean')
#  Dec <- x$Dec + attr(x$Dec,'mean')
  cline <- "merge.zoo("
  if (is.list(x)) {
    for (i in 1:length(x)) {
      ave <- switch(attr(x[[i]],'aspect'),
                    'original'=0,
                    'anomaly'=attr(x[[i]],'mean'))
      attr(x[[i]],'type') <- 'downscaled results'
      z <- x[[i]] + ave
      eval(parse(text=paste("ds.",i," <- z",sep="")))
      cline <- paste(cline,"ds.",i,",",sep="")
    }
 
    cline <- paste(substr(cline,1,nchar(cline)-1),')-> ALL')
    #print(cline)
    eval(parse(text=cline))
    y <- zoo(rowMeans(coredata(ALL),na.rm=TRUE),order.by=index(ALL))
    #print(names(y))
  } else {
    if (inherits(x,'ds')) {
      ave <- switch(attr(x,'aspect'),
                    'original'=0,
                    'anomaly'=attr(x,'mean'))
      attr(x,'type') <- 'downscaled results'
      z <- x + ave
      class(z) <- "zoo"
      y <- zoo(coredata(z),order.by=index(z))
      x <- list(x)
    }
  }
  
                     
#  merge.zoo(Jan,Feb,Mar,Apr,May,Jun,
#            Jul,Aug,Sep,Oct,Nov,Dec)-> ALL
  
  attr(y,'location') <- attr(x[[1]],'location')
  attr(y,'variable') <- attr(x[[1]],'variable')
  attr(y,'unit') <- attr(x[[1]],'unit')
  attr(y,'longitude') <- attr(x[[1]],'longitude')
  attr(y,'latitude') <- attr(x[[1]],'latitude') 
  attr(y,'altitude') <- attr(x[[1]],'altitude')
  attr(y,'country') <- attr(x[[1]],'country')
  attr(y,'longname') <- attr(x[[1]],'longname')
  attr(y,'station_id') <- attr(x[[1]],'station_id')
  attr(y,'quality') <- attr(x[[1]],'quality')
  attr(y,'calendar') <- attr(x[[1]],'calendar')
  attr(y,'source') <- attr(x[[1]],'source')
  attr(y,'URL') <- attr(x[[1]],'URL')
  #attr(y,'history') <- c('as.station.ds',attr(x$Jan,'history'))
  #attr(y,'date-stamp') <- attr(x,'date-stamp')
  attr(y,'type') <-  attr(x,'type')
  attr(y,'aspect') <- attr(x,'aspect')
  attr(y,'reference') <- attr(x,'reference')
  attr(y,'info') <- attr(x,'info')
  #attr(y,'date') <- date()
  #attr(y,'call') <- match.call()
  attr(y,'history') <- history.stamp(x)
  class(y) <- c("station",class(x[[1]]))
  return(y)
}

as.station.field <- function(x,is=NULL) {
  stopifnot(!missing(x),
           !zeros(inherits(x,c("field","zoo"),which=TRUE) ))
  if (!is.null(is)) y <- regrid.field(x,list(lon,lat)) else
                    y <- x
  y <- as.station.zoo(y,loc=NA,param=varid(y),unit=unit(y),
                      lon=lon(y),lat=lat(y),alt=NA,
                      cntr=NA,longname=NA,stid=NA,quality=NA,src=NA,
                      url=NA,reference=NA,info=NA, method="field",
                      type=NA,aspect=NA)
  attr(y,'history') <- history.stamp(x)
  class(y)[2] <- class(x)[2]
  if (dim(y)[2]==1) y <- subset(y,is=1)
  invisible(y)
}

as.station.spell <- function(x) {
  y <- coredata(x)
  ok1 <- is.finite(y[,1])
  above <- zoo(y[ok1,1],index(x)[ok1])
  ok2 <- is.finite(y[,2])
  below <- zoo(y[ok2,2],index(x)[ok2])
#  if (attr(x,'variable')=='t2m') suffix <- c("warm","cold") else
#                                 suffix <- c("wet","dry")
  suffix <- attr(x,'variable')
  y <- merge.zoo(above,below,suffixes=suffix,fill=0)
  #plot(y)
  y <- attrcp(x,y)
  attr(y,'unit') <- c("days","days")
  attr(y,'variable') <- suffix
  class(y) <- c('station',class(x))
  attr(y,'history') <- history.stamp(x)
  return(y)
}


as.station.eof <- function(x,pattern=1:10) {
  stopifnot(!missing(x),inherits(x,'eof'))
  z <- zoo(x[,pattern],order.by=index(x))
  y <- as.station.zoo(z,loc=paste('PC[',pattern,']',sep=''),
                      param=varid(x),unit=unit(x),
                      lon=lon(x),lat=lat(x),alt=NA,
                      cntr=NA,longname='principal component',stid=NA,quality=NA,
                      src=src(x),url=attr(x,'url'),reference=NA,info=NA, method="field",
                      type='index',aspect='anomaly')
  attr(y,'history') <- history.stamp(x)
  class(y)[2] <- class(x)[2]
  if (dim(y)[2]==1) y <- subset(y,is=1)
  invisible(y)
}


as.pca <- function(x) UseMethod("as.pca")

as.pca.ds <- function(x) {
  stopifnot(inherits(x,'pca'))
  class(x) <- class(x)[-1]
  return(x)
}

as.pca.station <- function(x) {
  y <- PCA(x)
  return(y)
}


as.ds <- function(x) {
  y <- zoo(X,order.by=index)
  attr(y,'location') <- attr(x,'location')
  attr(y,'variable') <- attr(x,'variable')
  attr(y,'unit') <- attr(x,'unit')
  attr(y,'longitude') <- attr(x,'longitude')
  attr(y,'latitude') <- attr(x,'latitude')
  attr(y,'altitude') <- attr(x,'altitude')
  attr(y,'country') <- attr(x,'country')
  attr(y,'longname') <- attr(x,'longname')
  attr(y,'station_id') <- attr(x,'station_id')
  attr(y,'quality') <- attr(x,'quality') 
  attr(y,'calendar') <- attr(x,'calendar')
  attr(y,'source') <- attr(x,'source')
  attr(y,'URL') <- attr(x,'URL')
  #attr(y,'history') <- attr(x,'history')
  #attr(y,'date-stamp') <- attr(x,'date-stamp')
  attr(y,'type') <- attr(x,'type')
  attr(y,'aspect') <- attr(x,'aspect')
  attr(y,'reference') <- attr(x,'reference')
  attr(y,'info') <- attr(x,'info')
  #attr(y,'call') <- match.call()
  attr(y,'history') <- history.stamp(x)
  class(y) <- c("station","month","zoo")
  return(y)
}

as.comb <- function(x,...) UseMethod("as.comb")

as.comb.eof <- function(x,...) {
  stopifnot(inherits(x,'eof'))
  #print("as.field.eof")
  y <- eof2field(as.eof(x))
  #print("HERE")
  n <- attr(x,'n.apps')
  if (!is.null(n)) {
    for ( i in 1:n ) {
      z <- as.eof(x,i)
      #plot(z)
      Z <- eof2field(z)
      eval(parse(text=paste("Z -> attr(y,'appendix.",i,"')",sep="")))
    }
  }
  n -> attr(y,'n.apps')
  attr(y,'history') <- history.stamp(x)
  attr(y,'quality') <- c(attr(x,'quality'),'EOF-filtered')
  return(y)
}


as.field <- function(x,...) UseMethod("as.field")

as.field.zoo <- function(x,lon,lat,param,unit,
                         longname=NA,quality=NA,src=NA,url=NA,
                         reference=NA,info=NA,calendar='gregorian',
                         greenwich=TRUE, method= NA,type=NA,aspect=NA) {
  #print("as.field.zoo")
  #print("lon"); print(lon); print("lat"); print(lat); print(param); print(unit)
  #print(c(length(lon),length(lat),length(index)))

  t <- index(x)
  #print(length(t))
  #dyr <- as.numeric(format(t[2],'%Y')) - as.numeric(format(t[1],'%Y')) 
  #dmo <- as.numeric(format(t[2],'%m')) - as.numeric(format(t[1],'%m')) 
  #dda <- as.numeric(format(t[2],'%d')) - as.numeric(format(t[1],'%d'))
  dyr <- diff(year(x))[1]
  dmo <- diff(month(x))[1]
  dda <- diff(day(x))[1]
  timescale <- "annual"
  if (dmo>0)  timescale <- "month"
  if (dmo==3)  timescale <- "season"
  if (dda>0)  timescale <- "day"
  #print(timescale)
  
# Add attributes to x
  attr(x,"variable") <- param
  attr(x,"longname") <- longname
  attr(x,"unit") <- unit
  attr(x,"source") <- src
  attr(x,"dimensions") <- c(length(lon),length(lat),length(t))
  attr(x,"longitude") <- lon
  attr(x,"latitude") <- lat
  attr(x,"greenwich") <- greenwich
  attr(x,"calendar") <- calendar 
  attr(x,'type') <- type
  attr(x,'aspect') <- aspect
  attr(x,'reference') <- reference
  attr(x,'info') <- attr(x,'info')
  attr(x,'history') <- history.stamp(x)
  class(x) <- c("field",timescale,"zoo")
  invisible(x)
}



as.field.default <- function(x,index,lon,lat,param,unit,
                         longname=NA,quality=NA,src=NA,url=NA,
                         reference=NA,info=NA,calendar='gregorian',
                         greenwich=TRUE, method= NA,type=NA,aspect=NA) {

#print("as.field.default")
#create a zoo object z
  z <- zoo(x=x,order.by=index)
  x <- as.field.zoo(z,lon=lon,lat=lat,param=param,unit=unit,
                    longname=longname,quality=quality,src=src,url=url,
                    reference=reference,info=info,calendar=calendar,
                    greenwich=greenwich, method=method,type=type,aspect=aspect)
  invisible(x)
}


as.field.field <- function(x,...) {
  if (inherits(x,'comb')) x <- as.field.comb(x,...)
  return(x)
}

as.field.comb <- function(x,iapp=NULL,...) {
  if (is.null(iapp)) {
    # Drop the appendend fields:
    n <- attr(x,'n.apps')
    for ( i in 1:n ) eval(parse(text=paste("attr(x,'appendix.",i,"') <- NULL",sep="")))
    attr(x,'n.apps') <- NULL
    class(x) <- class(x)[-1]
    return(x)
  }
  # Select one of the appended fields
  eval(parse(text=paste("y <- attr(x,'appendix.",iapp,"')",sep="")))
  attr(y,'longitude') <- attr(x,'longitude')
  attr(y,'latitude') <- attr(x,'latitude')
  attr(y,'history') <- history.stamp(x)
  return(y)  
}

as.field.eof <- function(x,iapp=NULL,...) {
  #print("as.field.eof")
  if (!inherits(x,'comb')) y <- eof2field(x) else {
   y <- as.eof(x,iapp)
   y <- eof2field(y)
 }
  return(y)
}

as.field.station <- function(x,lon=NULL,lat=NULL...) {
  if (is.null(lon)) lon <- seq(min(lon(x)),max(lon(x)),length=30)
  if (is.null(lat)) lat <- seq(min(lat(x)),max(lat(x)),length=30)
  y <- regrid.default(x,is=list(lon,lat))
  attr(y,'history') <- history.stamp(x)
  return(y)  
}






#as.annual <- function(x) {
#  year <- format(index), '%Y') 
#  return(year)
#  y <- annual(x)
#}
as.annual <- function(x, ...) UseMethod("as.annual")
#as.annual.default <- function(x, ...) as.annual(as.numeric(x))
as.annual.default <- function(x, ...) annual(x,...)
as.annual.numeric <- function(x, ...) annual(x)
as.annual.integer <- function(x, ...) structure(x, class = "annual")
as.annual.yearqtr <- function(x, frac = 0, ...) {
    if (frac == 0) annual(as.numeric(x)) else
    as.annual(as.Date(x, frac = frac), ...)
}
as.annual.station <- function(x, ...) annual.station(x,...)

as.annual.spell <- function(x, ...) annual.spell(x,...)

#as.annual <- function(x,FUN=mean,...) {
#  y <- annual(x,FUN=match.fun(FUN), ...)
#  return(y)
#}

as.monthly <- function(x,FUN='mean',...) {
  y <- aggregate(as.zoo(x),function(tt) as.Date(as.yearmon(tt)),FUN=FUN,...)
  y <- attrcp(x,y)
  attr(y,'history') <- history.stamp(x)
  class(y) <- class(x)
  class(y)[2] <- "month" 
  return(y)
}


# Not to confuse with season
# This function extracts a given seasonal interval and aggrigates a given statistic
as.seasons <- function(x,start='01-01',end='12-31',FUN='mean', ...) {
  IV <- function(x) sum(is.finite(x))
  yrs <- year(x); d <- dim(x)
  # ns = number of stations
  if (is.null(d)) ns <- 1 else ns <- d[2]
  years <- as.numeric(rownames(table(yrs))); n <- length(years)
  y <- matrix(rep(NA,n*ns),n,ns); k <- y
  start.1 <- as.numeric(as.Date(paste(years[1],start,sep='-')))
  end.1 <- as.numeric(as.Date(paste(years[1],end,sep='-')))
  if (start.1 > end.1) twoyears <- 1 else twoyears <- 0
  #browser()
  for (i in 1:n) {
    z <- coredata(window(x,start=as.Date(paste(years[i],start,sep='-')),
                           end=as.Date(paste(years[i]+twoyears,end,sep='-'))))
    k[i,] <- apply(matrix(z,length(z),ns),2,IV)
    y[i,] <- apply(matrix(z,length(z),ns),2,FUN, ...)
  }
  y <- zoo(y,order.by=as.Date(paste(years,start,sep='-')))
  y <- attrcp(x,y)
  attr(y,'history') <- history.stamp(x)
  if (twoyears==0) attr(y,'season.interval') <- paste(start,'to',end) else
                   attr(y,'season.interval') <- paste(start,'to',end,'the following year')
  attr(y,'n.valid') <- k
  class(y) <- class(x)
  class(y)[2] <- "annual"
  return(y)
}


## Original script : as.R
## Author          : Rasmus
## Major updates   : 2013-12-05 Abdelkader 

## Comments        : ! fourseasons() is no longer used 

as.4seasons <- function(x,...) UseMethod("as.4seasons")

# Aggregate for the seasons DJF (winter), MAM (spring), JJA (summer), and
# SON (autumn). Need to combine december with the following year's Jan-Feb.
# Best solution: fourseasons is a function applied to index(x)?



as.4seasons.default <- function(x,FUN='mean',slow=FALSE,...) {
  #print('as.4seasons.default')
  attr(x,'names') <- NULL

  d <- dim(coredata(x))
  #print(d)
  if (is.null(d)) d <- c(length(x),1)
  if (!slow) { 
    if (inherits(x,"month")) {
      if ( (is.null(d)) | (d[2]==1) )
        X <- c(NA,coredata(x)[1:length(x)-1]) # shift the coredata by 1 to start on December. This works only for monthly data !!!  
      else
        X <- rbind(rep(NA,d[2],1),coredata(x)[1:d[1]-1,])
      #print(dim(X))
      X <- zoo(X,order.by=index(x))
    ##yrseas <- fourseasons(ix)
    ##print(yrseas)
     #print('agrigate')
     #print(names(list(...)))

     y <- aggregate(x=as.zoo(X),by= as.yearqtr,FUN=match.fun(FUN),...)
    # convert yearqtr to yearmon
      y <- zoo(x=y,order.by=as.Date(as.yearmon(index(y))))
    } else y <- as.4seasons.day(x,FUN=FUN,...)
    #y <- as.4seasons.day(x,FUN=match.fun(FUN),...)
    
    #print(dim(y))
    ok <- length(index(y))
    #print(summary(c(coredata(y))))
  } else {
    yr <- sort(rep(as.integer(rownames(table(year(x)))),4))
    n <- length(yr)

    q <- rbind(c(12,1,2),3:5,6:8,9:11)
    X <- matrix(rep(NA,n*d[2]),n,d[2]) 
    t <-rep(NA,n)

    #print("start loop")
    for (i in 1:n) {
      iq <- (i-1) %% 4 + 1
      if (iq == 1) 
        ii <- (is.element(year(x),yr[i]) &
               is.element(month(x),c(1,2))) |
              (is.element(year(x),yr[i]-1) &
               is.element(month(x),12))  else
        ii <-  is.element(year(x),yr[i]) &
               is.element(month(x),q[iq,])
     if (d[2]==1) cline <- paste(FUN,"(coredata(x[ii,]),...)",sep="") else
                  cline <- paste("apply(coredata(x[ii,]),2,",FUN,", ...)",sep="")
      #print(cline)
      if ( (inherits(x,'day')) & (sum(ii)>=85) |
           (inherits(x,'month')) & (sum(ii)>=3) ) {
        X[i,] <- eval(parse(text=cline))
        t[i] <- yr[i] + (iq-1)/4
        print(c(i,yr[i],round(mean(X[i,]),2),iq,sum(ii),round(X[i],2),t[i]))
      }
    } 
    #print("end loop")
    ok <- is.finite(rowMeans(X,na.rm=TRUE)) & is.finite(t)
    #print(table(as.yearqtr(t[ok])))
    #print(summary(c(X)))
    y <- zoo(X[ok,],order.by=as.Date(as.yearqtr(t[ok])))
    #names(y) <- FUN
  }
  y <- attrcp(x,y)
  attr(y,'history') <- history.stamp(x)
  if (inherits(x,'field'))
    attr(y,'dimensions') <- c(attr(x,'dimensions')[1:2],sum(ok))
  class(y) <- class(x)
  class(y)[2] <- "season"
  #plot(y)
  return(y) 
}

as.4seasons.day <- function(x,FUN='mean',na.rm=TRUE,dateindex=TRUE,...) {

  IV <- function(x) sum(is.finite(x))

  #print('as.4seasons.day')
  attr(x,'names') <- NULL  
  t <- index(x)
  year <- as.numeric(format(t,'%Y'))
  month <- as.numeric(format(t,'%m'))
  day <- as.numeric(format(t,'%d'))
  #shift the time stamps by one month, sneaking December into the subsequent year
  month <- month + 1
  dec <- is.element(month,13)
  year[dec] <- year[dec] + 1
  month[dec] <- 1
  # Change the day to avoid warning that the calendar is wrong (e.g. due to
  # too many days in February). Since the data is aggregated, the exact day
  # in the month doesn't matter here.
  hour <- 12*(day - 2*trunc(day/2))
  day <- trunc(day/2) + 1
  #print(table(year)); print(table(month));
  #print(table(day)); print(table(hour));   
  #tshifted <- as.Date(paste(year,month,day,sep="-"))
  tshifted <-  ISOdate(year=year,month=month,day=day,hour=hour)
  #print(summary(tshifted))
  X <- zoo(coredata(x),order.by=tshifted)
  nd <- aggregate(X,as.yearqtr,FUN=IV)
  ok <- nd >= 85
  #browser()
  # Test for the presens of 'na.rm' in argument list - this is a crude fix and not a
  # very satisfactory one. Fails for FUN==primitive function.
  if (is.function(FUN)) test.na.rm <- FALSE else
      test.na.rm <- (sum(is.element(FUN,c('mean','min','max','sum','quantile')))>0)
  if ( (sum(is.element(names(formals(FUN)),'na.rm')==1)) | (test.na.rm) )
     y <- aggregate(X,as.yearqtr,FUN=match.fun(FUN),...,na.rm=na.rm) else
     y <- aggregate(X,as.yearqtr,FUN=match.fun(FUN),...)
  if (dateindex)
    y <- zoo(coredata(y[ok]),order.by=as.Date(index(y)[ok]))
  unit <- attr(y,'unit')
  y <- attrcp(x,y,ignore=c("unit","names"))
  unit -> attr(y,'unit')
  #str(y); print(unit)
  attr(y,'history') <- history.stamp(x)
  class(y) <- class(x)
  class(y)[2] <- "season"
  invisible(y)
}

as.4seasons.station <- function(x,FUN='mean',...) {
  #print('as.4seasons.station')
  y <- as.4seasons.default(x,FUN,...)
#  y <- attrcp(x,y)
#  attr(y,'history') <- history.stamp(x)
#  class(y) <- class(x)
#  class(y) <- gsub("month","season",class(x))
  return(y) 
}

as.4seasons.spell <- function(x,FUN='mean',...) {
  y <- as.4seasons.default(as.station(x),FUN,...)
#  y <- attrcp(x,y)
#  attr(y,'history') <- history.stamp(x)
#  class(y) <- class(x)
#  class(y) <- gsub("month","season",class(x))
  return(y) 
}


as.4seasons.field <- function(x,FUN='mean',...) {
  d <- attr(x,"dimensions")
  y <- as.4seasons.default(x,FUN,...)
  y <- attrcp(x,y)
  attr(y,'history') <- history.stamp(x)
  ## Update dimensions
  attr(y, "dimensions") <- c(d[1],d[2],dim(y)[1])
  class(y) <- class(x)
  class(y) <- gsub("month","season",class(x))
  return(y)
}



as.anomaly <- function(x,...) UseMethod("as.anomaly")

as.anomaly.default <- function(x,ref=NULL,monthly=NULL,na.rm=TRUE) {
# The argument monthly can be used to force the method to be
# julian-day regression-based or based on monthly mean

#  print('as.anomaly.default')
#  yr <- as.integer(format(index(x),'%Y'))
#  mon <- as.integer(format(index(x),'%m'))
#  dy <- as.integer(format(index(x),'%d'))
  yr <- year(x);  mon <- month(x);  dy <- day(x)
  if (is.null(ref))
    ref <- seq(min(yr),max(yr),by=1)
  # Check whether the object is a time series of monthly data.
  ndd <- length(table(dy))
  #print(ndd)
  if ( (ndd==1) & is.null(monthly) ) monthly <- TRUE else
                                     monthly <- FALSE
  y <- coredata(x)
  
  if (monthly) {
    # If monthly, subtract the
    #print(table(month))
    if (is.null(dim(x))) {
      #print("1D")
      clim <- rep(0,12)
      for (i in 1:12) {
        im <- is.element(mon,i)
        z <- mean(y[im & is.element(yr,ref)],na.rm=na.rm)
        clim[i] <- z
        y[im] <- y[im] - clim[i]
      }
      y <- zoo(y,order.by=index(x))
    } else {
      #print("2D")
      y <- t(y)
      clim <- rep(0,length(y[,1])*12); dim(clim) <- c(length(y[,1]),12)
      for (i in 1:12) {
        im <- is.element(mon,i)
        z <- rowMeans(y[,im & is.element(yr,ref)],na.rm=na.rm) 
        #print(c(length(z),length(clim[,1]))); plot(z)
        clim[,i] <- z
        y[,im] <- y[,im] - clim[,i]
        #print(c(i,sum(im),mean(clim[,i])))
      }
      y <- zoo(t(y),order.by=index(x))
   }
  } else {
    #print("daily")
    t0 <- julian(index(x)[is.element(yr,ref)]) -
          julian(as.Date(paste(yr[is.element(yr,ref)],"-01-01",sep="")))
    t <- julian(index(x)) -
         julian(as.Date(paste(yr,"-01-01",sep="")))
    c1 <- cos(pi*t0/365.25); s1 <- sin(pi*t0/365.25)
    c2 <- cos(2*pi*t0/365.25); s2 <- sin(2*pi*t0/365.25)
    c3 <- cos(3*pi*t0/365.25); s3 <- sin(3*pi*t0/365.25)
    c4 <- cos(4*pi*t0/365.25); s4 <- sin(4*pi*t0/365.25)
    C1 <- cos(pi*t/365.25); S1 <- sin(pi*t/365.25)
    C2 <- cos(2*pi*t/365.25); S2 <- sin(2*pi*t/365.25)
    C3 <- cos(3*pi*t/365.25); S3 <- sin(3*pi*t/365.25)
    C4 <- cos(4*pi*t/365.25); S4 <- sin(4*pi*t/365.25)
    cal <- data.frame(y=coredata(x),c1=c1,c2=c2,c3=c3,c4=c4,
                      s1=s1,s2=s2,s3=s3,s4=s4)
    pre <- data.frame(c1=C1,c2=C2,c3=C3,c4=C4,
                      s1=S1,s2=S2,s3=S3,s4=S4)
    i1 <- is.element(year(x),year(x)[1])
    pre1 <- data.frame(c1=C1[i1],c2=C2[i1],c3=C3[i1],c4=C4[i1],
                      s1=S1[i1],s2=S2[i1],s3=S3[i1],s4=S4[i1])
    acfit <- lm(y ~ c1 + s1 + c2 + s2 + c3 + s3 + c4 + s4,data=cal)
    clim <- predict(acfit,newdata=pre)
    y <- zoo(coredata(x) - clim,order.by=index(x))
    clim <-  predict(acfit,newdata=pre1)
  }
  #print("attributes")
  y <- attrcp(x,y)
  attr(y,"aspect") <- 'anomaly'
  attr(y,"anomaly_method") <- monthly
  attr(y,"climatology") <- clim
  #attr(y,"date") <- date()
  #attr(y,"call") <- match.call()
  attr(y,'history') <- history.stamp(x)
  class(y) <- class(x)
  invisible(y)
}

as.anomaly.zoo <- function(x,ref=NULL,monthly=NULL,na.rm=TRUE) {
  y <- as.anomaly.station(x,ref,monthly,na.rm)
  attr(y,'history') <- history.stamp(x)
  invisible(y)
}

as.anomaly.station <- function(x,ref=NULL,monthly=NULL,na.rm=TRUE) {
  y <- as.anomaly.default(x,ref,monthly,na.rm)
  attr(y,'history') <- history.stamp(x)
  invisible(y)
}

as.anomaly.field<- function(x,ref=NULL,monthly=NULL,na.rm=TRUE) {
   y <- as.anomaly.default(x,ref,monthly,na.rm)
   attr(y,'history') <- history.stamp(x)
   attr(y,'dimensions') <- attr(x,'dimensions')
   invisible(y)
}

# Handy conversion algorithms:
as.climatology <- function(x,...) {
  ya <- as.anomaly(x)
  y <- attr(ya,'climatology')
  y <- attrcp(x,y)
  attr(y,'history') <- history.stamp(x)
  class(y) <- class(x)
  invisible(y)
}

as.residual <- function(x) UseMethod("as.residual")

as.residual.ds <- function(x){
  y <- attr(x,'original_data') - attr(x,'fitted_values')
  y <- attrcp(x,y)
  attr(y,'aspect') <- 'residual'
  attr(y,'history') <- history.stamp(x)
  class(y) <- class(attr(x,'calibration_data'))
  invisible(y)
}

as.residual.station <- function(x){
  if (!is.null(attr(x,'calibration_data')))
    y <- as.residual.ds(x) else y <- NULL
  invisible(y)
}


as.calibrationdata <- function(x) UseMethod("as.calibrationdata")

as.calibrationdata.ds <- function(x) {
  y <- attr(x,'calibration_data')
  y <- attrcp(x,y)
  attr(y,'history') <- history.stamp(x)
  class(y) <- class(attr(x,'calibration_data'))
  invisible(y)
}

as.calibrationdata.station <- function(x) {
   if (!is.null(attr(x,'calibration_data')))
    y <- as.calibrationdata.ds(x) else y <- NULL
  invisible(y)
}

as.fitted.values <- function(x) UseMethod("as.fitted.values")

as.fitted.values.ds <- function(x) {
  y <- attr(x,'fitted_values')
  y <- attrcp(x,y)
  attr(y,'history') <- history.stamp(x)
  class(y) <- class(attr(x,'fitted_values'))
  invisible(y)
}

as.fitted.values.station <- function(x) {
   if (!is.null(attr(x,'fitted.values')))
    y <- as.fitted.values.ds(x) else y <- NULL
  invisible(y)
}



as.original.data <- function(x) UseMethod("as.original.data")

as.original.data.ds <- function(x) {
  y <- attr(x,'original_data')
  y <- attrcp(x,y)
  attr(y,'history') <- history.stamp(x)
  class(y) <- class(attr(x,'original_data'))
  invisible(y)
}

as.original.data.station <- function(x) {
  y <- as.original.data.ds(x)
  invisible(y)
}


as.pattern <- function(x,...) UseMethod("as.pattern")

as.pattern.ds <- function(x) {
  y <- attr(x,'pattern')
  attr(y,'history') <- history.stamp(x)
  invisible(y)
}

as.pattern.eof <- function(x) {
  y <- attr(x,'pattern')
  attr(y,'history') <- history.stamp(x)
  invisible(y)  
}

as.pattern.mvr <- function(x) {
  y <- attr(x,'pattern')
  invisible(y)  
}

as.pattern.cca <- function(x) {
  y <- attr(x,'pattern')
  attr(y,'history') <- history.stamp(x)
  invisible(y)  
}

as.pattern.trend <- function(x) {
  y <- attr(x,'pattern')
  attr(y,'history') <- history.stamp(x)
  invisible(y)  
}

as.pattern.field <- function(x,FUN=NULL,...) {
  if (!is.null(FUN)) {
    y <- apply(x,2,FUN,...)
    dim(y) <- attr(x,'dimension')[1:2]
  } else {
    y <- t(coredata(x))
    dim(y) <- attr(x,'dimension')
  }
  attr(y,'longitude') <- attr(x,'longitude')
  attr(y,'latitude') <- attr(x,'latitude')
  attr(y,'time') <- index(x) 
  attr(y,'history') <- history.stamp(x)
  invisible(y)
}

as.pattern.corfield <- function(x) {
  y <- coredata(x)
  dim(y) <- attr(x,'dimension')[1:2]
  y <- attrcp(x,y)
  attr(y,'longitude') <- lon(x)
  attr(y,'latitude') <- lat(x)
  attr(y,'time') <- attr(x,'time')
  class(y) <- 'matrix'
  
  invisible(y)
}

as.eof <- function(x,...) UseMethod("as.eof")

as.eof.zoo <- function(x,...) {
  class(x) <- c('eof','zoo')
  return(x)
}

as.eof.eof <-function(x,iapp=NULL) {
  if (inherits(x,'comb')) x <- as.eof.comb(x,iapp) 
  return(x)
}
  
as.eof.comb <- function(x,iapp=NULL) {
  #print("as.eof.comb")
  stopifnot(inherits(x,'comb'))

  # if x is a 'field'
  if (!inherits(x,'eof')) x <- EOF(x)

  # assume x from now on is an 'eof'
  if (!is.null(iapp)) {
    y <- as.eof.appendix(x,iapp)
    return(y)
  }
  class(x) <- class(x)[-grep('comb',class(x))]
  napps <- attr(x,'n.apps')
  for (i in napps) {
    eval(parse(text=paste("attr(x,'appendix.",i,"') <- NULL",sep="")))
  }
  attr(x,'n.apps') <- NULL
  attr(x,'history') <- history.stamp(x)
  return(x)
}

as.eof.field <- function(x,iapp=NULL,...) {
  y <- EOF(x,...)
  if (!is.null(iapp)) y <- as.eof.appendix(y,iapp)
  return(y)
}

as.eof.appendix <- function(x,iapp=1) {
  #print("as.eof.appendix")
  stopifnot(inherits(x,'comb'))
  y <- eval(parse(text=paste("attr(x,'appendix.",iapp,"')",sep="")))
  x <- as.eof.comb(x)
  y <- attrcp(x,y)
  attr(y,'history') <- history.stamp(x)
  class(y) <- class(x)
  return(y)
}


as.appended <- function(x,...) UseMethod("as.appended")

as.appended.ds.comb <- function(x,iapp=1) {
  eval(parse(text=paste("X <- attr(x,'appendix.",it,"')",sep="")))
  X <- attrcp(x,X,ignore='appendix')
  attr(X,'history') <- history.stamp(x)
  invisible(X)
}

as.appended.eof.comb <- function(x,iapp=1) {
  X <- as.appended.ds.comb(x,it=it)
  invisible(X)
}

as.appended.field.comb <- function(x,iapp=1) {
  X <- as.appended.ds.comb(x,it=it)
  invisible(X)
}

as.stand <- function(x,...) UseMethod("as.stand")

as.stand.station <- function(x) {
  if (is.precip(x)) {
    mu <- apply(x,2,mean,na.rm=na.rm)
    X <- 100*x/mu
    attr(X,'clim') <- mu
    attr(X,'aspect') <- 'proportional'
    attr(X,'unit') <- '%'
    attr(X,'oldunit') <- attr(x,'unit')
  } else if (is.T(x)) {
    mu <- apply(x,2,mean,na.rm=na.rm)
    sigma <- apply(x,2,sd,na.rm=na.rm)
    X <- (x - mu)/sigma
    attr(X,'mean') <- mu
    attr(X,'sigma') <- sigma
    attr(X,'aspect') <- 'standardised'    
  }
  attr(X,'history') <- history.stamp(x)
  return(X)
}



as.original <- function(x) UseMethod("as.original")

as.original.station <- function(x) {
  if (attr(x,'aspect')=='proportional') {
    X <- attr(x,'clim')*x/100
    attr(X,'clim') <- NULL
    attr(X,'unit') <- attr(x,'oldunit')
    attr(X,'oldunit') <- NULL
    attr(X,'aspect') <- 'original'
  } else if (attr(x,'aspect')=='standardised') {
    X <- x * attr(x,'sigma') + attr(x,'mean')
    attr(X,'mean') <- NULL
    attr(X,'sigma') <- NULL
    attr(X,'aspect') <- 'original'     
  } else X <- x
  attr(X,'history') <- history.stamp(x)
  return(X)
}
