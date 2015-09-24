events2field <- function(x,...) {
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

 

map.events <- function(x,it=NULL,is=NULL,
               projection="sphere",verbose=TRUE,...) {
  y <- subset(x,it=it,is=is)
  Y <- as.field(y)
  map(Y,...)  
}
