
trajectory2field <- function(x,dt='month',dx=2,dy=2,n=150,it=NULL,is=NULL) {
  stopifnot(is.trajectory(x))
  x <- subset(x,it=it,is=is)
  lons <- round(x[,colnames(x)=='lon']/dx)*dx
  lons <- seq(min(lons),max(lons),dx)
  lats <- round(x[,colnames(x)=='lat']/dy)*dy
  lats <- seq(min(lats),max(lats),dy)
  dates <- as.Date(strptime(x[,colnames(x)=='start'],format='%Y%m%d%H'))
  fn <- function(x,it=NULL) {
    x <- subset(x,it=it)
    X <- array(rep(0,),dim=c(length(lons),length(lats)))
    if(!is.null(dim(x)) & length(x)>0) {
      A <- density.trajectory(x,dx=dx,dy=dy,n=n)
      lat <- A$lat
      lon <- A$lon
      den <- A$density
      for(j in 1:length(lat)) {
        X[lons==lon[j],lats==lat[j]] <- den[j]	
      }
      dim(X) <- length(X)
    }
    invisible(X)
  }
  if (grepl('month',dt)) {   
    dvec <- seq(min(dates),max(dates),by='month')
    dall <- as.Date(as.yearmon(dates))
    X <- t(sapply(dvec,function(di) fn(x,it=(dall==di))))
    unit <- 'events/month/area'
  } else if (grepl('season',dt) | grepl('quarter',dt)) {
    dvec <- seq(min(dates),max(dates),by='quarter')
    dall <- as.Date(as.yearqtr(dates))
    X <- t(sapply(dvec,function(di) fn(x,it=(dall==di))))
    unit <- 'events/quarter/area'
  } else if (grepl('year',dt) | grepl('annual',dt)) {
    dvec <- seq(min(dates),max(dates),by='year')
    X <- t(sapply(year(dvec),function(yi) fn(x,it=c(yi,yi))))
    unit <- 'events/year/area'
  } else {
    print(paste("WARNING! invalid time resolution dt",dt))
    break
  }
  longname <- paste("trajectory density",attr(x,'longname'),sep=', ')
  param <- 'density'
  Y <- as.field(X,index=dvec,lon=lons,lat=lats,
          unit=unit,longname=longname,param=param,
          quality=attr(x,'quality'),src=attr(x,'source'),
          url=attr(x,'URL'),reference=attr(x,'reference'),
          info=attr(x,'info'),calendar=attr(x,'calendar'),
          method=attr(x,'method'),aspect=attr(x,'aspect'))
  invisible(Y)
}

density.trajectory <- function(x,it=NULL,is=NULL,dx=1,dy=1,n=150) {
  xbin <- bin.trajectory(x,it=it,is=is,dx=dx,dy=dy,n=n)
  lon <- as.numeric(levels(xbin$lon))[xbin$lon]
  lat <- as.numeric(levels(xbin$lat))[xbin$lat]
  freq <- as.numeric(xbin$Freq)
  A <- mapply(surfacearea,lon,lat,dx,dy)
  dens <- freq/A*1000
  dens[is.infinite(dens)] <- NA
  X <- data.frame(lon=lon,lat=lat,density=dens)
  invisible(X)
}

surfacearea <- function(lon,lat,dx,dy,R=6371) {
  return(R**2*( (lon+dx/2)*pi/180 - (lon-dx/2)*pi/180) *
       ( sin((lat+dy/2)*pi/180) - sin((lat-dy/2)*pi/180)) )
}

bin.trajectory <- function(x,it=NULL,is=NULL,dx=1,dy=1,n=150) {
  y <- subset(x,it=it,is=is)
  A <- apply(y,1,function(x) gridbin(x[colnames(y)=='lon'],
                      x[colnames(y)=='lat'],dx=dx,dy=dy,n=n))
  lon <- unlist(lapply(A,function(x) x$x))
  lat <- unlist(lapply(A,function(x) x$y))
  hits <- as.data.frame(table(lon,lat))
  invisible(hits[hits$Freq>0,])
}

gridbin <- function(x,y,dx=1,dy=1,n=150) {
  x <- approxlon(x,n=n)$y
  y <- approxlon(y,n=n)$y
  hit <- as.data.frame(table(round(x/dx)*dx,round(y/dy)*dy))
  hit <- hit[hit$Freq>0,]
  invisible(data.frame(x=hit$Var1,y=hit$Var2))
}
