
trajectory2field <- function(x,dt='month',dx=2,dy=2,n=150,
                             it=NULL,is=NULL,verbose=FALSE) {
  stopifnot(is.trajectory(x))
  x <- subset(x,it=it,is=is)
  lons <- round(x[,colnames(x)=='lon']/dx)*dx
  lons <- seq(min(lons),max(lons),dx)
  lats <- round(x[,colnames(x)=='lat']/dy)*dy
  lats <- seq(min(lats),max(lats),dy)
  dates <- as.Date(strptime(x[,colnames(x)=='start'],format='%Y%m%d%H'))
  fn <- function(x,it=NULL,verbose=FALSE) {
    x <- subset(x,it=it)
    X <- array(rep(0,),dim=c(length(lons),length(lats)))
    if(verbose)print(dim(x))
    if(!is.null(dim(x)) & length(x)>0) {
      A <- density.trajectory(x,dx=dx,dy=dy,n=n)
      lat <- A$lat
      lon <- A$lon
      den <- A$density
      for(j in 1:length(lat)) {
        X[lons==lon[j],lats==lat[j]] <- den[j]	
      }
      dim(X) <- length(X)
    } else {
        if(verbose) print("no storms")
    }
    invisible(X)
  }
  if (grepl('month',dt)) {   
    dvec <- seq(min(as.Date(as.yearmon(dates))),max(dates),by='month')
    dall <- as.Date(as.yearmon(dates))
    X <- t(sapply(dvec,function(di) fn(x,it=(dall==di))))
    unit <- 'events/month/area'
  } else if (grepl('season',dt) | grepl('quarter',dt)) {
    dvec <- seq(min(as.Date(as.yearmon(dates))),max(dates),by='quarter')
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

factor2numeric <- function(f) {
  if(!is.null(levels(f))) {return(as.numeric(levels(f))[f])
  } else return(as.numeric(f))
}

density.trajectory <- function(x,it=NULL,is=NULL,dx=2,dy=2,n=150,R=6378) {
  stopifnot(is.trajectory(x))
  x <- subset(x,it=it,is=is)
  xbin <- bin.trajectory(x,dx=dx,dy=dy,n=n)
  lon <- factor2numeric(xbin$lon)
  lat <- factor2numeric(xbin$lat)
  freq <- factor2numeric(xbin$Freq)
  A <- dx*(pi/180)*R**2*abs(sin((lat+dy/2)*pi/180)-
                            sin((lat-dy/2)*pi/180))
  dens <- freq/A
  dens[is.infinite(dens)] <- NA
  #if(any(lat>88)) {
  #  f90 <- sum(freq[lat>88])
  #  A90 <- 2*pi*R**2*(1-sin(88*pi/180))
  #  d90 <- f90/A90
  #  lon90 <- seq(-180,180,dx)
  #  lon <- c(lon90,lon[lat<=88])
  #  dens <- c(rep(d90,length(lon90)),freq[lat<=88])
  #  lat <- c(rep(89,length(lon90)),lat[lat<=88])
  #}
  X <- data.frame(lon=lon,lat=lat,density=dens)
  invisible(X)
}

## surfacearea <- function(lon,lat,dx,dy,R=6378,radians=FALSE) {
##   if(!radians) {
##     lat <- lat*pi/180
##     lon <- lon*pi/180
##     dx <- dx*pi/180
##     dy <- dy*pi/180
##   }
##   return(R**2*cos(lat)*dx*2*sin(dy/2))
## }

bin.trajectory <- function(x,it=NULL,is=NULL,dx=2,dy=2,n=150) {
  stopifnot(is.trajectory(x))
  y <- subset(x,it=it,is=is)
  A <- apply(y,1,function(x) gridbin(x[colnames(y)=='lon'],
                      x[colnames(y)=='lat'],dx=dx,dy=dy,n=n))
  lon <- unlist(lapply(A,function(x) factor2numeric(x$x)))
  lat <- unlist(lapply(A,function(x) factor2numeric(x$y)))
  hits <- as.data.frame(table(lon,lat))
  invisible(hits[hits$Freq>0,])
}

gridbin <- function(x,y,dx=1,dy=1,n=100) {
  xy <- approx.lonlat(x,y,n=n)
  xx <- xy[,1]
  yy <- xy[,2]
  xx <- round(xx/dx)*dx
  yy <- round(yy/dy)*dy
  xy <- unique(cbind(xx,yy))
  xx <- xy[,1]
  yy <- xy[,2]
  ## if(all(yy==90)) {
  ##   xx <- fnlon(mean)(xx)
  ##   yy <- 90
  ## } else if(any(yy==90) & any(yy<90)) {
  ##   xx <- c(xx[yy<90],fnlon(mean)(xx))
  ##   yy <- c(yy[yy<90],90)
  ## }
  hit <- as.data.frame(table(xx,yy))
  hit <- hit[hit$Freq>0,]
  invisible(data.frame(x=hit$xx,y=hit$yy))
}

approx.lonlat <- function(lon,lat,n=10,a=6.378e06,verbose=FALSE) {
  if (verbose) print("approx.lonlat")
  x <- a * cos( lat*pi/180 ) * cos( lon*pi/180 )
  y <- a * cos( lat*pi/180 ) * sin( lon*pi/180 )
  z <- a * sin( lat*pi/180 )
  xa <- approx(x,n=n)$y
  ya <- approx(y,n=n)$y
  za <- approx(z,n=n)$y
  lona <- atan2( ya, xa )*180/pi
  lata <- asin( za/sqrt( xa^2 + ya^2 + za^2 ))*180/pi
  invisible(cbind(lona,lata))
}
