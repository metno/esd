events2field <- function(x,verbose=FALSE,...) {
  if (verbose) print("events2field")
  y <- density.events(x,verbose=verbose,...)
  invisible(y)  
}

density.events <- function(x,dt="month",dx=2,dy=2,
                         lons=NULL,lats=NULL,it=NULL,is=NULL,
                         radius=TRUE,verbose=FALSE,...) {
  if (verbose) print("density.events")
  ok <- !is.na(x["time"][[1]]) & !is.na(x["lon"][[1]]) & !is.na(x["lat"][[1]])
  x <- subset(x,it=ok)
  y <- subset(x,it=it,is=is)
  if(is.null(lons)) lons <- y["lon"][[1]]
  if(is.null(lats)) lats <- y["lat"][[1]]
  if(is.null(dx)) dx <- min(diff(sort(unique(lons))))
  if(is.null(dy)) dy <- min(diff(sort(unique(lats))))
  lons <- round(lons/dx)*dx
  lons <- seq(min(lons),max(lons),dx)
  lats <- round(lats/dy)*dy
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
    d <- as.POSIXct(strptime(paste(y["date"][[1]],y["time"][[1]]),
                             format="%Y%m%d %H"))
    dh <- min(diff(sort(unique(y["time"][[1]]))))
    dvec <- seq(min(d),max(d),by=dh*60*60)
    unit <- paste('events/',dh,'hours/km^2',sep='')
  } else {
    print(paste("WARNING! invalid time resolution dt",dt))
    break
  }
  X <- array(rep(0,),dim=c(length(dvec),length(lons),length(lats)))
  if (verbose) print("looping...")
  t1 <- Sys.time()
  pb <- txtProgressBar(style=3)
  for (i in 1:length(dvec)) {
    setTxtProgressBar(pb,i/length(dvec)) 
    yi <- subset.events(y,it=(d==dvec[i]))
    if (!is.null(dim(yi)) & dim(yi)[1]>0) {
      if (radius) {
        di <- lonlatdensity(yi["lon"][[1]],yi["lat"][[1]],
         yi["radius"][[1]],dx=dx,dy=dy,verbose=verbose)
      } else {
         di <- lonlatdensity(yi["lon"][[1]],yi["lat"][[1]],
          dx=dx,dy=dy,verbose=verbose)
      }
      for(j in 1:length(di$lat)) {
        X[dvec==dvec[i],lons==di$lon[j],lats==di$lat[j]] <- di$density[j]
      }
    }
  }
  t2 <- Sys.time()
  if (verbose) print(paste('Calculating the density took',
                round(as.numeric(t2-t1,units="secs")),'s'))
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

lonlatdensity <- function(lons,lats,radius=NULL,dx=NULL,dy=NULL,
                          R=6371,verbose=FALSE,...) {
  if (!is.null(dx)) {
    lons <- round(lons/dx)*dx
  } else {
    dx <- min(diff(sort(unique(lons))))
  }
  if (!is.null(dy)) {
    lats <- round(lats/dy)*dy
  } else {
    dy <- min(diff(sort(unique(lats))))
  }
  fn <- function(x,y,r,dx,dy) {
    xvec <- seq(round(-180/dx)*dx,round(180/dx)*dx,dx)
    yvec <- seq(round(-90/dy)*dy,round(90/dy)*dy,dy)
    xx <- as.vector(sapply(xvec,function(x) rep(x,length(yvec))))
    yy <- rep(yvec,length(xvec))
    i <- mapply(function(a,b,d) distAB(a,b,xx,yy)<d,x,y,r)
    rx <- unlist(apply(i,2,function(j) xx[j]))
    ry <- unlist(apply(i,2,function(j) yy[j]))
    ok <- !is.na(rx) & !is.na(ry); rx <- rx[ok]; ry <- ry[ok]
    hits <- as.data.frame(table(lon=rx,lat=ry))
    hits <- hits[hits$Freq>0,]
    hx <- factor2numeric(hits$lon)
    hy <- factor2numeric(hits$lat)
    A <- dx*(pi/180)*R**2*abs(sin(min(hy+dy/2,90)*pi/180))-
                           sin((hy-dy/2)*pi/180)
    d <- hits$Freq/A
    d[is.infinite(d)] <- NA
    invisible(cbind(hx,hy,d))
  }
  if (is.null(radius)) radius <- 1E4
  dens <- fn(lons,lats,radius,dx,dy)
  hlons <- dens[,1]
  hlats <- dens[,2]
  d <- dens[,3]
  X <- data.frame(lon=hlons,lat=hlats,density=d)
  invisible(X)
}
