trajectory2field <- function(x,dt='month',dx=2,dy=2,radius=5E5,
                             it=NULL,is=NULL,verbose=FALSE) {
  if(verbose) print("trajectory2field")
  stopifnot(is.trajectory(x))
  x <- subset(x,it=it,is=is)
  lons <- round(x[,colnames(x)=='lon']/dx)*dx
  lons <- seq(min(lons),max(lons),dx)
  lats <- round(x[,colnames(x)=='lat']/dy)*dy
  lats <- seq(min(lats),max(lats),dy)
  dates <- strptime(x[,'start'],format='%Y%m%d%H')
  if(is.null(attr(x,"calendar"))) calendar <- "gregorian" else calendar <- attr(x,"calendar")
  if (requireNamespace("PCICt", quietly = TRUE)) {
    dates <- PCICt::as.PCICt(as.character(x[,'start']),format='%Y%m%d%H',cal=calendar)
  } else {
    dates <- as.POSIXct(as.character(x[,'start']),format='%Y%m%d%H')
  }
  fn <- function(x,it=NULL) {
    x <- subset(x,it=it)
    X.out <- array(rep(0,),dim=c(length(lons),length(lats)))
    if(length(x)>0) {
      if(verbose) {
        if(!is.null(dim(x))) {
          print(paste(range(x[,'start']),collapse="-"))
          print(dim(x))
        } else {
          print(x['start'])
        }
      }
      if(!is.null(dim(x))) {
        A <- density.trajectory(x,dx=dx,dy=dy,radius=radius)
        lat <- A$lat
        lon <- A$lon
        den <- A$density
        for(j in 1:length(lat)) {
          X.out[lons==lon[j],lats==lat[j]] <- den[j]
        }
        dim(X.out) <- c(length(lons),length(lats))#length(X)
      } else {
        if(verbose) print("no storms")
      }
    }
    invisible(X.out)
  }
  if(verbose) print("calculate trajectory density")
  if (grepl('month',dt)) {   
    if(verbose) print("monthly")
    dall <- as.Date(paste(format(dates,format="%Y-%m"),"01",sep="-"))
    dvec <- seq(min(dall),max(dall),by='month')
    unit <- 'events/month/area'
    fn2 <- function(di) fn(x,it=(dall==di))
    #X <- t(sapply(dvec,function(di) fn2))
  } else if (grepl('season',dt) | grepl('quarter',dt)) {
    if(verbose) print("seasonal")
    dall <- as.Date(as.yearqtr(as.Date(paste(format(dates,format="%Y-%m"),"01",sep="-"))))
    dvec <- seq(min(dall),max(dall),by='quarter')
    fn2 <- function(di) fn(x,it=(dall==di))
    #X <- t(sapply(dvec,fn2))
    unit <- 'events/quarter/area'
  } else if (grepl('year',dt) | grepl('annual',dt)) {
    if(verbose) print("annual")
    dall <- as.Date(paste(format(dates,format="%Y"),"01-01",sep="-"))
    dvec <- seq(min(dall),max(dall),by='year')
    unit <- 'events/year/area'
    fn2 <- function(di) fn(x,it=c(year(di),year(di)))
    #X <- t(sapply(dvec,fn2)
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
    di <- try(fn2(dvec[i]))
    X[i,,] <- di
  }
  t2 <- Sys.time()
  if (verbose) print(paste('Calculating the density took',
                round(as.numeric(t2-t1,units="secs")),'s'))
  d <- dim(X)
  dim(X) <- c(d[1],d[2]*d[3])
  longname <- paste("trajectory density",attr(x,'longname'),sep=', ')
  param <- 'density'
  if(verbose) print("transform to field")
  Y <- as.field(X,index=dvec,lon=lons,lat=lats,
          unit=unit,longname=longname,param=param,
          quality=attr(x,'quality'),src=attr(x,'source'),
          url=attr(x,'URL'),reference=attr(x,'reference'),
          info=attr(x,'info'),calendar=attr(x,'calendar'),
          method=attr(x,'method'),aspect=attr(x,'aspect'))
  if(verbose) print("done")
  invisible(Y)
}

## density.trajectory <- function(x,it=NULL,is=NULL,dx=2,dy=2,R=6378) {
##   stopifnot(is.trajectory(x))
##   x <- subset(x,it=it,is=is)
##   xbin <- bin.trajectory(x,dx=dx,dy=dy,n=n)
##   lon <- factor2numeric(xbin$lon)
##   lat <- factor2numeric(xbin$lat)
##   freq <- factor2numeric(xbin$Freq)
##   A <- dx*(pi/180)*R**2*abs(sin((lat+dy/2)*pi/180)-
##                             sin((lat-dy/2)*pi/180))
##   dens <- freq/A
##   dens[is.infinite(dens)] <- NA
##   X <- data.frame(lon=lon,lat=lat,density=dens)
##   invisible(X)
## }

## bin.trajectory <- function(x,it=NULL,is=NULL,dx=2,dy=2,n=20) {
##   stopifnot(is.trajectory(x))
##   y <- subset(x,it=it,is=is)
##   A <- apply(y,1,function(x) gridbin(x[colnames(y)=='lon'],
##                       x[colnames(y)=='lat'],dx=dx,dy=dy,n=n))
##   lon <- unlist(lapply(A,function(x) factor2numeric(x$x)))
##   lat <- unlist(lapply(A,function(x) factor2numeric(x$y)))
##   hits <- as.data.frame(table(lon,lat))
##   invisible(hits[hits$Freq>0,])
## }

## gridbin <- function(x,y,dx=1,dy=1,n=20) {
##   xy <- approx.lonlat(x,y,n=n)
##   xx <- xy[,1]
##   yy <- xy[,2]
##   xx <- round(xx/dx)*dx
##   yy <- round(yy/dy)*dy
##   xy <- unique(cbind(xx,yy))
##   xx <- xy[,1]
##   yy <- xy[,2]
##   hit <- as.data.frame(table(xx,yy))
##   hit <- hit[hit$Freq>0,]
##   invisible(data.frame(x=hit$xx,y=hit$yy))
## }



