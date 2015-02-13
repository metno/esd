## Author 	 Kajsa Parding
## Last update   13.02.2015
## Require 	 geoborders.rda

# map.storm: Storm tracks that pass the date line
# are excluded because they look terrible. Should they be
# split and plotted as two separate parts on each side of dateline?

map.storm <- function(x,it=NULL,is=NULL,FUN=NULL,
      projection="sphere",lonR=10,latR=90,
      col='red',colmap='rainbow',alpha=0.3,pfit=FALSE,
      xlim=NULL,ylim=NULL,new=TRUE) {

  y <- subset.storm(x,it=it,is=is)
  if (pfit) {
    lats <- t(polyfit.storm(y))
    y[,colnames(y)=='lat'] <- lats
  }

  if (!is.null(FUN)) {
    z <- round(FUN(y))
    uz <- sort(unique(z))
    dz <- min(uz[2:length(uz)]-uz[1:(length(uz)-1)])
    zlist <- seq(min(z),max(z),dz)
    clist <- colscal(n=length(zlist),col=colmap)
    col <- sapply(z,function(x) clist[zlist==x])
  }
  
  if (projection=="sphere" | projection=="np" | projection=="sp") {
    if (projection=="np") latR <- 90
    if (projection=="sp") latR <- -90
    sphere.storm(y,new=new,
    lonR=lonR,latR=latR,col=col,alpha=alpha)
  } else if (projection=="latlon" | projection=="lonlat") {
    lonlat.storm(y,new=new,
    xlim=xlim,ylim=ylim,col=col,alpha=alpha)
  }
}


lonlat.storm <- function(x,
    xlim=NULL,ylim=NULL,col='blue',alpha=0.1,
    size=1/16,lty=1,lwd=1,new=TRUE) {
  
  x0 <- x
  lons <- x[,colnames(x)=='lon']
  lats <- x[,colnames(x)=='lat']
  if (is.null(xlim)) xlim <- range(lons)
  if (is.null(ylim)) ylim <- range(lats)

  data("geoborders",envir=environment())
  ok <- is.finite(geoborders$x) & is.finite(geoborders$y)
  mlon <- geoborders$x[ok]
  mlat <- geoborders$y[ok]
  
  if (new) dev.new()
  par(bty="n",xaxt="n",yaxt="n",new=TRUE)
  plot(mlon,mlat,pch=".",col="white",
    xlab="lon",ylab="lat",xlim=xlim,ylim=ylim)

  OK <- apply(lons,1,function(x) !(any(x < -90) & any(x > 90)))
  matlines(t(lons[OK,]),t(lats[OK,]),lty=lty,lwd=lwd,
           col=adjustcolor(col,alpha.f=alpha))
  #matlines(t(lons[!OK,]),t(lats[!OK,]),col='red',lty=1)

  points(mlon,mlat,pch=".")
}


sphere.rotate <- function(lon,lat,lonR=0,latR=90) {
  theta <- pi*lon/180
  phi <- pi*lat/180
  x <- sin(theta)*cos(phi)
  y <- cos(theta)*cos(phi)
  z <- sin(phi)
  a <- rotM(x=0,y=0,z=lonR) %*% rbind(x,y,z)
  a <- rotM(x=latR,y=0,z=0) %*% a
  invisible(a)
}


sphere.storm <- function(x,
    xlim=NULL,ylim=NULL,col='blue',alpha=0.1,
    lty=1,lwd=1,lonR=0,latR=90,new=TRUE) {
  
  x0 <- x
  ilons <- colnames(x)=='lon'
  ilats <- colnames(x)=='lat'
  fn <- function(x) sphere.rotate(x[ilons],x[ilats],lonR=lonR,latR=latR)
  A <- apply(x,1,fn)
  n <- length(ilons[ilons])
  X <- A[seq(1,3*n,3),]
  Y <- A[seq(2,3*n,3),]
  Z <- A[seq(3,3*n,3),]
  X[Y<=0] = NA; Z[Y<=0] <- NA

  data("geoborders",envir=environment())
  ok <- is.finite(geoborders$x) & is.finite(geoborders$y)
  mlon <- geoborders$x[ok]
  mlat <- geoborders$y[ok]
  a <- sphere.rotate(mlon,mlat,lonR=lonR,latR=latR)
  x <- a[1,]; y <- a[2,]; z <- a[3,]
    
  if (new) dev.new()
  par(bty="n",xaxt="n",yaxt="n",new=TRUE)
  plot(x[y>0],z[y>0],pch=".",type="n",xlab="",ylab="")

  OK <- apply(x0[,ilons],1,function(x) !(any(x < -90) & any(x > 90)))
  j <- 1
  lines(X[,OK][,j],Z[,OK][,j],lty=lty,lwd=lwd,
        col=adjustcolor(col,alpha.f=alpha))
  matlines(X[,OK],Z[,OK],lty=lty,lwd=lwd,
           col=adjustcolor(col,alpha.f=alpha))
  # crossing the dateline messes everything up: 
  #matlines(X[,!OK],Z[,!OK],col='black',lty=1)

  points(x[y>0],z[y>0],pch=".")
  lines(cos(pi/180*1:360),sin(pi/180*1:360),col="black")
}


petals.storm <- function(x,dlon=5,dlat=2,digits=6) {
  x0 <- x
  lats <- x[,colnames(x)=='lat']
  lons <- x[,colnames(x)=='lon']
  D <- dim(lons)
  lons <- matrix(lons,1,D[1]*D[2])
  lats <- matrix(lats,1,D[1]*D[2])
  lons <- round(lons/dlon)*dlon
  lats <- round(lats/dlat)*dlat
  xy <- xy.coords(lons,lats,'lon','lat','')
  tt <- xyTable(xy, digits=digits)
  invisible(tt)  
}


sunflower.storm <- function(x, it = NULL, is = NULL,
             xlab = NULL, ylab = NULL, xlim = NULL, ylim = NULL,
             add = FALSE, rotate = FALSE, leg=TRUE,
             projection='latlon', petalsize = 10, dlon = 5, dlat = 2,
             pch = 1, cex = 0.8, cex.fact =  1.5,
             col = 'red', col.small = 'grey',
             bg = 'white', size = 1/24, seg.lwd = 1.5,
             alpha = 0.5, lonR = 0, latR = 90, ...)
{
  
  y <- subset.storm(x,it=it,is=is)
  tt <- petals.storm(y,dlon=dlon,dlat=dlat)
  
  if(!add & (projection=='latlon' | projection=='lonlat')) {
     xlab <- if (is.null(xlab)) 'lon' else xlab
     ylab <- if (is.null(ylab)) 'lat' else ylab
     xlim <- if (is.null(xlim)) range(tt$x[is.finite(tt$x)]) else xlim
     ylim <- if (is.null(ylim)) range(tt$y[is.finite(tt$y)]) else ylim
  }

  lon <- tt$x
  lat <- tt$y
  number <- tt$number
  number.original <- number
  if (is.numeric(petalsize) & (petalsize>1)) {
    number <- sapply(number,function(x) max(c(1,floor(x/petalsize)))) }

  if (projection=='sphere' | projection=='np' | projection=='sp') {
    A <- sphere.rotate(lon,lat)
    X <- A[1,]; Y <- A[2,]; Z <- A[3,]
    lon <- X[Y>0]; lat <- Z[Y>0]
    xlim <- NULL; ylim <- NULL
    xlab <- ''; ylab <- ''
    par(bty="n",xaxt="n",yaxt="n",new=TRUE)
  } 
  
  if(!add) {
    plot(lon, lat, xlab = xlab, ylab = ylab,
      xlim=xlim, ylim=ylim, type="n",frame.plot=F)
  }

  n <- length(lon)
  n.is1 <- number == 1
  if(any(n.is1))
    points(lon[ n.is1], lat[ n.is1], pch=pch, col=col.small, bg=bg, cex= cex)
  if(any(!n.is1)) {
    points(lon[!n.is1], lat[!n.is1], pch=pch,
           col=adjustcolor(col,alpha.f=alpha),bg=bg,
           cex= cex/cex.fact,)
    i.multi <- (1L:n)[number > 1]
    ppin <- par("pin")
    pusr <- par("usr")
    xr <- size * abs(pusr[2L] - pusr[1L])/ppin[1L]
    yr <- size * abs(pusr[4L] - pusr[3L])/ppin[2L]

    i.rep <- rep.int(i.multi, number[number > 1])
    z <- numeric()
    for(i in i.multi)
       z <- c(z, 1L:number[i] + if(rotate) stats::runif(1) else 0)
    deg <- (2 * pi * z)/number[i.rep]
    segments(lon[i.rep], lat[i.rep],
             lon[i.rep] + xr * sin(deg),
             lat[i.rep] + yr * cos(deg),
             col=adjustcolor(col,alpha.f=alpha), lwd = seg.lwd)
  }

  data("geoborders",envir=environment())
  ok <- is.finite(geoborders$x) & is.finite(geoborders$y)
  mlon <- geoborders$x[ok]
  mlat <- geoborders$y[ok]
  
  if (projection=='sphere' | projection=='np' | projection=='sp') {
    A <- sphere.rotate(mlon,mlat)
    X <- A[1,]; Y <- A[2,]; Z <- A[3,]
    mlon <- X[Y>0]; mlat <- Z[Y>0]
  }  

  points(mlon,mlat,pch=".")
  lines(cos(pi/180*1:360),sin(pi/180*1:360),col="black")

  if (leg) legend("bottomleft",inset=c(0,0),pch=3,lty=0,
    legend=paste('1 petal =',as.character(petalsize),'obs'),
        col=adjustcolor(col,alpha.f=alpha),
        lwd=seg.lwd,bg='white')
}



#library(esd)
#source('as.storm.R')
#source('subset.storm.R')
#source('stormtools.R')
#
# fname="/vol/fou/klima/IMILAST/ERAinterim_0.75_NH_M07_19790101_20091231.txt"
# m07 <- read.fwf(fname,width=c(2,7,4,11,5,3,3,3,8,8,9),
#         col.names=c("Code99","cyclone","timeStep","Date","Year",
#         "Month","Day","Time","Lon","Lat","Pressure"))
# x <- as.storm(m07)
