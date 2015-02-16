## Author 	 Kajsa Parding
## Last update   16.02.2015
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
  par(bty="n",xaxt="n",yaxt="n")
  plot(mlon,mlat,pch=".",col="white",
    xlab="lon",ylab="lat",xlim=xlim,ylim=ylim)

  OK <- apply(lons,1,function(x) !((max(x)-min(x))>180))
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

map.hexbin.storm <- function(x,dx=6,dy=2,Nmax=NULL,
      xgrid=NULL,ygrid=NULL,add=FALSE,leg=TRUE,
      xlim=NULL,ylim=NULL,col='red',border='firebrick4') {

  lon <- x[,colnames(x)=='lon']
  lat <- x[,colnames(x)=='lat']
  lon <- matrix(lon,length(lon),1)
  lat <- matrix(lat,length(lat),1)
  if (is.null(xlim)) xlim <- range(lon)
  if (is.null(ylim)) ylim <- range(lat)

  data("geoborders",envir=environment())
  ok <- is.finite(geoborders$x) & is.finite(geoborders$y)
  mlon <- geoborders$x[ok]
  mlat <- geoborders$y[ok]

  if(!add) {
    if(leg) par(bty="n",mar=c(5.0,4.0,3.0,5.3))
    else par(bty="n",mar=c(4.4,4.0,1.0,1.0))
    plot(lon, lat, xlab="lon", ylab="lat",
         xlim=xlim,ylim=ylim,type="n",frame.plot=F)
  }
  
  OK <- (findInterval(lon,xlim)==1 & findInterval(lat,ylim)==1)
  binscatter.hex(lon[OK],lat[OK],dx=dx,dy=dy,xgrid=xgrid,ygrid=ygrid,
                 new=FALSE,leg=leg,col=col,border=border,Nmax=Nmax)

  OK <- (findInterval(mlon,xlim)==1 & findInterval(mlat,ylim)==1)
  points(mlon[OK],mlat[OK],pch=".",col='grey60')
}

map.sunflower.storm <- function(x,dx=6,dy=2,petalsize=7,
      xgrid=NULL,ygrid=NULL,leg=TRUE,leg.loc=2,
      xlim=NULL,ylim=NULL,rotate=TRUE,alpha=0.6) {

  lon <- x[,colnames(x)=='lon']
  lat <- x[,colnames(x)=='lat']
  lon <- matrix(lon,length(lon),1)
  lat <- matrix(lat,length(lat),1)
  if (is.null(xlim)) xlim <- range(lon)
  if (is.null(ylim)) ylim <- range(lat)

  data("geoborders",envir=environment())
  ok <- is.finite(geoborders$x) & is.finite(geoborders$y)
  mlon <- geoborders$x[ok]
  mlat <- geoborders$y[ok]

  par(bty="n",mar=c(4.4,4.0,1.0,1.0))
  OK <- (findInterval(lon,xlim)==1 & findInterval(lat,ylim)==1)
  binscatter.sunflower(lon[OK],lat[OK],petalsize=petalsize,
              dx=dx,dy=dy,xlab='lon',yla='lat',
              xgrid=xgrid,ygrid=ygrid,leg=leg,leg.loc=leg.loc,
              xlim=xlim,ylim=ylim,rotate=rotate,alpha=alpha)

  OK <- (findInterval(mlon,xlim)==1 & findInterval(mlat,ylim)==1)
  if (leg) {
    if (leg.loc==1) {
      xbox <- c(max(xlim)-0.35*(max(xlim)-min(xlim)),max(xlim))
      ybox <- c(max(ylim)-0.1*(max(ylim)-min(ylim)),max(ylim))
    } else if (leg.loc==2 | is.null(leg.loc)) {
      xbox <- c(min(xlim),min(xlim)+0.35*(max(xlim)-min(xlim)))
      ybox <- c(max(ylim)-0.1*(max(ylim)-min(ylim)),max(ylim))
    } else if (leg.loc==3) {
      xbox <- c(min(xlim),min(xlim)+0.35*(max(xlim)-min(xlim)))
      ybox <- c(min(ylim),min(ylim)+0.1*(max(ylim)-min(ylim)))
    } else if (leg.loc==4) {
      xbox <- c(max(xlim)-0.35*(max(xlim)-min(xlim)),max(xlim))
      ybox <- c(min(ylim),min(ylim)+0.1*(max(ylim)-min(ylim)))
    }
    OK <- OK & !(findInterval(mlon,xbox)==1 & findInterval(mlat,ybox)==1)
  }
  points(mlon[OK],mlat[OK],pch=".",col='grey20')
}

