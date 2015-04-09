## Author 	 Kajsa Parding
## Last update   16.02.2015
## Require 	 geoborders.rda

map.storm <- function(x,it=NULL,is=NULL,
      projection="sphere",lonR=10,latR=90,
      col='red',colmap='rainbow',alpha=0.3,pfit=FALSE,
      main=NULL,xlim=NULL,ylim=NULL,
      verbose=FALSE,new=TRUE) {

  y <- subset.storm(x,it=it,is=is)
  if (pfit) {
    lats <- t(polyfit.storm(y))
    y[,colnames(y)=='lat'] <- lats
  }

  if (projection=="sphere" | projection=="np" | projection=="sp") {
    if (projection=="np") latR <- 90
    if (projection=="sp") latR <- -90
    sphere.storm(y,new=new,verbose=verbose,
    lonR=lonR,latR=latR,col=col,alpha=alpha,main=main)
  } else if (projection=="latlon" | projection=="lonlat") {
    lonlat.storm(y,new=new,verbose=verbose,
    xlim=xlim,ylim=ylim,col=col,alpha=alpha,main=main)
  }
}


lonlat.storm <- function(x,
    xlim=NULL,ylim=NULL,col='blue',alpha=0.1,
    lty=1,lwd=1,main=NULL,new=TRUE,verbose=FALSE) {
  
  x0 <- x
  lons <- x[,colnames(x)=='lon']
  lats <- x[,colnames(x)=='lat']
  if (is.null(xlim)) xlim <- range(lons)
  if (is.null(ylim)) ylim <- range(lats)
  if(verbose) print(paste('xlim',paste(xlim,collapse="-"),
                          ', ylim',paste(ylim,collapse="-")))

  data("geoborders",envir=environment())
  ok <- is.finite(geoborders$x) & is.finite(geoborders$y)
  mlon <- geoborders$x[ok]
  mlat <- geoborders$y[ok]

  if (new) dev.new(width=8,height=7)
  par(bty="n")
  plot(mlon,mlat,pch=".",col="white",main=main,
    xlab="lon",ylab="lat",xlim=xlim,ylim=ylim)

  OK <- apply(lons,1,function(x) !((max(x)-min(x))>180))
  if(verbose) print(paste(dim(lons)[1],'storms,',sum(!OK),'crossing dateline'))
  matlines(t(lons[OK,]),t(lats[OK,]),lty=lty,lwd=lwd,
           col=adjustcolor(col,alpha.f=alpha))

  # storms crossing the dateline plotted in two parts
  if (sum(!OK)>0) {
    fn <- function(lon,lat) {
      lon[lon<0] <- lon[lon<0]+360
      xy <- approx(lon,lat,sort(c(lon,180)))
      lon <- xy$x; lat <- xy$y
      lines(lon[lon<=180],lat[lon<=180],
          lty=lty,lwd=lwd,col=adjustcolor(col,alpha.f=alpha))
      lines(lon[lon>=180]-360,lat[lon>=180],
          lty=lty,lwd=lwd,col=adjustcolor(col,alpha.f=alpha))
    }
    for (i in 1:sum(!OK)) fn(lons[!OK,][i,],lats[!OK,][i,])
  }

  # draw coastlines
  points(mlon,mlat,pch=".",col='grey20',cex=1.4)
  
  # box marking the spatial subset
  slon <- attr(x0,'longitude')
  slat <- attr(x0,'latitude')
  if(verbose) print(paste('subset','lon',paste(slon,collapse="-"),
                          'lat',paste(slat,collapse="-")))
  if (any(!is.na(c(slat,slon)))) {
    if(verbose) print('draw subset box')
    if (sum(is.na(attr(x0,'longitude')))==0) {
      xlim <- attr(x0,'longitude')
    } else {
      xlim <- c(min(x0[,colnames(x0)=='lon']),
                max(x0[,colnames(x0)=='lon']))
    }
    if (sum(is.na(attr(x0,'latitude')))==0) {
      ylim <- attr(x0,'latitude')
    } else {
      ylim <- c(min(x0[,colnames(x0)=='lat']),
                max(x0[,colnames(x0)=='lat']))
    }
    if(verbose) print(paste('xlim',paste(xlim,collapse="-"),
                            'ylim',paste(ylim,collapse="-")))
    xbox <- c(xlim[1],xlim[2],xlim[2],xlim[1],xlim[1])
    ybox <- c(ylim[1],ylim[1],ylim[2],ylim[2],ylim[1])
    lines(xbox,ybox,lty=1,col='grey20',lwd=1.0)
  }
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
    lty=1,lwd=1,lonR=0,latR=90,main=NULL,
    verbose=FALSE,new=TRUE) {
  
  x0 <- x
  ilons <- colnames(x)=='lon'
  ilats <- colnames(x)=='lat'
  lon <- x[,ilons]
  lon[lon<0] <- lon[lon<0]+360
  x[,ilons] <- lon
  
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
  par(bty="n",xaxt="n",yaxt="n")
  plot(x[y>0],z[y>0],pch=".",type="n",xlab="",ylab="",main=main)

  matlines(X,Z,lty=lty,lwd=lwd,col=adjustcolor(col,alpha.f=alpha))
  
  points(x[y>0],z[y>0],pch=".",col='grey30')
  lines(cos(pi/180*1:360),sin(pi/180*1:360),col="black")

  # box marking the spatial subset
  slon <- attr(x0,'longitude')
  slat <- attr(x0,'latitude')
  if(verbose) print(paste('subset','lon',paste(slon,collapse="-"),
                          'lat',paste(slat,collapse="-")))
  if (any(!is.na(c(slat,slon)))) {
    if(verbose) print('draw subset box')
    if (sum(is.na(attr(x0,'longitude')))==0) {
      xlim <- attr(x0,'longitude')
    } else {
      xlim <- c(min(x0[,colnames(x0)=='lon']),
                max(x0[,colnames(x0)=='lon']))
    }
    if (sum(is.na(attr(x0,'latitude')))==0) {
      ylim <- attr(x0,'latitude')
    } else {
      ylim <- c(min(x0[,colnames(x0)=='lat']),
                max(x0[,colnames(x0)=='lat']))
    }
    if(verbose) print(paste('xlim',paste(xlim,collapse="-"),
                            'ylim',paste(ylim,collapse="-")))
    xbox <- c(xlim[1],xlim[2],xlim[2],xlim[1],xlim[1])
    ybox <- c(ylim[1],ylim[1],ylim[2],ylim[2],ylim[1])
    xbox <- approx(xbox,n=200)$y
    ybox <- approx(ybox,n=200)$y
    a <- sphere.rotate(xbox,ybox,lonR=lonR,latR=latR)
    x <- a[1,]; y <- a[2,]; z <- a[3,]
    lines(x,z,lty=1,col='grey20',lwd=1.0)
  }
}


map.hexbin.storm <- function(x,dx=6,dy=2,it=NULL,is=NULL,Nmax=NULL,
          xgrid=NULL,ygrid=NULL,add=FALSE,leg=TRUE,
          xlim=NULL,ylim=NULL,col='red',border='firebrick4',
          colmap='heat.colors',scale.col=TRUE,scale.size=FALSE,
          main=NULL,new=TRUE) {

  x <- subset.storm(x,it=it,is=is)
  ilon <- colnames(x)=='lon'
  ilat <- colnames(x)=='lat'
  ilen <- colnames(x)=='n'
  lon <- unlist(apply(x,1,function(x) approx(x[ilon],n=x[ilen])$y))
  lat <- unlist(apply(x,1,function(x) approx(x[ilat],n=x[ilen])$y))
  if (is.null(xlim)) xlim <- range(lon)
  if (is.null(ylim)) ylim <- range(lat)
  data("geoborders",envir=environment())
  ok <- is.finite(geoborders$x) & is.finite(geoborders$y)
  mlon <- geoborders$x[ok]
  mlat <- geoborders$y[ok]
  if(new) dev.new(width=8,height=7)
  if(leg) par(bty="n",mar=c(5.0,4.0,3.0,5.3))
  else par(bty="n",mar=c(4.4,4.0,1.0,1.0))
  if(!add) plot(lon, lat, xlab="lon", ylab="lat", main=main,
                xlim=xlim,ylim=ylim,type="n",frame.plot=F)
  OK <- (findInterval(lon,xlim)==1 & findInterval(lat,ylim)==1)
  scatter.hexbin(lon[OK],lat[OK],dx=dx,dy=dy,xgrid=xgrid,ygrid=ygrid,
                 new=FALSE,leg=leg,col=col,border=border,Nmax=Nmax,
                 scale.col=scale.col,scale.size=scale.size,colmap=colmap)
  OK <- (findInterval(mlon,xlim)==1 & findInterval(mlat,ylim)==1)
  points(mlon[OK],mlat[OK],pch=".",col='grey20',cex=1.4)
  # box marking the spatial subset
  slon <- attr(x,'longitude')
  slat <- attr(x,'latitude')
  if(verbose) print(paste('subset','lon',paste(slon,collapse="-"),
                          'lat',paste(slat,collapse="-")))
  if (any(!is.na(c(slat,slon)))) {
    if(verbose) print('draw subset box')
    if (sum(is.na(attr(x,'longitude')))==0) {
      xlim <- attr(x,'longitude')
    } else {
      xlim <- c(min(x[,colnames(x)=='lon']),
                max(x[,colnames(x)=='lon']))
    }
    if (sum(is.na(attr(x,'latitude')))==0) {
      ylim <- attr(x,'latitude')
    } else {
      ylim <- c(min(x[,colnames(x)=='lat']),
                max(x[,colnames(x)=='lat']))
    }
    if(verbose) print(paste('xlim',paste(xlim,collapse="-"),
                            'ylim',paste(ylim,collapse="-")))
    xbox <- c(xlim[1],xlim[2],xlim[2],xlim[1],xlim[1])
    ybox <- c(ylim[1],ylim[1],ylim[2],ylim[2],ylim[1])
    lines(xbox,ybox,lty=1,col='grey20',lwd=1.0)
  }
}

map.sunflower.storm <- function(x,it=NULL,is=NULL,
      dx=6,dy=2,petalsize=7,
      xgrid=NULL,ygrid=NULL,leg=TRUE,leg.loc=2,
      xlim=NULL,ylim=NULL,rotate=TRUE,alpha=0.6,
      main=NULL,new=TRUE) {

  x <- subset.storm(x,it=it,is=is)
  ilon <- colnames(x)=='lon'
  ilat <- colnames(x)=='lat'
  ilen <- colnames(x)=='n' 
  lon <- unlist(apply(x,1,function(x) approx.lon(x[ilon],n=x[ilen])$y))
  lat <- unlist(apply(x,1,function(x) approx(x[ilat],n=x[ilen])$y))
  if (is.null(xlim)) xlim <- range(lon)
  if (is.null(ylim)) ylim <- range(lat)

  data("geoborders",envir=environment())
  ok <- is.finite(geoborders$x) & is.finite(geoborders$y)
  mlon <- geoborders$x[ok]
  mlat <- geoborders$y[ok]

  if(new) dev.new(width=8,height=7)
  par(bty="n",mar=c(4.4,4.0,1.0,1.0))
  OK <- (findInterval(lon,xlim)==1 & findInterval(lat,ylim)==1)
  scatter.sunflower(lon[OK],lat[OK],petalsize=petalsize,
           dx=dx,dy=dy,xlab='lon',yla='lat',
           xgrid=xgrid,ygrid=ygrid,leg=leg,leg.loc=leg.loc,
           xlim=xlim,ylim=ylim,rotate=rotate,alpha=alpha,
           main=main,new=FALSE)

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
  points(mlon[OK],mlat[OK],pch=".",col='grey20',cex=1.4)

  # box marking the spatial subset
  slon <- attr(x,'longitude')
  slat <- attr(x,'latitude')
  if(verbose) print(paste('subset','lon',paste(slon,collapse="-"),
                          'lat',paste(slat,collapse="-")))
  if (any(!is.na(c(slat,slon)))) {
    if(verbose) print('draw subset box')
    if (sum(is.na(attr(x,'longitude')))==0) {
      xlim <- attr(x,'longitude')
    } else {
      xlim <- c(min(x[,colnames(x)=='lon']),
                max(x[,colnames(x)=='lon']))
    }
    if (sum(is.na(attr(x,'latitude')))==0) {
      ylim <- attr(x,'latitude')
    } else {
      ylim <- c(min(x[,colnames(x)=='lat']),
                max(x[,colnames(x)=='lat']))
    }
    if(verbose) print(paste('xlim',paste(xlim,collapse="-"),
                            'ylim',paste(ylim,collapse="-")))
    xbox <- c(xlim[1],xlim[2],xlim[2],xlim[1],xlim[1])
    ybox <- c(ylim[1],ylim[1],ylim[2],ylim[2],ylim[1])
    lines(xbox,ybox,lty=1,col='grey20',lwd=1.0)
  }
}


mean.lon <- function(lon) {
  if (!any(lon>0)|!any(lon<0)|(mean(lon[lon>0])-mean(lon[lon<0]))<120){
    x <- mean(lon)
  } else {
    lon[lon<0] <- lon[lon<0]+360
    x <- mean(lon)
    if (x>180) x <- x-360
  }
  return(x)
}
  
map.pca.storm <- function(X,projection="sphere",lonR=NULL,latR=NULL,
      xlim=NULL,ylim=NULL,main=NULL,m=2,param=c('lon','lat')) {

  stopifnot(!missing(X), inherits(X,"storm"))
  if (inherits(X,'pca')) {
    pca <- X; X <- pca2storm(pca)
  } else pca <- PCA.storm(X,param=param)

  U <- attr(pca,'pattern')
  V <- coredata(pca)
  W <- attr(pca,'eigenvalues')
  R2 <- round(100*attr(pca,'eigenvalues')^2/attr(pca,'tot.var'),2)

  if (!is.null(m)) m <- min(m,dim(U)[2])
  else m <- sum(R2>=5)
  
  colvec <- c('red3','mediumblue', 'chartreuse3',
              'darkorange','darkturquoise')

  if (is.null(latR)) latR <- 90
  if (is.null(lonR)) lonR <- mean.lon(X[,colnames(X)=='lon'])
  map.storm(X,projection=projection,lonR=lonR,latR=latR,
    col='grey20',alpha=0.1,xlim=xlim,ylim=ylim,main=main,new=TRUE)
 
  for (i in 1:m) { 
    X.PC.max <- max(V[,i]) * (U[,i]*W[i])
    X.PC.min <- min(V[,i]) * (U[,i]*W[i])
    if (any(aspect(pca)=='anomaly')) {
      for (j in 1:length(attr(pca,'mean'))) {
        if ((names(attr(pca,'mean'))[j])=='lon') {
          mj <- mean.lon(unlist(attr(pca,'mean')[j]))
        } else {
          mj <- mean(unlist(attr(pca,'mean')[j]))
        }
        X.PC.max[attr(pca,'colnames')==names(attr(pca,'mean'))[j]] <-
          X.PC.max[attr(pca,'colnames')==names(attr(pca,'mean'))[j]] + mj
        X.PC.min[attr(pca,'colnames')==names(attr(pca,'mean'))[j]] <-
          X.PC.min[attr(pca,'colnames')==names(attr(pca,'mean'))[j]] + mj
      }
    }

    lon.max <- X.PC.max[attr(pca,'colnames')=='lon']
    lat.max <- X.PC.max[attr(pca,'colnames')=='lat']
    lon.min <- X.PC.min[attr(pca,'colnames')=='lon']
    lat.min <- X.PC.min[attr(pca,'colnames')=='lat']
    if (any(projection %in% c('sphere','np','sp'))) {
      # rotate lon.max and lat.max
      a <- sphere.rotate(lon.max,lat.max,lonR=lonR,latR=latR)
      x <- a[1,]; y <- a[2,]; z <- a[3,]
      lon.max <- x[y>0]; lat.max <- z[y>0]
      # rotate lon.min and lat.min
      a <- sphere.rotate(lon.min,lat.min,lonR=lonR,latR=latR)
      x <- a[1,]; y <- a[2,]; z <- a[3,]
      lon.min <- x[y>0]; lat.min <- z[y>0] 
    }
    points(lon.max,lat.max,col=colvec[i],type='l',lwd=2,lty=1)#,pch=19)
    points(lon.max[1],lat.max[1],col=colvec[i],type='p',pch=19)
    points(lon.min,lat.min,col=colvec[i],type='l',lwd=2,lty=2)
    points(lon.min[1],lat.min[1],col=colvec[i],type='p',pch=19)
  }

  invisible(pca)
}




