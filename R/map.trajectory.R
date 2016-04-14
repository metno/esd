## Author 	 Kajsa Parding
## Last update   16.02.2015
## Require 	 geoborders.rda

map.trajectory <- function(x,it=NULL,is=NULL,type="paths",
                           projection="sphere",verbose=TRUE,...) {
  if (verbose) print("map.trajectory")
  stopifnot(is.trajectory(x))
  y <- subset.trajectory(x,it=it,is=is)
  if (type=='paths' | is.null(type)) {
    if (projection=="sphere" | projection=="np" | projection=="sp") {
      if (projection=="np") latR <- 90
      if (projection=="sp") latR <- -90
      sphere.trajectory(y,...)
    } else if (projection=="latlon" | projection=="lonlat") {
      lonlat.trajectory(y,...)
    }
  } else if (type=='density') {
    map.density.trajectory(y,projection=projection,verbose=verbose,...)
  } else if (type=='shapes') {
    map.anomaly.trajectory(y,projection=projection,verbose=verbose,...)
  } else if (type=='colors') {
    segments.trajectory(y,verbose=verbose,...)
  } else print("unkown map type")
}

map.anomaly.trajectory <- function(x,col=NULL,alpha=NULL,
  main=NULL,xlim=NULL,ylim=NULL,lty=1,lwd=1,verbose=FALSE,new=TRUE) {
  if (verbose) print('map.anomaly.trajectory')
  stopifnot(is.trajectory(x))
  if(!('anomaly' %in% attr(x,'aspect'))) x <- anomaly(x)
  if(is.null(alpha)) alpha <- 0.01 + min(10/(dim(x)[1]),0.5)
  if(is.null(col)) col <- "red"
  if(new) dev.new()
  par(bty="n")
  lons <- x[,colnames(x)=='lon']
  lats <- x[,colnames(x)=='lat']
  plot(lons,lats,type='.',cex=1,col=adjustcolor(col,alpha.f=alpha),
       main=main,xlim=xlim,ylim=ylim)
  matlines(t(lons),t(lats),lty=lty,lwd=lwd,
         col=adjustcolor(col,alpha.f=alpha))
}

segments.trajectory <- function(x,param="month",
      xlim=NULL,ylim=NULL,colbar=list(pal='rainbow',rev=FALSE,n=10,
      breaks=NULL,type="p",cex=2,h=0.6, v=1,pos=0.1,show=TRUE),
      show.start=FALSE,show.end=FALSE,show.segment=TRUE,
      alpha=0.1,cex=0.5,lty=1,lwd=3,main=NULL,new=TRUE,projection="lonlat",
      verbose=FALSE,...) {
  if(verbose) print("segments.trajectory")
  x0 <- x
  if(is.null(dim(x0))) {
    dim(x) <- c(1,length(x0))
    colnames(x) <- names(x0)
  }
  lons <- x[,colnames(x)=='lon']
  lats <- x[,colnames(x)=='lat']
  if(is.null(dim(x0))) {
    dim(lons) <- c(1,length(lons))
    dim(lats) <- c(1,length(lons))
  }
  if (is.null(xlim)) xlim <- range(lons)
  if (is.null(ylim)) ylim <- range(lats)
  if(verbose) print(paste('xlim:',paste(round(xlim),collapse=" - "),
                          ', ylim:',paste(round(ylim),collapse=" - ")))
  lab.breaks <- NULL
  if(is.character(param)) {
    if (tolower(param)=="nao") {
      param <- NAO()
    } else if (tolower(param)=="amo") {
      param <- AMO()
    } else if (tolower(param)=="t2m") {
      param <- HadCRUT4()
    } else if (param %in% colnames(x)) {
      if(sum(colnames(x)==param)==1) {
        p <- matrix(rep(x[,colnames(x)==param],sum(colnames(x)=='lon')),dim(lons))
      } else {
        p <- x[,colnames(x)==param]
      }
    } else if (param %in% c("date","year","month","season")) {
      n <- sum(colnames(x)=="lon")
      p <- t(apply(x,1,function(y) strftime(seq(strptime(y[colnames(x)=="start"],"%Y%m%d%H"),
                                   strptime(y[colnames(x)=="end"],"%Y%m%d%H"),
                                   length.out=n),"%Y%m%d%H") ))
      p <- matrix(as.numeric(p),dim(p))
      if(param=="date") {
        p <- matrix(round(p*1E-2),dim(p))
      } else if(param=="year") {
        p <- matrix(year(as.Date(strptime(p,"%Y%m%d%H"))),dim(p))
      } else if(param=="month") {
        p <- matrix(month(as.Date(strptime(p,"%Y%m%d%H"))),dim(p))
        colbar$breaks <- seq(1,13)
        lab.breaks <- c(month.abb,"")
      } else if(param=="season") {
        fn <- function(x) as.numeric(switch(x, "12"="1", "1"="1", "2"="1",
                               "3"="2", "4"="2", "5"="2",
                               "6"="3", "7"="3", "8"="3",
                               "9"="4", "10"="4", "11"="4"))
        p <- matrix(month(as.Date(strptime(p,"%Y%m%d%H"))),dim(p))
        p <- apply(p, c(1,2), fn)
        colbar$breaks <- seq(1,5)
        lab.breaks <- c("djf","mam","jja","son","")
      }
    }
  }
  if(inherits(param,c("zoo","station"))) {
    n <- sum(colnames(x)=="lon")
    d <- t(apply(x,1,function(y) as.numeric(strftime(seq(strptime(y[colnames(x)=="start"],"%Y%m%d%H"),
                                 strptime(y[colnames(x)=="end"],"%Y%m%d%H"),
                                 length.out=n),"%Y%m%d")) ))
    param <- subset(param,it=as.Date(strptime(c(min(d),max(d)),format="%Y%m%d")))
    tp <- as.numeric(strftime(index(param),format="%Y%m%d"))
    if(inherits(param,"month")) {
      d <- round(d*1E-2)*1E2+1
    } else if(inherits(param,"annual")) {
      d <- round(d*1E-4)
      tp <- as.numeric(index(param))
    }
    p <- matrix(rep(NA,length(d)),dim(d))
    for (i in seq(length(tp))) p[d==tp[i]] <- param[tp==tp[i]]
  }
  #if() {
  #  print(paste("unknown input param =",param))
  #  print(paste("possible choices:"))
  #  print(paste(unique(colnames(x),collapse=",")))
  #  print("year, month, yearmon, season")
  #}
   
  lon0 <- lons[,1:(dim(lons)[2]-1)]
  lon1 <- lons[,2:dim(lons)[2]]
  lat0 <- lats[,1:(dim(lats)[2]-1)]
  lat1 <- lats[,2:dim(lats)[2]]
  pcol <- 0.5*(p[,1:(dim(p)[2]-1)] + p[,2:dim(p)[2]])

  colbar <- colbar.ini(p,FUN=NULL,colbar=colbar,verbose=verbose)
  #colbar$col <- adjustcolor(colbar$col,alpha.f=alpha)
  icol <- apply(pcol,2,findInterval,colbar$breaks)
  icol[icol==0] <- 1
  icol[icol>colbar$n] <- colbar$n
  col <- matrix(colbar$col[icol],dim(pcol))

  data("geoborders",envir=environment())
  #ok <- is.finite(geoborders$x) & is.finite(geoborders$y)
  #if(!is.null(xlim)) ok <- ok & geoborders$x >= min(xlim) & geoborders$x <= max(xlim)
  mlon <- geoborders$x#[ok]
  mlat <- geoborders$y#[ok]

  if (new) dev.new(width=8,height=7)
  par0 <- par()
  par(bty="n",fig=c(0,1,0.1,1))
  plot(mlon,mlat,pch=".",col="white",main=main,
    xlab="lon",ylab="lat",xlim=xlim,ylim=ylim)

  OK <- apply(lons,1,function(x) !((max(x)-min(x))>180))
  if(verbose) print(paste(dim(lons)[1],'trajectories,',
                          sum(!OK),'crossing dateline'))
  
  if(show.segment & sum(OK)>0) {
    segments(lon0[OK,],lat0[OK,],lon1[OK,],lat1[OK,],
             col=adjustcolor(col[OK,],alpha.f=alpha),lty=lty,lwd=lwd)
  }
  if(show.start) {
    points(lon0[OK,1],lat0[OK,1],pch=19,cex=cex,col=adjustcolor(col[OK,1],
           alpha.f=alpha))
  }

  if(show.end) {
    points(lon0[OK,dim(lon0)[2]],lat0[OK,dim(lon0)[2]],pch=18,cex=cex,
           col=adjustcolor(col[OK,dim(lon0)[2]],alpha.f=alpha))
  }
  
  # draw coastlines
  #points(mlon,mlat,pch=".",col='grey20',cex=1.4)
  lines(mlon,mlat,lty=1,col='grey40',lwd=1.4)
  
  par(fig=par0$fig,new=TRUE)
  image.plot(breaks=colbar$breaks,lab.breaks=lab.breaks,horizontal = TRUE,
             legend.only = T, zlim = range(colbar$breaks),
             col = adjustcolor(colbar$col,alpha.f=0.8), legend.width = 1,
             axis.args = list(cex.axis = 0.8,
              xaxp=c(range(colbar$breaks),n=colbar$n)),
             border = FALSE)

  # trajectories crossing the dateline plotted in two parts
  ## if (sum(!OK)>0) {
  ##   fn <- function(lon,lat) {
  ##     lon[lon<0] <- lon[lon<0]+360
  ##     xy <- approx(lon,lat,sort(c(lon,180)))
  ##     lon <- xy$x; lat <- xy$y
  ##     lines(lon[lon<=180],lat[lon<=180],
  ##         lty=lty,lwd=lwd,col=adjustcolor(col,alpha.f=alpha))
  ##     lines(lon[lon>=180]-360,lat[lon>=180],
  ##         lty=lty,lwd=lwd,col=adjustcolor(col,alpha.f=alpha))
  ##   }
  ##   if (sum(!OK)==1) {
  ##     fn(lons[!OK,],lats[!OK,])
  ##   } else {
  ##     for (i in 1:sum(!OK)) fn(lons[!OK,][i,],lats[!OK,][i,])
  ##   }
  ## }
  
}

lonlat.trajectory <- function(x,show.start=TRUE,
    xlim=NULL,ylim=NULL,col='blue',alpha=0.05,cex=1,
    lty=1,lwd=2,main=NULL,add=FALSE,new=TRUE,verbose=FALSE,...) {
  if (verbose) print("lonlat.trajectory")
  x0 <- x
  if(is.null(dim(x0))) {
    dim(x) <- c(1,length(x0))
    colnames(x) <- names(x0)
  }
  lons <- x[,colnames(x)=='lon']
  lats <- x[,colnames(x)=='lat']
  if(dim(x)[1]==1) {
    dim(lons) <- c(1,length(lons))
    dim(lats) <- c(1,length(lats))
  }
  if (is.null(xlim)) xlim <- range(lons)
  if (is.null(ylim)) ylim <- range(lats)
  if(verbose) print(paste('xlim',paste(xlim,collapse="-"),
                          ', ylim',paste(ylim,collapse="-")))
  
  data("geoborders",envir=environment())
  #ok <- is.finite(geoborders$x) & is.finite(geoborders$y)
  mlon <- geoborders$x#[ok]
  mlat <- geoborders$y#[ok]

  if(is.null(dev.list())) add <- FALSE
  if(add) new <- FALSE
  
  if(new) dev.new(width=8,height=7)
  par(bty="n")
  if(!add) plot(mlon,mlat,pch=".",col="white",main=main,
                xlab="lon",ylab="lat",xlim=xlim,ylim=ylim)
  
  OK <- apply(lons,1,function(x) !((max(x)-min(x))>180))
  if(verbose) print(paste(dim(lons)[1],'trajectories,',
                          sum(!OK),'crossing dateline'))

  if(sum(OK)>1) {
    matlines(t(lons[OK,]),t(lats[OK,]),lty=lty,lwd=lwd,
           col=adjustcolor(col,alpha.f=alpha))
  } else if(sum(OK)==1) {
    lines(t(lons[OK,]),t(lats[OK,]),lty=lty,lwd=lwd,
           col=adjustcolor(col,alpha.f=alpha))
  }
  #points(lons[OK,],lats[OK,],pch='.',cex=0.5,
  #       col=adjustcolor(col,min(1,alpha+0.1)))
  if(show.start) points(lons[OK,1],lats[OK,1],pch=19,cex=cex,
                        col=adjustcolor(col,alpha.f=alpha))

  # trajectories crossing the dateline plotted in two parts
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
    if (sum(!OK)==1) {
      fn(lons[!OK,],lats[!OK,])
    } else {
      for (i in 1:sum(!OK)) fn(lons[!OK,][i,],lats[!OK,][i,])
    }
  }

  # draw coastlines
  #points(mlon,mlat,pch=".",col='grey20',cex=1.4)
  lines(mlon,mlat,lty=1,col='grey40',lwd=1.4)
  
  # box marking the spatial subset
  slon <- attr(x0,'longitude')
  slat <- attr(x0,'latitude')
    
  if(verbose & !is.null(slon)) print(paste('subset','lon',paste(slon,collapse="-"),
                          'lat',paste(slat,collapse="-")))
  if (any(!is.null(c(slat,slon)))) {
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


sphere.trajectory <- function(x,
    xlim=NULL,ylim=NULL,col='blue',alpha=0.05,cex=0.5,
    lty=1,lwd=2,lonR=0,latR=90,main=NULL,add=FALSE,
    show.start=TRUE,verbose=FALSE,new=TRUE,...) {
  
  x0 <- x
  if(is.null(dim(x0))) {
    dim(x) <- c(1,length(x0))
    colnames(x) <- names(x0)
  }

  ## KMP 2015-12-07: apply xlim and ylim
  is <- NULL
  if (!is.null(xlim)) is$lon <- xlim
  if (!is.null(ylim)) is$lat <- ylim
  x <- subset(x,is=is)

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
  gx <- geoborders$x
  gy <- geoborders$y
  ok <- is.finite(gx) & is.finite(gy)
  if (!is.null(xlim)) ok <- ok & gx>=min(xlim) & gx<=max(xlim)
  if (!is.null(ylim)) ok <- ok & gy>=min(ylim) & gy<=max(ylim)
  a <- sphere.rotate(gx[ok],gy[ok],lonR=lonR,latR=latR)
  x <- a[1,]; y <- a[2,]; z <- a[3,]

  if(is.null(dev.list())) add <- FALSE
  if(add) new <- FALSE
      
  if(new) dev.new()
  par(bty="n",xaxt="n",yaxt="n")
  if(!add) plot(x[y>0],z[y>0],pch=".",type="n",xlab="",ylab="",main=main)

  matlines(X,Z,lty=lty,lwd=lwd,col=adjustcolor(col,alpha.f=alpha))
  #points(X,Z,pch='.',cex=0.8,col=adjustcolor(col,alpha.f=min(1,alpha+0.2)))
  if(show.start) {
    if(is.null(dim(x0))) {
      points(X[1],Z[1],pch=19,cex=cex,col=adjustcolor(col,alpha.f=alpha))
    } else {
      points(X[1,],Z[1,],pch=19,cex=cex,col=adjustcolor(col,alpha.f=alpha))
    }
  }
  points(x[y>0],z[y>0],pch=".",col='grey30')
  lines(cos(pi/180*1:360),sin(pi/180*1:360),col="black")

  # box marking the spatial subset
  slon <- attr(x0,'longitude')
  slat <- attr(x0,'latitude')
  if(verbose) print(paste('subset','lon',paste(slon,collapse="-"),
                          'lat',paste(slat,collapse="-")))
  if (any(!is.null(c(slat,slon)))) {
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


map.density.trajectory <- function(x,dx=4,dy=2,it=NULL,is=NULL,
      colbar=list(pal='precip',rev=TRUE,breaks=NULL,cex=2,h=0.6,v=1),
      projection='sphere',latR=90,lonR=10,gridlines=FALSE,...) {
  stopifnot(inherits(x,c("trajectory","field")))
  x <- subset(x,it=it,is=is)
  if (!inherits(x,"field")) {
    X <- trajectory2field(x,dt='year',dx=dx,dy=dy)
  } else X <- x
  map(X,colbar=colbar,projection=projection,latR=latR,
      lonR=lonR,gridlines=gridlines,...)
}


map.hexbin.trajectory <- function(x,dx=6,dy=2,it=NULL,is=NULL,Nmax=NULL,
          xgrid=NULL,ygrid=NULL,add=FALSE,leg=TRUE,
          xlim=NULL,ylim=NULL,col='red',border='firebrick4',
          colmap='heat.colors',scale.col=TRUE,scale.size=FALSE,
          main=NULL,new=TRUE,verbose=FALSE) {

  x <- subset.trajectory(x,it=it,is=is)
  ilon <- colnames(x)=='lon'
  ilat <- colnames(x)=='lat'
  ilen <- colnames(x)=='n'
  lon <- unlist(apply(x,1,function(x) approx(x[ilon],n=x[ilen])$y))
  lat <- unlist(apply(x,1,function(x) approx(x[ilat],n=x[ilen])$y))
  if (is.null(xlim)) xlim <- range(lon)
  if (is.null(ylim)) ylim <- range(lat)
  data("geoborders",envir=environment())
  mlon <- geoborders$x
  mlat <- geoborders$y
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
  lines(mlon[OK],mlat[OK],lty=1,col='grey20',lwd=1.0)
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

map.sunflower.trajectory <- function(x,it=NULL,is=NULL,
      dx=6,dy=2,petalsize=7,
      xgrid=NULL,ygrid=NULL,leg=TRUE,leg.loc=2,
      xlim=NULL,ylim=NULL,rotate=TRUE,alpha=0.6,
      main=NULL,new=TRUE,verbose=FALSE) {

  x <- subset.trajectory(x,it=it,is=is)
  ilon <- colnames(x)=='lon'
  ilat <- colnames(x)=='lat'
  ilen <- colnames(x)=='n' 
  lon <- unlist(apply(x,1,function(x) approx.lon(x[ilon],n=x[ilen])$y))
  lat <- unlist(apply(x,1,function(x) approx(x[ilat],n=x[ilen])$y))
  if (is.null(xlim)) xlim <- range(lon)
  if (is.null(ylim)) ylim <- range(lat)

  data("geoborders",envir=environment())
  mlon <- geoborders$x
  mlat <- geoborders$y

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
  lines(mlon[OK],mlat[OK],lty=1,col='grey20',lwd=1)

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
  
map.pca.trajectory <- function(X,projection="sphere",lonR=NULL,latR=NULL,
      xlim=NULL,ylim=NULL,main=NULL,m=2,param=c('lon','lat')) {

  stopifnot(!missing(X), inherits(X,"trajectory"))
  if (inherits(X,'pca')) {
    pca <- X; X <- pca2trajectory(pca)
  } else pca <- PCA.trajectory(X,param=param)

  if (any('anomaly' %in% attr(X,'aspect'))) X <- anomaly2trajectory(X)
  
  U <- attr(pca,'pattern')
  V <- coredata(pca)
  W <- attr(pca,'eigenvalues')
  R2 <- round(100*attr(pca,'eigenvalues')^2/attr(pca,'tot.var'),2)

  if (!is.null(m)) m <- min(m,dim(U)[2])
  else m <- sum(R2>=5)
  
  colvec <- c('red3','mediumblue', 'chartreuse3',
              'darkorange','darkturquoise')
  mean.lon <- fnlon(mean)
  if (is.null(latR)) latR <- 90
  if (is.null(lonR)) lonR <- mean.lon(X[,colnames(X)=='lon'])
  map.trajectory(X,projection=projection,lonR=lonR,latR=latR,
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




