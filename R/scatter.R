## Scatter

# Binned scatterplot with sunflowers

scatter <- function(x,y,type='heat', verbose=FALSE,...) {
  if (verbose) print(match.call())
  if (tolower(type)=='heat') scatter.heat(x,y,verbose=verbose,...)
  if (tolower(type)=='sunflower') scatter.sunflower(x,y,verbose=verbose,...)
  if (tolower(type)=='hexbin') scatter.hexbin(x,y,verbose=verbose,...)
}

scatter.heat <- function(x,y,xlim=NULL,ylim=NULL,breaks=NULL,main='Scatter',xlab='',ylab='',sub='',
                         col=NULL,log=FALSE,dig=NULL, fig = c(0.65,0.85,0.22,0.32), verbose=FALSE) {
  if (verbose) print('scatter.heat')
  par(bty='n',mar=c(5.1, 5.1, 4.1, 2.1))
  if (is.null(dig)) dig <- max(c(-2*log(max(c(x,y,na.rm=TRUE)))/log(10),0),na.rm=TRUE)
  if (verbose) print(paste(dig,'digits. Max:',max(c(x,y,na.rm=TRUE))))
  txy <- table(round(x,dig),round(y,dig))
  if (log) txy <- log(txy)/log(10)
  if (is.null(breaks) & is.null(col)) {
    breaks <- pretty(as.numeric(txy))
    col <- heat.colors(length(breaks)-1)
  }
  if (is.null(breaks)) {
    breaks <- pretty(as.numeric(txy),n=length(col)-1)
    if (length(breaks) != length(col)-1)
    breaks <- seq(min(as.numeric(txy),na.rm=TRUE), max(as.numeric(txy),na.rm=TRUE),length=length(col)+1)
  }
  if (is.null(col)) col <- heat.colors(length(breaks) - 1)
  if (verbose) print(c(length(col),length(breaks)))
  image(as.numeric(rownames(txy)),as.numeric(colnames(txy)),txy,breaks=breaks,
        main=main,xlab=xlab,ylab=ylab,sub=sub,col=col)
  grid()
  lines(c(0,0.6),c(0,0.6),lty=2)
  if (!is.null(fig)) { 
    par(fig=fig,yaxt='n',yaxt='n',mar=c(2,0,0,0),new=TRUE,cex.axis=0.75,cex.lab=0.75,col.axis='white')
    image(breaks,1:2,cbind(c(breaks),c(breaks)),breaks=breaks,col=heat.colors(length(breaks)-1))
    par(xaxt='s',col.axis='grey')
    if (log) axis(1,at=breaks[seq(1,length(breaks),by=2)],labels=paste0('10^',
                                                                        breaks[seq(1,length(breaks),by=2)]),cex=0.75) else
                                                                          axis(1,at=breaks[seq(1,length(breaks),by=2)],labels=
                                                                                 as.character(breaks[seq(1,length(breaks),by=2)]),cex=0.75)       
    par(fig=c(0,1,0,1),yaxt='s',yaxt='s',c(5.1, 5.1, 4.1, 2.1),col.axis='black',new=FALSE)
  }
}

scatter.sunflower <- function(x,y,petalsize=7,dx=NULL,dy=NULL,
                              xgrid=NULL,ygrid=NULL,xlim=NULL,ylim=NULL,
                              xlab=NULL,ylab=NULL,main=NULL,leg=TRUE,rotate=TRUE,
                              alpha=0.6,leg.loc=2,new=TRUE, verbose=FALSE) {
  
  stopifnot(is.numeric(x) & is.numeric(y) & length(x)==length(y))
  i <- !(is.na(x) | is.na(y))
  x <- x[i]; y <- y[i]
  
  # Define grid
  if (is.null(dx) & length(xgrid)<=2) dx <- (max(x)-min(x))/20
  if (is.null(dy) & length(ygrid)<=2) dy <- (max(y)-min(y))/20
  
  if (is.null(xgrid)) {
    xgrid <- seq(min(x),max(x),dx*(1+sin(pi/6)))
  } else if (length(xgrid)==2) {
    xgrid <- seq(min(xgrid),max(xgrid),dx*(1+sin(pi/6)))
  } else if (length(xgrid)>2) {
    dx <- mean(xgrid[2:length(xgrid)]-xgrid[1:(length(xgrid)-1)])/(1+sin(pi/6))
  }
  
  if (is.null(ygrid)) {
    ygrid <- seq(min(y),max(y),2*dy*cos(pi/6))
  } else if (length(ygrid)==2) {
    ygrid <- seq(min(ygrid),max(ygrid),2*dy*cos(pi/6))
  } else if (length(ygrid)>2) {
    dy <- mean(ygrid[2:length(ygrid)]-ygrid[1:(length(ygrid)-1)])/(2*cos(pi/6))
  }
  
  Y <- replicate(length(xgrid),ygrid)
  X <- t(replicate(length(ygrid),xgrid))
  fn <- function(x) {
    dx <- x[2:length(x)]-x[1:(length(x)-1)]
    x[1:(length(x)-1)] <- x[1:(length(x)-1)]+dx/2
    x[length(x)] <- x[length(x)]+dx[length(dx)]/2
    return(x)
  }
  Y[,seq(2,dim(Y)[2],2)] <- apply(Y[,seq(2,dim(Y)[2],2)],2,fn)
  
  # Count observations in each grid point
  XYN <- bin(x,y,X,Y)
  X <- XYN[,1]; Y <- XYN[,2]; N <- XYN[,3]
  
  # Define stuff for plot
  dx <- 0.9*dx; dy <- dy*0.9
  if (is.null(xlim)) xlim <- c(min(x)-dx/2,max(x)+dx/2)
  if (is.null(ylim)) ylim <- c(min(y)-dy/2,max(y)+dy/2)
  xr <- 0.8*dx
  yr <- 0.8*dy
  n <- length(X)
  
  # Generate figure
  if(new) dev.new()
  plot(X,Y,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main,type='n')
  
  # Grid points with 1 observation
  if (any(N==1)) symbols(X[N==1],Y[N==1],
                         circles=rep(1,sum(N==1)),fg='blue',bg=F,
                         inches=xr/(max(xlim)-min(xlim)),add=T)
  
  # Grid points with few observations
  if (any(N>1) & any(N<petalsize*2)) {
    i.multi <- (1L:n)[( N>1 & N<petalsize*2 )]
    # Plot hexagons
    mapply(polygon.fill,X[i.multi],Y[i.multi],dx,dy,
           col=adjustcolor('khaki1',alpha.f=alpha),
           border=adjustcolor('khaki3',alpha.f=alpha),n=6)
    # Draw sunflowers
    i.rep <- rep.int(i.multi, N[i.multi])
    z <- numeric()
    for(i in i.multi)
      z <- c(z, 1L:N[i] + if(rotate) stats::runif(1) else 0)
    deg <- (2 * pi * z)/N[i.rep]
    segments(X[i.rep], Y[i.rep],
             X[i.rep] + xr * sin(deg),
             Y[i.rep] + yr * cos(deg),
             col="orange", lwd=1.5)
  }
  
  # Grid points with many observations
  if (any(N>=(petalsize*2))) {
    N2 <- floor(N/petalsize)
    i.multi <- (1L:n)[( N2>1 )]
    # Plot hexagons
    mapply(polygon.fill,X[i.multi],Y[i.multi],dx,dy,
           col=adjustcolor('coral',alpha.f=alpha),
           border=adjustcolor('coral2',alpha.f=alpha),n=6)
    # Draw sunflowers
    i.rep <- rep.int(i.multi, N2[i.multi])
    z <- numeric()
    for(i in i.multi)
      z <- c(z, 1L:N2[i] + if(rotate) stats::runif(1) else 0)
    deg <- (2 * pi * z)/N2[i.rep]
    segments(X[i.rep], Y[i.rep],
             X[i.rep] + xr * sin(deg),
             Y[i.rep] + yr * cos(deg),
             col="tomato4", lwd=1.5)
  }
  
  # Legend
  if (leg) {
    
    dx <- (max(xlim)-min(xlim))/20
    dy <- (max(ylim)-min(ylim))/20
    
    if (any( leg.loc %in% c(1,'upper right'))) {
      xy <- polygon.vertex(
        max(xlim)-3.2*dx,max(ylim)-dy/2,5.3*dx,dy*1.5,n=4,rot=pi/4)
    } else if (any( leg.loc %in% c(2,'upper left',NULL))) {
      xy <- polygon.vertex(
        min(xlim)+3.2*dx,max(ylim)-dy/2,5.3*dx,dy*1.5,n=4,rot=pi/4)
    } else if (any( leg.loc %in% c(3,'lower left'))) {
      xy <- polygon.vertex(
        min(xlim)+3.2*dx,min(ylim)+dy/2,5.3*dx,dy*1.5,n=4,rot=pi/4)
    } else if (any( leg.loc %in% c(4,'lower right'))) {
      xy <- polygon.vertex(
        max(xlim)-3.2*dx,min(ylim)+dy/2,5.3*dx,dy*1.5,n=4,rot=pi/4)
    }
    
    polygon(xy[,1],xy[,2],col='white',border='gray')
    polygon.fill(xy[2,1]+dx,xy[2,2]-dy*0.6,dx*0.35,dy*0.35,
                 col=adjustcolor('khaki1',alpha.f=alpha),
                 border=adjustcolor('khaki3',alpha.f=alpha),n=6)
    polygon.fill(xy[3,1]+dx,xy[3,2]+dy*0.6,dx*0.35,dy*0.35,
                 col=adjustcolor('coral',alpha.f=alpha),
                 border=adjustcolor('coral2',alpha.f=alpha),n=6)
    points(xy[2,1]+dx,xy[2,2]-dy*0.6,pch=3,col="orange",lwd=1)
    points(xy[3,1]+dx,xy[3,2]+dy*0.6,pch=3,col="tomato4",lwd=1)
    text(xy[2,1]+dx*1.5,xy[2,2]-dy*0.7,'1 petal = 1 obs',pos=4)
    text(xy[3,1]+dx*1.5,xy[3,2]+dy*0.5,paste('1 petal = ',
                                             as.character(petalsize),' obs'),pos=4)
  }
}

# Binned scatterplot with hexagons
scatter.hexbin <- function(x,y,new=TRUE,Nmax=NULL,
                           dx=NULL,dy=NULL,xgrid=NULL,ygrid=NULL,
                           xlim=NULL,ylim=NULL,
                           xlab=NULL,ylab=NULL,main=NULL,
                           leg=TRUE,col='blue',border='white',
                           colmap='gray.colors',
                           scale.col=TRUE,scale.size=FALSE, verbose=FALSE) {
  
  stopifnot(is.numeric(x) & is.numeric(y) & length(x)==length(y))
  i <- !(is.na(x) | is.na(y))
  x <- x[i]; y <- y[i]
  
  # Define grid
  if (is.null(dx) & length(xgrid)<=2) dx <- (max(x)-min(x))/20
  if (is.null(dy) & length(ygrid)<=2) dy <- (max(y)-min(y))/20
  
  if (is.null(xgrid)) {
    xgrid <- seq(min(x),max(x),dx*(1+sin(pi/6)))
  } else if (length(xgrid)==2) {
    xgrid <- seq(min(xgrid),max(xgrid),dx*(1+sin(pi/6)))
  } else if (length(xgrid)>2) {
    dx <- mean(xgrid[2:length(xgrid)]-xgrid[1:(length(xgrid)-1)])/(1+sin(pi/6))
  }
  
  if (is.null(ygrid)) {
    ygrid <- seq(min(y),max(y),2*dy*cos(pi/6))
  } else if (length(ygrid)==2) {
    ygrid <- seq(min(ygrid),max(ygrid),2*dy*cos(pi/6))
  } else if (length(ygrid)>2) {
    dy <- mean(ygrid[2:length(ygrid)]-ygrid[1:(length(ygrid)-1)])/(2*cos(pi/6))
  }
  
  Y <- replicate(length(xgrid),ygrid)
  X <- t(replicate(length(ygrid),xgrid))
  fn <- function(x) {
    dx <- x[2:length(x)]-x[1:(length(x)-1)]
    x[1:(length(x)-1)] <- x[1:(length(x)-1)]+dx/2
    x[length(x)] <- x[length(x)]+dx[length(dx)]/2
    return(x)
  }
  Y[,seq(2,dim(Y)[2],2)] <- apply(Y[,seq(2,dim(Y)[2],2)],2,fn)
  
  # Count observations in each grid point
  XYN <- bin(x,y,X,Y)
  X <- XYN[,1]; Y <- XYN[,2]; N <- XYN[,3]
  if(is.null(Nmax)) Nmax <- max(N)
  Nf <- sapply(N/Nmax,function(x) 0.2+0.8*min(1,x))
  X <- X[N>0]; Y <- Y[N>0]
  Nf <- Nf[N>0]; N <- N[N>0]
  
  # Plot
  if (scale.col) {
    colorscale <- colscal(n=10,colmap)[seq(10,1,-1)]
    col <- colorscale[round(Nf*10)]
    border <- colorscale[round(Nf*10)]
  }
  if (scale.size) {
    dX <- dx*Nf
    dY <- dy*Nf
  } else {
    dX <- dx
    dY <- dy
  }
  if (is.null(xlim)) xlim <- c(min(X)-dx,max(X)+dx)
  if (is.null(ylim)) ylim <- c(min(Y)-dy,max(Y)+dy)
  par(bty='n',xpd=NA,mar=c(5.1,4.1,4.1,5.4))
  if(new) plot(x,y,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main,type='n')
  mapply(polygon.fill,X,Y,dX,dY,n=6,col=col,border=border)
  
  # Legend
  if (leg) {
    dn <- round(Nmax/7)
    if (dn>5 & dn<50) { dn <- round(dn/5)*5 
    } else if (dn>50) { dn <- signif(dn,digits=(nchar(dn)-1))}
    nticks <- seq(Nmax,1,-dn)
    if (length(dX)>1) {
      szticks <- sapply(nticks/Nmax,function(x) 0.2+0.8*min(1,x))
    } else {
      szticks <- rep(1,length(nticks))
    }
    if (length(col)>1) {
      cticks <- colorscale[round(nticks/Nmax*10)]
      bticks <- colorscale[round(nticks/Nmax*10)]
    } else {
      cticks <- rep(col,length(nticks))
      bticks <- rep(border,length(nticks))  
    }              
    x0 <- max(xlim)
    y0 <- max(ylim)
    dy0 <- max(szticks)*dy*2.1
    dx0 <- (max(xlim)-min(xlim))/10+max(szticks)*dy/2
    text(x0+dx0,y0+dy0/2+(max(ylim)-min(ylim))/30,"count")
    polygon.fill(x0+dx0,y0,dx*szticks[1],dy*szticks[1],
                 n=6,col=cticks[1],border=bticks[1])
    if (max(N)>Nmax) {
      text(x0+dx0+max(szticks)*dx,y0,
           paste("\u2265",as.character(nticks[1])),pos=4)
    } else {
      text(x0+dx0+max(szticks)*dx,y0,as.character(nticks[1]),pos=4)
    }
    
    ivec <- 2:10
    j <- 1
    for (i in ivec) {
      polygon.fill(x0+dx0,y0-j*dy0,
                   dx*szticks[i],dy*szticks[i],
                   n=6,col=cticks[i],border=bticks[i])
      text(x0+dx0+max(szticks)*dx,y0-j*dy0,
           as.character(nticks[i]),pos=4)
      j <- j+1
    }
  }
}

# Count observations (x,y) in grid points (X,Y)
bin <- function(x,y,X,Y) {
  fn <- function(x,y) {
    d <- sqrt( (X-x)**2 + (Y-y)**2 )
    imin <- which(d==min(d),arr.ind=T)
    if (length(imin)>2) imin <- imin[1,]
    invisible(imin)
  }
  ibin <- mapply(fn,x,y)
  N <- matrix(rep(0,nrow(X)*ncol(X)),nrow(X),ncol(X))
  for (k in 1:(dim(ibin)[2])) {
    i <- ibin[1,k]
    j <- ibin[2,k]
    N[i,j] <- N[i,j]+1
  }
  N <- array(N,length(N))
  X <- array(X,length(X))
  Y <- array(Y,length(Y))
  invisible(cbind(X,Y,N))
}

# Vertices of polygon with n sides, width dx, height dy, center (x0,y0)
polygon.vertex <- function(x0,y0,dx,dy=NULL,n=6,rot=0) {
  if (is.null(dy)) dy <- dx
  if (!findInterval(rot,c(-2*pi,2*pi))) rot <- 0
  i <- seq(0,n-1,1)
  x <- x0 + dx*cos(2*pi*i/n + rot)
  y <- y0 + dy*sin(2*pi*i/n + rot)
  invisible(cbind(x,y))
}

# Plot polygon defined by polygon.vertex
polygon.fill <- function(x0,y0,dx,dy,n=6,col='white',border='black',rot=0) {
  xy <- polygon.vertex(x0,y0,dx,dy,n=n,rot=rot)
  polygon(xy[,1],xy[,2],col=col,border=border)
}

