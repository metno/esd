# Visualise - different type of plotting... 'Infographics' type.

vis <- function(x,...) UseMethod("vis")

vis.station <- function(x,...) {
  if (is.precip(x)) vis.station.precip(x,...) else
  if (is.T(x)) vis.station.t2m(x,...)
}

vis.station.precip <- function(x,p=c(seq(0.1,0.95,0.05),0.97,0.98,0.99),threshold=1,...) {
  # From qqplotter:
  x[x < threshold] <- NA
  if (is.null(dim(x))) {
    qp <- quantile(x,prob=p,na.rm=TRUE)
    qmu <- -log(1-p)*mean(x,na.rm=TRUE)
  } else {
    qp <- apply(coredata(x),2,quantile,prob=p,na.rm=TRUE)
    qmu <- -log(1-p)%o%apply(coredata(x),2,wetmean,na.rm=TRUE)
    #fw <- round(100*apply(coredata(x),2,wetfreq))
  }
  plot(qp,qmu,pch=19,col=rgb(0,0,1,0.2))
  lines(range(qp,qmu),range(qp,qmu))
  grid()
}

vis.station.t2m <- function(x,p=c(0.01,0.02,0.03,0.04,seq(0.1,0.95,0.05),
                                0.97,0.98,0.99),...) {
  d <- dim(x); if (is.null(d)) d <- c(length(x),1)
  if (is.null(dim(x))) {
    qp <- quantile(x,prob=p,na.rm=TRUE)
    qmu <-  qnorm(p=p,mean=mean(coredata(x),na.rm=TRUE),sd=sd(coredata(x),na.rm=TRUE))
  } else {
    qp <- apply(coredata(x),2,quantile,prob=p,na.rm=TRUE)
    qmu <- qp + NA
    for (i in 1:length(p))
      qmu[i,] <- qnorm(p=p[i],mean=apply(coredata(x),2,mean,na.rm=TRUE),
                       sd=apply(coredata(x),2,sd,na.rm=TRUE))
  }
  plot(qp,qmu,pch=19,col=rgb(1,0,0,0.2))
  lines(range(qp,qmu),range(qp,qmu))
  grid() 
}

vis.field <- function(x,...) {
}

vis.eof <- function(x,...) {
}

vis.spell <- function(x,...) {
}

vis.cca <- function(x,...) {
}

vis.mvr <- function(x,...) {
}

vis.ds <- function(x,...) {
}

vis.trends <- function(x,unitlabel="unit",varlabel="",
 pmax=0.01,minlen=15,lwd=NA,vmax=NA,new=TRUE) {

  T <- calculate.trends(x,minlen=minlen)
  trends <- T$trends*10
  p <- T$p
  cols <- as.numeric(colnames(trends))
  rows <- as.numeric(rownames(trends))
  significant <- ifelse(p<pmax,trends,NA)
  
  ticks <- seq(1,length(cols),signif(length(cols)/10,1))
  if (is.na(lwd)) lwd <- max(3-0.05*length(cols),0.2)
   
  if (is.na(vmax) | vmax=="q995") vmax <- q995(abs(trends))
  if (vmax=="max") vmax <- max(abs(trends),na.rm=T)
  if (vmax<1) vmax <- signif(vmax,1)
  if (vmax>1) vmax <- signif(vmax,2)
  dv <- signif(vmax/8,1)
  v0 <- 0#signif(dv/2,1)
  vstep <- seq(v0,vmax,dv)
  vstep <- unique(c(-1*vstep,vstep))
  vstep <- vstep[order(vstep)]
  cticks <- vstep[2:length(vstep)]-dv/2

  #cstep <- colscal(n=length(vstep)-1,col="t2m")
  cmin <- rgb(239,138,98,max=255) # blue
  cmid <- rgb(247,247,247,max=255) # white
  cmax <- rgb(103,169,207,max=255) # red
  rgb.palette <- colorRampPalette(c(cmax,cmid,cmin),space="rgb")
  cstep <- rgb.palette(n=length(vstep)-1)
  
  # Plot trend as color
  if (new) dev.new()
  image(rows,cols,t(trends),breaks=vstep,col=cstep,
        xlab='end year',ylab="start year",
        main=paste(c(varlabel," trend (",unitlabel,"/decade)"),collapse=""))

  trends.plus <- t(trends)
  trends.plus[trends.plus<max(vstep)] <- NA
  image(rows,cols,trends.plus,col=cstep[length(cstep)],add=TRUE)
  trends.minus <- t(trends)
  trends.minus[trends.minus>min(vstep)] <- NA
  image(rows,cols,trends.minus,col=cstep[1],add=TRUE)

  # Mark significant trends with dark borders
  i <- which((is.finite(t(p)) & t(p)<pmax))
  x <- rep(rows,nrow(p))[i]
  y <- array(sapply(cols,function(x) rep(x,nrow(p))),length(p))[i]
  matlines(rbind(x-1/2,x+1/2),rbind(y-1/2,y-1/2),col='black',lwd=lwd,lty=1)
  matlines(rbind(x-1/2,x+1/2),rbind(y+1/2,y+1/2),col='black',lwd=lwd,lty=1)
  matlines(rbind(x-1/2,x-1/2),rbind(y-1/2,y+1/2),col='black',lwd=lwd,lty=1)
  matlines(rbind(x+1/2,x+1/2),rbind(y-1/2,y+1/2),col='black',lwd=lwd,lty=1)

  colbar(cticks,cstep,fig=c(0.2,0.25,0.65,0.85))
}
 
calculate.trends <- function(x,minlen=10){
  # Calculate trends of time series x
  stopifnot(inherits(x,'zoo'))
  xm <- aggregate(x,by=as.yearmon(index(x)),FUN="mean")
  xy <- aggregate(xm,by=strftime(index(xm),"%Y"),FUN="mean")
  ny <- aggregate(xm,by=strftime(index(xm),"%Y"),FUN="nv")
  xy <- xy[ny==max(ny)] # exclude years with missing months 
  year <- as.numeric(index(xy))
  firstyear <- min(year):(max(year)-minlen+1)
  lastyear <- firstyear+minlen-1
  n <- length(firstyear)
  trends <- matrix(NA,n,n)
  colnames(trends) <- lastyear
  rownames(trends) <- firstyear
  p <- trends
  # speed up with apply?
  for (i in firstyear) {
    jvec <- (i+minlen-1):(max(year)+1)
    for (j in jvec) {
      ij <- which(year %in% i:j)
      ij.model <- lm(xy[ij]~year[ij])
      #ij.kendall <- Kendall(x[ij],year[ij])
      iout <- as.numeric(colnames(trends))==j
      jout <- as.numeric(rownames(trends))==i
      trends[jout,iout] <- ij.model$coefficients[2]
      p[jout,iout] <- anova(ij.model)$Pr[1]#ij.kendall$sl[1]
    }
  }  
  return(list("trends"=trends,"p"=p))
}

# Binned scatterplot with sunflowers
scatter.sunflower <- function(x,y,petalsize=7,dx=NULL,dy=NULL,
                          xgrid=NULL,ygrid=NULL,xlim=NULL,ylim=NULL,
                          xlab=NULL,ylab=NULL,main=NULL,leg=TRUE,rotate=TRUE,
                          alpha=0.6,leg.loc=2,new=TRUE) {

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
           border=adjustcolor('khaki3',alpha=alpha),n=6)
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
                           scale.col=TRUE,scale.size=FALSE) {

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


diagram <- function(x,...) UseMethod("diagram")

diagram.dsensemble <- function(x,it=0,...) {
  stopifnot(inherits(x,'dsensemble'))
  #print("subset") 
  if (!inherits(attr(x,'station'),'annual')) z <- subset(x,it=it) else
                                             z <- x
  y <- attr(z,'station')
  #print("diagnose")
  pscl <- c(0.9,1.1)
  if (max(coredata(z),na.rm=TRUE) < 0) pscl <- rev(pscl)
  #print("...")
  plot(y,type="b",pch=19,
       xlim=range(year(z)),
       ylim=pscl*range(coredata(z),na.rm=TRUE))
  grid()
  usr <- par()$usr; mar <- par()$mar; fig <- par()$fig
  t <- year(z); n <- dim(z)[2]
  col <- rgb(seq(1,0,length=n)^2,sin(seq(0,pi,length=n))^2,seq(0,1,length=n)^2,0.2)
  for (i in 1:n) lines(t,z[,i],col=col[i],lwd=2)
  points(y,pch=19,lty=1)
}

diagram.ds <- function(x,...) {
}

# Show the temperatures against the day of the year. Use
# different colours for different year.
diagram.station <- function(x,it=NULL,new=TRUE,...) {
  yrs <- as.numeric(rownames(table(year(x))))
  d <- dim(x)
  #print(yrs)
  ny <- length(yrs)
  j <- 1:ny
  col <- rgb(j/ny,abs(sin(pi*j/ny)),(1-j/ny),0.2)
  class(x) <- "zoo"

  if ( (attr(x,'unit') == "deg C") | (attr(x,'unit') == "degree Celsius") )
      unit <- expression(degree*C) else
      unit <- attr(x,'unit')
  eval(parse(text=paste("main <- expression(paste('Seasonal evaution: ',",
               attr(x,'variable'),"))")))
  if (new) dev.new()
  par(bty="n")
  z <- coredata(x)
  plot(c(0,365),1.25*range(z,na.rm=TRUE),
       type="n",xlab="",
       main=main,
       sub=attr(x,'location'),ylab=unit)
  grid()
  for (i in 1:ny) {
    y <- window(x,start=as.Date(paste(yrs[i],'-01-01',sep='')),
                    end=as.Date(paste(yrs[i],'-12-31',sep='')))
    t <- julian(index(y)) - julian(as.Date(paste(yrs[i],'-01-01',sep='')))
    if (is.null(d)) points(t,coredata(y),lwd=2,col=col[i],pch=19,cex=0.5) else
                    points(rep(t,d[2]),coredata(y),lwd=2,col=col[i],pch=19,cex=0.5)
  }
  if (!is.null(it)) {
    y <- window(x,start=as.Date(paste(it,'-01-01',sep='')),
                    end=as.Date(paste(it,'-12-31',sep='')))
    t <- julian(index(y)) - julian(as.Date(paste(it,'-01-01',sep='')))
  }
  if (is.null(d)) points(t,coredata(y),col="black",cex=0.7) else
                  points(rep(t,d[2]),coredata(y),col="black",cex=0.7)

  par(new=TRUE,fig=c(0.70,0.85,0.70,0.85),mar=c(0,3,0,0),
      cex.axis=0.7,yaxt="s",xaxt="n",las=1)
  colbar <- rbind(1:ny,1:ny)
  image(1:2,yrs,colbar,col=col)
}

#seNorge
#nodata	10000
#nobits	16
#header1	*Temperatur
#header2	*Siste døgn (18-18 UTC)
#legend	*Grader Celsius
#From	To	R	G	B	Forklaring
#2931	10000	204	0	0	*Over 20 
#2881	2931	255	25	0	*15 - 20
#2831	2881	255	102	0	*10 - 15
#2781	2831	255	179	0	*5 - 10
#2761	2781	255	230	77	*3 - 5
#2741	2761	255	255	64	*1 - 3
#2731	2741	255	255	190	*0 - 1
#2721	2731	217	255	255	*÷1 - 0
#2701	2721	179	255	255	*÷3 - ÷1
#2681	2701	128	235	255	*÷5 - ÷3
#2631	2681	64	204	255	*÷10 - ÷5
#2581	2631	0 	153	255	*÷15 - ÷10
#2531	2581	0 	25	255	*÷20 - ÷15
#0	2531	0	0	153	*Under ÷20


#nodata	10000
#nobits	16
#header1	*Nedbør
#header2	*Siste døgn (06-06 UTC)
#legend	*mm
#From	To	R	G	B	Forklaring
#1500	10000	0	0	153	*Over 150 
#750	1500	0 	25	255	*75 - 150 
#500	750	0 	153	255	*50 - 75
#300	500	64	204	255	*30 - 50
#200	300	128	235	255	*20 - 30
#100	200	179	255	255	*10 - 20
#1	100	217	255	255	*Under 10
#0	1	229	229	229	*Ikke nedbør




# Show the cumulative sum of station value from January 1st. Use
# different colours for different year.
cumugram <- function(x,it=NULL,prog=FALSE,verbose=FALSE,...) {
  stopifnot(!missing(x),inherits(x,"station"))
  
  #print("cumugram")
  yrs <- as.numeric(rownames(table(year(x))))
  today <- Sys.Date(); yesterday <- seq(today, length.out=2, by=-1)[2]
  
  #print(yrs)
  ny <- length(yrs)
  j <- 1:ny
  col <- rgb(j/ny,abs(sin(pi*j/ny)),(1-j/ny),0.3)
  class(x) <- "zoo"

  if ( (attr(x,'unit') == "deg C") | (attr(x,'unit') == "degree Celsius") )
      unit <- expression(degree*C) else
      unit <- attr(x,'unit')
  eval(parse(text=paste("main <- expression(paste('Running cumulative mean of  ',",
               attr(x,'variable'),sep=' ',"))")))
  dev.new()
  par(bty="n")
  z <- coredata(x)
  ylim <- c(NA,NA)

  #print('Find the y-range')
  y.rest <- rep(NA,ny); y2n <- y.rest
  ylim <- max(coredata(x),na.rm=TRUE) # to avoid getting warnings with empty vectors.
  for (i in 1:ny) {
    y <- window(x,start=as.Date(paste(yrs[i],'-01-01',sep='')),
                    end=as.Date(paste(yrs[i],'-12-31',sep='')))
    y.rest[i] <- mean(coredata(window(x,start=as.Date(paste(yrs[i],format(Sys.Date(),'-%m-%d'),sep='')),
                                      end=as.Date(paste(yrs[i],'-12-31',sep='')))))
    y2n[i] <- mean(coredata(window(x,end=as.Date(paste(yrs[i],format(Sys.Date()-1,'-%m-%d'),sep='')),
                                     start=as.Date(paste(yrs[i],'-01-01',sep='')))))                                  
    t <- julian(index(y)) - julian(as.Date(paste(yrs[i],'-01-01',sep='')))
    z <- cumsum(coredata(y))/1:length(y)
    ok <- is.finite(z)
    #rint(c(i,yrs[i],range(z[ok],na.rm=TRUE),ylim))
    ylim[!is.finite(ylim)] <- NA
    ylim[1] <- min(c(ylim,z[ok]),na.rm=TRUE)
    ylim[2] <- max(c(ylim,z[ok]),na.rm=TRUE)
  }
  #print(ylim)
  names(y2n) <- yrs
  y2n <- round(sort(y2n,decreasing=TRUE),2)
  
  plot(c(0,365),ylim,
       type="n",xlab="",
       main=main,sub=attr(x,'location'),ylab=unit,...)
  grid()

  cm <- rep(NA,ny)
  
  #browser()
  for (i in 1:ny) {
    y <- window(x,start=as.Date(paste(yrs[i],'-01-01',sep='')),
                    end=as.Date(paste(yrs[i],'-12-31',sep='')))
    t <- julian(index(y)) - julian(as.Date(paste(yrs[i],'-01-01',sep='')))
    z <- cumsum(coredata(y))/1:length(y)

    mm <- format(yesterday, "%m")
    dd <- as.numeric(yesterday, "%d")
    
    cm[i] <- mean(coredata(window(x,
            start=as.Date(paste(yrs[i],'-01-01',sep='')),
            end=as.Date(paste(yrs[i],mm,dd,sep='-')))))
    lines(t,z,lwd=2,col=col[i])
    print(c(i,yrs[i],range(z[ok],na.rm=TRUE),ylim))
  }
  if (is.null(it)) {
    lines(t,z,lwd=5,col="black")
    lines(t,z,lwd=2,col=col[i])
  } else {
    y <- window(x,start=as.Date(paste(it,'-01-01',sep='')),
                    end=as.Date(paste(it,'-12-31',sep='')))
    t <- julian(index(y)) - julian(as.Date(paste(it,'-01-01',sep='')))
    z <- cumsum(coredata(y))/1:length(y)   
    lines(t,z,lwd=5,col="black")
    lines(t,z,lwd=2,col=col[i])
  }
  tn <- t[length(t)]; 
  tm <- julian(as.Date('1900-12-31')) - julian(as.Date('1900-01-01'))
  zn <- coredata(z[length(z)-1])
  n <- 365
  #browser()
  zp <- length(z)/n * zn + (1-length(z)/n) * quantile(y.rest,0.95,na.rm=TRUE)
  zm <- length(z)/n * zn + (1-length(z)/n) * quantile(y.rest,0.05,na.rm=TRUE)
  zz <- length(z)/n * zn + (1-length(z)/n) * mean(y.rest,na.rm=TRUE)
  if (prog) {
    polygon(c(tn,rep(tm,2),tn),c(zn,zp,zm,zn),
            col=rgb(0.5,0.5,0.5,0.1),border=rgb(0.5,0.5,0.5,0.2),lwd=2)
    lines(c(tn,tm),c(zn,zz),col=rgb(0.3,0.3,0.3,0.1),lwd=3)
    text(tm,zp,round(zp,1),pos=3,cex=0.5,col='grey40')
    text(tm,zm,round(zm,1),pos=1,cex=0.5,col='grey40')
    text(tm,zz,round(zz,1),pos=4,cex=0.75)
    print(paste('Prognosis for end-of-year: ',round(zz,1),' (',round(zm,1),',',round(zp,1),')',sep=''))
  }

  if (!is.precip(x))
    par(new=TRUE,fig=c(0.70,0.85,0.20,0.35),mar=c(0,3,0,0),
        cex.axis=0.7,yaxt="s",xaxt="n",las=1)
  else
    par(new=TRUE,fig=c(0.70,0.85,0.70,0.85),mar=c(0,3,0,0),
        cex.axis=0.7,yaxt="s",xaxt="n",las=1)
  colbar <- rbind(1:ny,1:ny)
  image(1:2,yrs,colbar,col=col)

  srt <- order(cm,decreasing=TRUE)
  invisible(cbind(yrs[srt],cm[srt]))
  if (verbose) print(y2n)
  invisible(y2n)
}

# Estimate how the variance varies with season 
# sd from inter-annual variability of daily values

climvar <- function(x,FUN='sd',plot=TRUE,...) {
  yrs <- as.numeric(rownames(table(year(x))))
  #print(yrs)
  ny <- length(yrs)
  X <- x; class(X) <- "zoo"
  
  if ( (attr(x,'unit') == "deg C") | (attr(x,'unit') == "degree Celsius") )
      unit <- expression(degree*C) else
      unit <- attr(x,'unit')
  eval(parse(text=paste("main <- expression(paste('seasonal ",
#               deparse(substitute(FUN))," of ',",
               FUN," of ',",attr(x,'variable'),"))")))
  Z <- matrix(rep(NA,ny*365),ny,365)
  
  for (i in 1:ny) {
    y <- window(X,start=as.Date(paste(yrs[i],'-01-01',sep='')),
                    end=as.Date(paste(yrs[i],'-12-31',sep='')))
    t <- julian(index(y)) - julian(as.Date(paste(yrs[i],'-01-01',sep='')))
    Z[i,] <- approx(t,y,1:365)$y
  }
  z <- apply(Z,2,FUN,na.rm=TRUE,...)
  wt <- 2*pi*(1:365)/365
  s1 <- sin(wt); c1 <- cos(wt); s2 <- sin(2*wt); c2 <- cos(2*wt)
  s3 <- sin(3*wt); c3 <- cos(3*wt); s4 <- sin(4*wt); c4 <- cos(4*wt)
  acfit <- predict(lm(z ~s1 + c1 + s2 + c2 + s3 + c3 + c4 + s4))
    
  if (plot) {
    dev.new()
    par(bty="n")
    plot(c(0,365),range(z,na.rm=TRUE),
         type="n",xlab="",
         main=main,
        sub=attr(x,'location'),ylab=unit)
    grid()
    lines(z,lwd=5)
    lines(z,lwd=3,col="grey")
    lines(acfit,lwd=5)
    lines(acfit,lwd=3,col="red")

    par(new=TRUE,fig=c(0.15,0.35,0.70,0.90),mar=c(0,0,0,0),
        yaxt="n",xaxt="n",las=1)
    plot(c(0,1),c(0,1),type="n",xlab="",ylab="")
    legend(0,1,c("raw data","harmonic fit"),lwd=3,col=c("grey","red"),bty="n",cex=0.6)  
  }
  
  acfit <- attrcp(x,acfit)
  attr(acfit,'raw_data') <- z
  attr(acfit,'history') <- history.stamp(x)
  return(z)
}



histwet <- function(x,breaks=NULL,threshold=1) {
  if (is.null(breaks)) breaks=seq(0,1.1*max(x,na.rm=TRUE),by=5)
  h <- hist(x[x > threshold],breaks=breaks,plot=FALSE)
  return(h$density)
}

## Some additional infographics: show p.d.f. for each year    
visprob <- function(x,...) UseMethod("visprob")

visprob.default <- function(x,...) {

}

visprob.station <- function(x,y=NULL,is=1,dy=0.01,verbose=FALSE,...) {
  if (is.precip(x)) visprob.station.precip(x,y=y,is=is,
                                           dy=dy,verbose=verbose,...) 
}

visprob.station.precip <- function(x,y=NULL,is=1,threshold=1,dy=0.005,
                                   breaks=NULL,pdf=FALSE,verbose=FALSE,...) {
## Plot the histogram for each year in different colours, depending on y. Iy y
## is NULL, use the year to set the colour
  if (verbose) print('visprob.station.precip')
  ## If y is provided, synchronise the two time series:
  if (!is.null(y)) {
    y <- subset(y,it=c(start(x),end(x)))
    x <- subset(x,it=c(start(y),end(y)))
    y <- annual(y)
  } else y <- year(annual(x))
  mu <- aggregate(x,year,FUN='wetmean',threshold=threshold)
  fw <- aggregate(x,year,FUN='wetfreq',threshold=threshold)
  col <- colscal(n=length(y),alpha=coredata(fw))
  srtc <- order(y)
  col <- col[srtc]
  if (is.null(breaks))
    breaks <- seq(floor(min(x,na.rm=TRUE)),ceiling(max(x,na.rm=TRUE))+5,by=5)
  z <- aggregate(x,year,FUN='histwet',breaks=breaks,threshold=threshold)
  if (verbose) print(c(dim(z),length(y)))
  dy <- abs(max(z,na.rm=TRUE)*dy)
  mids <- 0.5*(breaks[-1] + breaks[-length(breaks)])
  par(bty='n',yaxt='n')
  plot(range(breaks),c(0,max(z,na.rm=TRUE) + length(y)*dy),type='n',
  ylab='f(x)',xlab=paste(varid(x),unit(x)),
       main=paste('Statistical distribution for',loc(x)),...)
  for (i in seq(length(y),1,by=-1)) {
    lines(mids,z[i,]+dy*i,col="grey",lwd=5)
    lines(mids,z[i,]+dy*i,col=col[i],lwd=4)
  }
  if (pdf) {
    if (verbose) print('add pdfs')
    for (i in 1:length(y)) {
      lines(mids,dy*i + exp(-mids/coredata(mu[i]))/coredata(mu[i]),
            col='black',lty=2)
    }
  }
  if (!is.null(loc(y)))
    text(par()$xaxp[2],par()$yaxp[2],loc(y),pos=2)
}


vis.dsensemble <- function(x,...) {
  stopifnot(inherits(x,"dsensemble"))
  if (inherits(x,"list")) {
    vis.dsensemble.list(x,...)
  } 
}

vis.dsensemble.list <- function(X,verbose=FALSE,FUN='trend',
    colbar=NULL,legend.shrink=1,n=11,plim=0.01,...) {

  if (verbose) print('vis.dsensemble.list')
  stopifnot(inherits(X,"dsensemble") & inherits(X,"list"))

  if (inherits(X,"pca")) {
    if(verbose) print("as.station.dsensemble.pca")
    Y <- as.station(X)
  } else {
    Y <- X
  }

  if (is.null(attr(Y,"unit"))) attr(Y,"unit") <- attr(Y[[1]],"unit")
  if (is.null(attr(Y,"variable"))) attr(Y,"variable") <- attr(Y[[1]],"variable")
  
  gcms <- attr(Y[[1]],"model_id")
  lons <- sapply(Y,lon)
  lats <- sapply(Y,lat)

  if (FUN=='trend') {
    FUN <- trend.coef
    FUN2 <- trend.pval
    label_fun <- "trend"
    attr(X,"unit") <- paste(attr(X,"unit"),"/decade",sep="")
  } else {
    FUN <- trend.coef
    FUN2 <- NULL
    label_fun <- FUN
  }
  if(verbose) print(paste("calculate",label_fun))
  Z <- matrix(rep(NA,length(gcms)*length(Y)),length(Y),length(gcms))
  for (i in 1:length(Y)) {
     Z[i,] <- apply(Y[[i]],2,FUN)
  }
  z <- apply(Z,1,mean)
  z.q5 <- apply(Z,1,q5)
  z.q95 <- apply(Z,1,q95)
  
  if (!is.null(FUN2)) {
    if(verbose) print(paste("calculate significance"))
    P <- matrix(rep(NA,length(gcms)*length(Y)),length(Y),length(gcms))
    for (i in 1:length(Y)) {
      P[i,] <- apply(Y[[i]],2,FUN2)
    }
    p <- apply(P,1,mean)
  } else p <- rep(0,length(z))
   
  if (is.null(colbar)) {
      colbar=list(palette='t2m',rev=FALSE,n=n,breaks=NULL,
          type="p",cex=2,h=0.6,v=1,pos=0.1,show=TRUE)
  }
  if (is.null(colbar$breaks)) {
      colbar$breaks <- c(-max(abs(z.q95)),max(abs(z.q95)))
  }

  if (inherits(X,"pca")) {
    xval <- lapply(X[3:length(X)],function(x) attr(x,"evaluation"))
    r2 <- sapply(xval,function(x) 100*cor(x[,1],x[,2])^2)
    if(verbose) {
      print(paste("calibration period",
        paste(range(year(index(xval[[1]]))),collapse="-")))
      print(paste("x-validation: r2 = ",round(mean(r2)),"%",
        " (",round(q5(r2)),"%,",round(q95(r2)),"%)",sep=""))
    }
  }
  
  if(verbose) print("diagnose")
  d <- diagnose(Y,verbose=verbose,plot=FALSE)
  delta <- abs(2*(0.5-pnorm(d$deltaobs,mean=mean(d$deltagcm),
                             sd=sd(d$deltagcm))))
  outside <- abs(2*(0.5-pbinom(d$outside,size=d$N,prob=0.1)))
  quality <- 1-0.5*(delta + outside)
  cex <- 1 + 2*quality
  
  pch <- rep(21,length(p))
  pch[p>0.05] <- 23
  cex[pch==23] <- cex[pch==23]*0.9
  
  if(verbose) print('shape - significance of trend')
  if(verbose) print('size - quality of fit (magnitude & trend)')

  colbar <- colbar.ini(z,colbar=colbar,verbose=verbose)
  colbar$breaks <- signif(colbar$breaks,digits=2)

  icol <- apply(as.matrix(z),2,findInterval,colbar$breaks)
  col <- colbar$col[icol]
  icol.q5 <- apply(as.matrix(z.q5),2,findInterval,colbar$breaks)
  col.q5 <- colbar$col[icol.q5]
  icol.q95 <- apply(as.matrix(z.q95),2,findInterval,colbar$breaks)
  col.q95 <- colbar$col[icol.q95]
  col[as.matrix(z)>max(colbar$breaks)] <- colbar$col[n-1]
  col[as.matrix(z)<min(colbar$breaks)] <- colbar$col[1]
  col.q5[as.matrix(z.q5)>max(colbar$breaks)] <- colbar$col[n-1]
  col.q5[as.matrix(z.q5)<min(colbar$breaks)] <- colbar$col[1]
  col.q95[as.matrix(z.q95)>max(colbar$breaks)] <- colbar$col[n-1]
  col.q95[as.matrix(z.q95)<min(colbar$breaks)] <- colbar$col[1]
  
  #alpha <- p
  #alpha[p>0.05] <- 0.5
  #col <- mapply(function(a,b) adjustcolor(a,alpha=b),col,alpha)
  #col.q5 <- mapply(function(a,b) adjustcolor(a,alpha=b),col.q5,alpha)
  #col.q95 <- mapply(function(a,b) adjustcolor(a,alpha=b),col.q95,alpha)

  data("geoborders",envir=environment())
  dx <- max(diff(range(lons))/10,2)
  dy <- max(diff(range(lats))/10,2)
  xrange <- range(lons) + c(-dx,dx)
  yrange <- range(lats) + c(-0.5*dy,1.5*dy)
  dlon <- diff(xrange)/6
  dlat <- diff(yrange)/5

  if(verbose) print("plot map of trends")
  dev.new()
  par0 <- par()
  par(mgp=c(0,0.5,0),lheight=.9)
  plot(xrange,yrange,type="n",xlab=NA,ylab=NA,axes=FALSE)
  title(paste("dsensemble ",label_fun,", ",paste(year(range(index(Y[[1]]))),
        collapse="-"),sep=""),line=3)
  axis(2)
  axis(3)
  lines(geoborders$x,geoborders$y,col="darkblue")
  lines(attr(geoborders,'borders')$x,attr(geoborders,'borders')$y,col="pink")
  lines(geoborders$x+360,geoborders$y,col="darkblue")
  points(lons,lats,pch=pch,col=col.q95,bg=col,cex=cex,lwd=cex)
  points(lons,lats,pch=pch,col=col,bg=col.q5,cex=cex/3,lwd=0.5)
  text(mean(xrange),min(yrange),cex=1,labels=paste(attr(X,"var"),
    " ",label_fun," (",attr(X,"unit"),")",sep="")) 
  if (colbar$show) {
    if(verbose) print("add colorbar")
    par(fig=par0$fig,new=TRUE)
    image.plot(lab.breaks=colbar$breaks,horizontal = TRUE,
               legend.only = T, zlim = range(colbar$breaks),
               col = colbar$col, legend.width = 1,
               border = FALSE,legend.shrink=legend.shrink,
               axis.args = list(cex.axis = 0.8,
               xaxp=c(range(colbar$breaks),n=colbar$n)))
  }

  ## Add legends
  if(verbose) print("add legends")
  rect(min(xrange)-0.1*dlon,max(yrange)-0.9*dlat,
       mean(xrange)+0.3*dlon,max(yrange)+0.2*dlat,
       col=adjustcolor("white",alpha=0.8),border="grey50")
  ## inner and outer layers
  points(min(xrange)+dlon/3,max(yrange)-dlat/3,pch=21,
         col="grey70",bg="grey80",cex=7,lwd=7)
  points(min(xrange)+dlon/3,max(yrange)-dlat/3,pch=21,col="grey95",
         bg="grey95",cex=2,lwd=0.5)
  text(min(xrange)+dlon/3,max(yrange)-dlat/3,"q5",cex=0.7)
  text(min(xrange)+dlon/3,max(yrange)-dlat*0.2,"mean",cex=0.7)
  text(min(xrange)+dlon/3,max(yrange)-dlat*0.05,"q95",cex=0.7)
  ## size
  points(min(xrange)+2.0*dlon,max(yrange)-0.3*dlat,pch=21,
         bg="grey70",col="grey70",cex=3)
  points(min(xrange)+1.3*dlon,max(yrange)-0.3*dlat,
         pch=23,bg="grey70",col="grey70",cex=3*0.9)
  points(min(xrange)+2.0*dlon,max(yrange)-0.55*dlat,
         pch=21,bg="grey70",col="grey70",cex=2)
  points(min(xrange)+1.3*dlon,max(yrange)-0.55*dlat,
         pch=23,bg="grey70",col="grey70",cex=2*0.9)
  points(min(xrange)+2.0*dlon,max(yrange)-0.75*dlat,
         pch=21,bg="grey70",col="grey70",cex=1)
  points(min(xrange)+1.3*dlon,max(yrange)-0.75*dlat,
         pch=23,bg="grey70",col="grey70",cex=0.9)
  text(min(xrange)+2.0*dlon,max(yrange),"sign. trend\n(95%-level)",cex=0.7)
  text(min(xrange)+1.3*dlon,max(yrange),"not sign.\n trend",cex=0.7)
  text(min(xrange)+2.3*dlon,max(yrange)-0.3*dlat,
       "quality of esd:\nexcellent",pos=4,cex=0.7)
  text(min(xrange)+2.3*dlon,max(yrange)-0.55*dlat,"ok",pos=4,cex=0.7)
  text(min(xrange)+2.3*dlon,max(yrange)-0.75*dlat,"bad",pos=4,cex=0.7)

  if(verbose) print("finished!")
  invisible(d)
}
