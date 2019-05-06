vis.map <- function(x,col='red',map.type=NULL,
                    xrange=NULL,yrange=NULL,cex=1,
                    add.text=FALSE,cex.axis=NULL,
                    map.insert=TRUE,verbose=FALSE,
                    usegooglemap=TRUE,zoom=NULL,...) {
  if(verbose) {print('vis.map'); print(lon(x)); print(lat(x)); print(zoom)}
  ## KMP 2017-06-07 Weird problem: cex.axis is not found even though it is an argument to the function.
  ## It looks like cex.axis exists but when applying 'print' the following error message shows up: 
  ## 'Warning: restarting interrupted promise evaluation. Error in print(cex.axis) : object 'cex.axis' not found'
  cex.axis <- 0.7  # Temporary fix
  
  ## REB 2016-11-25: for dsensemble object
  if (is.null(lon(x))) attr(x,'longitude') <- lon(attr(x,'station'))
  if (is.null(lat(x))) attr(x,'latitude') <- lat(attr(x,'station'))
  if(is.null(xrange)) xrange <- range(lon(x)) + c(-5,5)
  if(is.null(yrange)) yrange <- range(lat(x)) + c(-2,2)
  if(!map.insert) new <- TRUE else new <- FALSE
  
  if (is.null(map.type)) {
    if( inherits(x,"field") | length(lon(x))!=length(lat(x)) |
        (length(lon(x))==2 & length(lat(x))==2) ) {
      map.type <- "rectangle"
    } else {
      map.type <- "points"
    }
  }
  
  ## REB: 2016-10-12 - add the possibility to use google maps
  ## KMP 2018-10-31: Don't use require inside the esd package. 
  ## Instead check if it the external package is installed and then 
  ## call it explicitly, e.g., RgoogleMaps::GetMap().
  ## Also add the package under 'Suggested' in the DESCRIPTION file.
  if (!requireNamespace("RgoogleMaps", quietly = TRUE)) {
    usegooglemap <- FALSE
  }
  
  if(usegooglemap) {
    if (is.null(zoom)) {
      if (verbose) print('zoom not defined')
      if (length(lon(x))==1) {
        zoom <- 5 
      } else {
        ## zoom = 12 is very local, zoom = 1 is the world
        mxdst <- max(diff(range(lat(x))),diff(range(lon(x))))
        zoom <- 1 - floor(0.75*log(mxdst/360))
      }
    }
    if (!is.finite(zoom)) zoom <- 5
    if (verbose) print(paste('zoom=',zoom))
    bgmap <- try(RgoogleMaps::GetMap(center=c(lat=mean(lat(x)),lon=mean(lon(x))),
                                     destfile = "map.station.esd.png",
                                     maptype = "mobile", zoom=zoom))
    if(inherits(bgmap,"try-error")) {
      usegooglemap <- FALSE
    } else {
      if(map.insert) {
        par(fig=c(0.75,0.95,0.75,0.95),new=TRUE,
            mar=c(0,0,0,0),xpd=NA,col.main="grey",bty="n")
      }
      if(map.type=="rectangle") {
        xx <- c(rep(max(lat(x)),2), rep(min(lat(x)),2), max(lat(x)))
        yy <- c(range(lon(x)), rev(range(lon(x))), min(lon(x)))
        RgoogleMaps::plotmap(xx, yy, bgmap, pch=19, col=col, cex=0.25)
        RgoogleMaps::PlotOnStaticMap(bgmap, lat=xx, lon=yy, lwd=1, col=col, FUN=lines, add=TRUE)
      } else {
        RgoogleMaps::plotmap(lat(x), lon(x), bgmap, pch=19, col=col, cex=2)
      }
    }
  }
  
  if(!usegooglemap) {
    if (verbose) {
      print('basic map')
      print(cex.axis)
    }
    data("geoborders", envir = environment())
    lon <- geoborders$x
    lat <- geoborders$y
    ok <- lon>(min(xrange)-1) & lon<(max(xrange)+1) &
      lat>(min(yrange)-1) & lat<(max(yrange)+1) &
      is.finite(lon) & is.finite(lat)
    lon2 <- attr(geoborders,"borders")$x
    lat2 <- attr(geoborders,"borders")$y
    ok2 <- lon2>(min(xrange)-1) & lon2<(max(xrange)+1) &
      lat2>(min(yrange)-1) & lat2<(max(yrange)+1) &
      is.finite(lon2) & is.finite(lat2)
    if (verbose) {print(sum(ok)); print(range(lon[ok])); print(range(lat[ok]))}
    if(map.insert) {
      par(fig=c(0.76,0.97,0.76,0.97),new=TRUE,
          mar=c(0,0,0,0),xpd=NA,col.main="grey",bty="n")
    } else {
      dev.new()
    }
    plot(lon[ok],lat[ok],lwd=1,col="black",type="p",pch='.',cex=2,
         #type='l', KMP 2016-03-16 problem with lines in map
         xlab=NA,ylab=NA,axes=FALSE,new=new,
         xlim=xrange,ylim=yrange)
    #xlim=range(c(lon[ok],lon2[ok2]),na.rm=TRUE),
    #ylim=range(c(lat[ok],lat2[ok2]),na.rm=TRUE))
    par(xpd=FALSE)
    lines(lon,lat) ## REB: 2016-11-25 need more solid lines.
    axis(1,mgp=c(3,0.5,0.3),cex.axis=cex.axis)
    axis(2,mgp=c(2,0.5,0.3),cex.axis=cex.axis)
    lines(lon2,lat2,col = "pink",lwd=1)
    #lines(lon2[ok2],lat2[ok2],col = "pink",lwd=1)
    if (verbose) print(map.type)
    if (map.type=="points") {
      if (verbose) {print(c(lon(x),lat(x),cex)); print(col)}
      points(lon(x),lat(x),pch=21,cex=cex,col=col,bg=col,lwd=1)
      if (add.text) text(lon(x),lat(x),labels=loc(x),col=col) 
    } else if (map.type=="rectangle") {
      rect(min(lon(x)),min(lat(x)),max(lon(x)),max(lat(x)),
           border="black",lwd=1,lty=2)
    }
  }
  if(verbose) print("exit vis.map")
}


vis.pca <- function(x,...,cex=1.5,new=TRUE,verbose=FALSE) {
  if(verbose) print("vis.pca")
  y <- x # quick fix
  col <- colscal(col=varid(y)); nc <- length(col)
  #if (is.precip(y)) col <- rev(col)
  lon <- attr(y,'longitude') 
  lat <- attr(y,'latitude') 
  N <- length(lon)
  R2 <- round(100*attr(y,'eigenvalues')^2/attr(y,'tot.var'),2)
  
  #print(N); print(length(attr(y,'mean')))
  m <- min(3,dim(attr(y,'pattern'))[2])
  
  # Set scale for colour scheme
  #str(y)
  a.T <- matrix(rep(NA,4*N),4,N)
  ax <- quantile(abs(attr(y,'mean')),0.99,na.rm=TRUE)
  if (min(attr(y,'mean'))<0) scale0 <- seq(-ax,ax,length=nc) else
    scale0 <- seq(0,ax,length=nc)
  ax <- quantile(abs(attr(y,'pattern')),0.99,na.rm=TRUE)
  scale <- seq(-ax,ax,length=nc)
  
  #print("here")
  for (i in 1:N) {
    a.T[1,i] <-  sum(attr(y,'mean')[i] > scale0)
    for (j in 1:m) 
      a.T[j+1,i] <-  sum(attr(y,'pattern')[i,j] > scale)
  }
  a.T[a.T < 1] <- 1; a.T[a.T > 100] <- 100
  
  if (new) dev.new(width=5,height=7)
  par(mfrow=c(3,2),mar=c(3.5,3,3.5,3),bty="n",xaxt="n",yaxt="n")
  
  plot(lon,lat,
       main="Climatology",
       col=col[a.T[1,]],pch=19,xlab="",ylab="",cex=cex)
  points(lon,lat,cex=cex)
  data("geoborders",envir=environment())
  lines(geoborders,col='grey40')
  lines(geoborders$x - 360,geoborders$y,col='grey40')
  points(lon,lat,cex=cex,col=col[a.T[1,]],pch=19)
  
  plot(lon,lat,
       main=paste("EOF #1:",R2[1],"% of variance"),
       col=col[a.T[2,]],pch=19,xlab="",ylab="",cex=cex)
  points(lon,lat,cex=cex)
  lines(geoborders)
  lines(geoborders$x - 360,geoborders$y)
  points(lon,lat,cex=cex,col=col[a.T[2,]],pch=19)
  
  plot(lon,lat,
       main=paste("EOF #2:",R2[2],"% of variance"),
       col=col[a.T[3,]],pch=19,xlab="",ylab="",cex=cex)
  points(lon,lat,cex=cex)
  lines(geoborders,col='grey40')
  lines(geoborders$x - 360,geoborders$y,col='grey40')
  points(lon,lat,cex=cex,col=col[a.T[3,]],pch=19)
  
  plot(lon,lat,
       main=paste("EOF #3:",R2[3],"% of variance"),
       col=col[a.T[4,]],pch=19,xlab="",ylab="",cex=cex)
  points(lon,lat,cex=cex)
  lines(geoborders,col='grey40')
  lines(geoborders$x - 360,geoborders$y,col='grey40')
  points(lon,lat,cex=cex,col=col[a.T[4,]],pch=19)
  
  par(mar=c(1,0,0,0),fig=c(0.1,0.3,0.665,0.695),new=TRUE,cex.axis=0.6)
  image(cbind(1:nc,1:nc),col=col)
  nl <- pretty(scale0)
  par(xaxt="s")
  axis(1,at=seq(0,1,length=length(nl)),labels=nl)
  
  par(mar=c(1,0,0,0),fig=c(0.1,0.3,0.32,0.35),new=TRUE,cex.axis=0.6,xaxt="n")
  image(cbind(1:nc,1:nc),col=col)
  nl <- pretty(scale)
  par(xaxt="s")
  axis(1,at=seq(0,1,length=length(nl)),labels=nl)
  
  par(mar=c(1,0,0,0),fig=c(0.6,0.8,0.665,0.695),new=TRUE,cex.axis=0.6,xaxt="n")
  image(cbind(1:nc,1:nc),col=col)
  nl <- pretty(scale)
  par(xaxt="s")
  axis(1,at=seq(0,1,length=length(nl)),labels=nl)
  
  par(mar=c(1,0,0,0),fig=c(0.6,0.8,0.32,0.35),new=TRUE,cex.axis=0.6,xaxt="n")
  image(cbind(1:nc,1:nc),col=col)
  nl <- pretty(scale)
  par(xaxt="s")
  axis(1,at=seq(0,1,length=length(nl)),labels=nl)
  
  par(mfcol=c(1,1),fig=c(0,1,0,0.33),new=TRUE,xaxt="s",yaxt="n",bty="n",
      mar=c(2,2,1,1))
  ylim <- 2*range(coredata(y[,1:m]),na.rm=TRUE)
  plot(y[,1]+0.5*ylim[2],lwd=2,ylim=ylim)
  grid()
  col <- c("red","blue")
  for (j in 1:m) lines(y[,j+1]+(1-j)*0.5*ylim[2],lwd=2,col=col[j])
  legend(index(y)[1],ylim[1],c('PC 1','PC 2','PC 3'),
         col=c('black','red','blue'),bty='n',lwd=2)
  invisible(a.T)
}