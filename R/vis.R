vis <- function(x,...) UseMethod("vis")

vis.default <- function(x,...,it=NULL,img=NULL,verbose=FALSE,
                        ref=c(as.Date('1961-01-01'),as.Date('1990-12-31'))) {
  if(verbose) print("vis.default")
  if (!is.null(img)) {
    if (is.character(img)) {
      ## KMP 2018-11-08: Don't use require in the R-scripts.
      ## Call the function explicitly and add the package under 'Suggests' in the DESCRIPTION file.
      ## Also check if the package is installed and add an alternative or error message if it isn't.
      #require(jpeg)
      if (requireNamespace("jpeg", quietly = TRUE)) {
        img <- jpeg::readJPEG(img)
      } else{
        stop("Package \"jpeg\" needed to read a jpeg file. Please install it.")
      }
    }
  }
  if (!is.null(it)) x <- subset(x,it=it)
  y <- attr(x,'station')
  dev.new()
  if (!is.null(img)) {
    par0 <- par()
    par(mar=rep(0,4))
    plot(c(0,1),c(0,1),type='n')
    rasterImage(img, -0.05, -0.05, 1.05, 1.05)
    par(new=TRUE,col.axis='white',col.lab='white',xaxt='s',yaxt='s',mar=par0$mar)  
  }
  par(bty='n')
  plot(zoo(y),lwd=5,col='black',ylim=range(y,na.rm=TRUE)+ c(-1,5),xlim=range(index(x)),
       ylab=ylab(y),xlab='Time')
  for (i in 1:dim(x)[2]) lines(x[,i],lwd=7,col=rgb(1,0.7,0.7,0.1))
  
  lines(y,lty=2,lwd=3)
  balls(y) 
}

vis.station <- function(x,...,new=FALSE,col=NULL,n=NULL,main=NULL,log.precip=TRUE,
                        plot=TRUE,verbose=FALSE) {
  if(verbose) print("vis.station")
  yrs <- as.numeric(rownames(table(year(x))))
  ny <- length(yrs)
  if (is.null(n)) n <- 50
  if (is.null(col)) col <- colscal(n,varid(x))
  if (is.precip(x)) x[coredata(x)==0] <- NA
  
  if ( (attr(x,'unit')[1] == "deg C") | (attr(x,'unit')[1] == "degree Celsius") )
    unit <- expression(degree*C) else
      unit <- attr(x,'unit')
  if (is.null(main)) eval(parse(text=paste("main <- expression(paste('Annual+seasonal evaluation of daily ',",
                                           attr(x,'variable'),"))")))
  
  z <- matrix(rep(NA*366*ny),366,ny)
  for (i in 1:ny) {
    y <- window(x,start=as.Date(paste(yrs[i],'-01-01',sep='')),
                end=as.Date(paste(yrs[i],'-12-31',sep='')))
    t <- julian(index(y)) - julian(as.Date(paste(yrs[i],'-01-01',sep=''))) +1
    it <- is.element(1:366,as.numeric(t))
    if (is.precip(x) & log.precip) z[it,i] <- log(y) else
      z[it,i] <- y
  }
  if (plot) { 
    if (new) dev.new()
    par(bty="n")
    image(1:366,yrs,z,main=main,xlab='',ylab='year',col=col,sub=loc(x),...)
    grid()
  }
  attr(z,'x') <- 1:366
  attr(z,'y') <- yrs
  invisible(z)
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

vis.map <- function(x,...,col='red',map.type=NULL,
                    xrange=NULL,yrange=NULL,cex=1,
                    add.text=FALSE,cex.axis=NULL,
                    map.insert=TRUE,verbose=FALSE,
                    usegooglemap=TRUE,zoom=NULL) {
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

vis.dsensemble <- function(x,...) {
  stopifnot(inherits(x,"dsensemble"))
  if (inherits(x,"list")) {
    vis.dsensemble.list(x,...)
  } else vis.default(x,...)
}

vis.dsensemble.list <- function(X,...,verbose=FALSE,FUN='trend',
                                colbar=NULL,legend.shrink=1,n=11,plim=0.01) {
  
  if (verbose) print('vis.dsensemble.list')
  stopifnot(inherits(X,"dsensemble") & inherits(X,"list"))
  
  if (inherits(X,"pca")) {
    if(verbose) print("as.station.dsensemble.pca")
    Y <- as.station(X)
  } else {
    Y <- X
  }
  
  if (is.null(attr(X,"unit"))) attr(X,"unit") <- attr(Y[[1]],"unit")
  if (is.null(attr(X,"variable"))) attr(X,"variable") <- attr(Y[[1]],"variable")
  
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
  d <- dim(Y[[1]])
  Z <- sapply(Y,function(x) {
    y <- apply(x,2,FUN=FUN)
    invisible( rbind(q5(y),median(y),q95(y)) )} )
  z.q5 <- Z[1,]
  z <- Z[2,]
  z.q95 <- Z[3,]
  
  if (!is.null(FUN2)) {
    if(verbose) print(paste("calculate significance"))
    p <- sapply(Y,function(x) mean(apply(x,2,FUN=FUN2)) )   
  } else {
    p <- rep(0,length(z))
  }
  
  if (is.null(colbar)) {
    colbar=list(palette='t2m',rev=FALSE,n=n,breaks=NULL,
                type="p",cex=2,h=0.6,v=1,pos=0.1,show=TRUE)
  }
  ## REB: 2015-12-10: drop these after colbar.ini has been revised
  #  if (is.null(colbar$breaks)) {
  #      colbar$breaks <- round(c(-max(abs(z.q95)),max(abs(z.q95))),ndig(z.q95))
  #  }
  
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
  #  colbar$breaks <- signif(colbar$breaks,digits=2)
  
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
  #col <- mapply(function(a,b) adjustcolor(a,alpha.f=b),col,alpha)
  #col.q5 <- mapply(function(a,b) adjustcolor(a,alpha.f=b),col.q5,alpha)
  #col.q95 <- mapply(function(a,b) adjustcolor(a,alpha.f=b),col.q95,alpha)
  
  data("geoborders",envir=environment())
  dx <- max(diff(range(lons))/10,2)
  dy <- max(diff(range(lats))/10,2)
  xrange <- range(lons) + c(-dx,dx)
  yrange <- range(lats) + c(-0.5*dy,2*dy)
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
  ##browser()
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
       col=adjustcolor("white",alpha.f=0.8),border="grey50")
  ## inner and outer layers
  points(min(xrange)+dlon/3,max(yrange)-dlat/3,pch=21,
         col="grey70",bg="grey80",cex=7,lwd=7)
  points(min(xrange)+dlon/3,max(yrange)-dlat/3,pch=21,col="grey95",
         bg="grey95",cex=2,lwd=0.5)
  text(min(xrange)+dlon/3,max(yrange)-dlat/3,"Low",cex=0.7)#"q5",cex=0.7)
  text(min(xrange)+dlon/3,max(yrange)-dlat*0.2,"Median",cex=0.7)#"mean",cex=0.7)
  text(min(xrange)+dlon/3,max(yrange)-dlat*0.05,"High",cex=0.7)#"q95",cex=0.7)
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
  text(min(xrange)+2.3*dlon,max(yrange)-0.55*dlat,pos=4,cex=0.7,"satisfactory")#"ok")
  text(min(xrange)+2.3*dlon,max(yrange)-0.75*dlat,pos=4,cex=0.7,"not good")#"bad")
  
  if(verbose) print("finished!")
  invisible(d)
}

vis.trends <- function(x,...,unitlabel="unit",varlabel="",is=1,
                       pmax=0.01,minlen=15,lwd=NA,vmax=NA,new=TRUE,
                       show.significance=TRUE,verbose=FALSE) {
  if(verbose) print("vis.trends")
  T <- calculate.trends(x,minlen=minlen,is=is,verbose=verbose)
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
  cmin <- rgb(239,138,98,maxColorValue=255) # blue
  cmid <- rgb(247,247,247,maxColorValue=255) # white
  cmax <- rgb(103,169,207,maxColorValue=255) # red
  rgb.palette <- colorRampPalette(c(cmax,cmid,cmin),space="rgb")
  cstep <- rgb.palette(n=length(vstep)-1)
  # Plot trend as color
  if (new) dev.new()
  image(cols,rows,t(trends),breaks=vstep,col=cstep,
        xlab='start year',ylab='length of period (years)',
        main=paste(c(varlabel," trend (",unitlabel,"/decade)"),collapse=""))
  
  trends.plus <- t(trends)
  trends.plus[trends.plus<max(vstep)] <- NA
  image(cols,rows,trends.plus,col=cstep[length(cstep)],add=TRUE)
  trends.minus <- t(trends)
  trends.minus[trends.minus>min(vstep)] <- NA
  image(cols,rows,trends.minus,col=cstep[1],add=TRUE)
  
  # Mark significant trends with dark borders
  if(show.significance) {
    if(verbose) print(paste("mark significant trends (p<",pmax,")",sep=""))
    i <- which((is.finite(t(p)) & t(p)<pmax))
    x <- array(sapply(cols,function(x) rep(x,nrow(p))),length(p))[i]
    y <- rep(rows,nrow(p))[i]
    matlines(rbind(x-1/2,x+1/2),rbind(y-1/2,y-1/2),col='black',lwd=lwd,lty=1)
    matlines(rbind(x-1/2,x+1/2),rbind(y+1/2,y+1/2),col='black',lwd=lwd,lty=1)
    matlines(rbind(x-1/2,x-1/2),rbind(y-1/2,y+1/2),col='black',lwd=lwd,lty=1)
    matlines(rbind(x+1/2,x+1/2),rbind(y-1/2,y+1/2),col='black',lwd=lwd,lty=1)
  }
  colbar(cticks,cstep,fig=c(0.85,0.9,0.65,0.85))
}

calculate.trends <- function(x,minlen=15,is=1,verbose=FALSE){
  # Calculate trends of time series x
  if(verbose) print("calculate.trends - calculate trends for all subperiods")
  stopifnot(inherits(x,'zoo'))
  if(!is.null(dim(x))) {
    x <- subset(x,is=is)
    if(!is.null(dim(x))) {
      x <- apply(x,2,mean,na.rm=TRUE)
    }
  }
  if(!inherits(x,c("annual","season"))) {
    xm <- aggregate(x,by=as.yearmon(index(x)),FUN="mean")
    xy <- aggregate(xm,by=strftime(index(xm),"%Y"),FUN="mean")
    ny <- aggregate(xm,by=strftime(index(xm),"%Y"),FUN="nv")
    xy <- xy[ny==max(ny)] # exclude years with missing months
  } else xy <- x
  year <- as.numeric(index(xy))
  firstyear <- min(year):(max(year)-minlen+1)
  trendlen <- minlen:(diff(range(year))+1)
  #lastyear <- firstyear+minlen-1
  n <- length(firstyear)
  trends <- matrix(NA,n,n)
  rownames(trends) <- trendlen#firstyear
  colnames(trends) <- firstyear#lastyear
  p <- trends
  # speed up with apply?
  for (i in firstyear) {
    jvec <- i+trendlen[trendlen<=(max(year)-i+1)]-1#(i+minlen-1):(max(year)+1)
    for (j in jvec) {
      if(!is.na(xy[year==i]) & !is.na(xy[year==j])) {
        #if(verbose) print(paste(i,j))
        ij <- which(year %in% i:j & !is.na(xy))
        ij.model <- lm(xy[ij]~year[ij])
        #ij.kendall <- Kendall(x[ij],year[ij])
        iout <- firstyear==i#as.numeric(colnames(trends))==j
        jout <- trendlen==(j-i+1)#as.numeric(rownames(trends))==i
        trends[jout,iout] <- ij.model$coefficients[2]
        p[jout,iout] <- anova(ij.model)$Pr[1]#ij.kendall$sl[1]
      }
    }
  }  
  return(list("trends"=trends,"p"=p))
}