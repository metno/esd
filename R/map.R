# R.E. Benestad
# Plot a map of the station locations, fields, EOFs, CCA results, correlation, composites, ...

#require(zoo)

map <- function(x,it=NULL,is=NULL,new=TRUE,...) UseMethod("map")

map.default <- function(x,it=NULL,is=NULL,new=TRUE,projection="lonlat",
                        xlim=NULL,ylim=NULL,zlim=NULL,n=15,
                        col=NULL,breaks=NULL,
                        what=NULL,gridlines=FALSE,
                        lonR=NULL,latR=-90,axiR=NULL,verbose=FALSE,...) {

  # default with no arguments will produce a map showing the station data in the esd package.

  if (verbose) print('map.default')
  x <- subset(x,it=it,is=is)
  X <- attr(x,'pattern')

  ## if zlim is specified, then mask data outside this range
  if (!is.null(zlim)) {
    d <- dim(X)
    mask <- (X < min(zlim)) | (X > max(zlim))
    X[mask] <- NA
    dim(X) <- d
    if (verbose) {print(zlim); print(dim(X)); print(sum(mask))}
  }
  attr(X,'longitude') <- lon(x)
  attr(X,'latitude') <- lat(x)
  attr(X,'variable') <- attr(x,'variable')
  attr(X,'unit') <- attr(x,'unit')[1]
  attr(X,'source') <- attr(x,'source')
  attr(X,'variable') <- varid(x)
  if (inherits(X,'zoo')) attr(X,'time') <- range(index(x)) else
  if (!is.null(attr(x,'time'))) attr(X,'time') <- attr(x,'time')
  if (projection=="lonlat") lonlatprojection(x=X,xlim=xlim,ylim=ylim,n=n,
                                             col=col,breaks=breaks,
                                             what=what,new=new,
                                             gridlines=gridlines,...) else
  if (projection=="sphere") map2sphere(x=X,lonR=lonR,latR=latR,axiR=axiR,
                                       what=what,gridlines=gridlines,
                                       col=col,new=new,...) else
  if (projection=="np") map2sphere(X,lonR=lonR,latR=latR,axiR=axiR,
                                   what=what,gridlines=gridlines,
                                   col=col,new=new,...) else
  if (projection=="sp") map2sphere(X,lonR=lonR,latR=latR,axiR=axiR,new=new,
                                   what=what,gridlines=gridlines,col=col,...)
  
  invisible(X)
}
 
map.matrix <- function(x,new=TRUE,projection="lonlat",...) {
  
# If x is provided, map only x...

  # default with no arguments will produce a map showing the station data in the esd package.

#  image(lon(x),lat(x),x)
  if (inherits(x,'zoo')) attr(x,'time') <- range(index(x)) 
  if (projection=="lonlat") lonlatprojection(x=x,new=new,...)  else
  if (projection=="sphere") map2sphere(x=x,new=new,...) else
  if (projection=="np") map2sphere(x,new=new,...) else
  if (projection=="sp") map2sphere(x,new=new,...)
  
  #map.station(NULL,...)
}

map.array <- function(x,pattern=1,new=TRUE,verbose=FALSE,
                      projection="lonlat",...) {
  if (verbose) print('map.array')
  z <- x[,,pattern]
  attr(z,'longitude') <- lon(x)
  attr(z,'latitude') <- lat(x)
  attr(z,'variable') <- varid(x)
  attr(z,'unit') <- unit(x)[1]

  map(z,new=new,projection=projection,verbose=verbose,...)
}


map.comb <- function(x,it=NULL,is=NULL,new=TRUE,
                     xlim=NULL,ylim=NULL,zlim=NULL,
                     pattern=1,n=15,
                     projection="lonlat",col=NULL,breaks=NULL,
                     lonR=NULL,latR=NULL,axiR=0,what=c("fill","contour"),
                     gridlines=TRUE,verbose=FALSE,...) {
  if (verbose) print('map.comb')
  stopifnot(inherits(x,'eof'))
  x <- subset(x,it=it,is=is)
  projection <- tolower(projection)
  if (is.null(col)) col <- colscal(n=n-1) else
  if (length(col)==1) {
     palette <- col
     col <- colscal(col=palette,n=n-1)
  }
  if (is.null(varid(x))) attr(x,'variable') <- 'NA'
  if (tolower(varid(x))=='precip') col <- rev(col) 
  
  map.eof(x=x,xlim=xlim,ylim=ylim,zlim=zlim,pattern=pattern,
          n=n,projection=projection,col=col,new=new,
          breaks=breaks,lonR=lonR,latR=latR,axiR=axiR,what=what,
          gridlines=gridlines,verbose=verbose,...) -> result
  invisible(result)
 }

map.eof <- function(x,it=NULL,is=NULL,new=TRUE,pattern=1,
                    xlim=NULL,ylim=NULL,zlim=NULL,n=15,
                    projection="lonlat",col=NULL,
                    breaks=NULL,lonR=NULL,latR=NULL,axiR=0,
                    what=c("fill","contour"),gridlines=TRUE,verbose=FALSE,...) {
  if (verbose) print('map.eof')
  stopifnot(inherits(x,'eof'))
  #x <- subset(x,it=it,is=is)
  projection <- tolower(projection)
  tot.var <- attr(x,'tot.var')
  D <- attr(x,'eigenvalues')
  var.eof <- 100* D^2/tot.var
  X <- attr(x,'pattern')[,,pattern]

  ## if zlim is specified, then mask data outside this range
  if (!is.null(zlim)) {
    d <- dim(X)
    mask <- (X < min(zlim)) | (X > max(zlim))
    X[mask] <- NA
    dim(X) <- d
    if (verbose) {print(zlim); print(dim(X)); print(sum(mask))}
  }
  #str(x)
  attr(X,'longitude') <- attr(x,'longitude')
  attr(X,'latitude') <- attr(x,'latitude')
  attr(X,'variable') <- attr(x,'variable')
  attr(X,'unit') <- attr(x,'unit')[1]
  attr(X,'source') <- attr(x,'source')
  attr(X,'time') <- range(index(x))
  if ( (pattern==1) & !is.null(attr(x, "area.mean.expl")) )
    if (attr(x, "area.mean.expl")) what="fill"
  if (projection=="lonlat") lonlatprojection(x=X,xlim=xlim,ylim=ylim,
                             n=n,col=col,breaks=breaks,new=new,
                             what=what,gridlines=gridlines,
        verbose=verbose,...) else
  if (projection=="sphere") map2sphere(x=X,lonR=lonR,latR=latR,axiR=axiR,
                                       what=what,gridlines=gridlines,
                                       col=col,new=new,verbose=verbose,...) else
  if (projection=="np") map2sphere(X,lonR=lonR,latR=90,axiR=axiR,
                                       what=what,gridlines=gridlines,
                                       col=col,new=new,verbose=verbose,...) else
  if (projection=="sp") map2sphere(X,lonR=lonR,latR=-90,axiR=axiR,
                                       what=what,gridlines=gridlines,
                                       col=col,new=new,verbose=verbose,...)
  invisible(X)
}


map.ds <- function(x,it=NULL,is=NULL,new=TRUE,xlim=NULL,ylim=NULL,
                   what=c("fill","contour"),pattern=1,
                   n=15,projection="lonlat",
                   lonR=NULL,latR=NULL,axiR=0,gridlines=TRUE,
                   col=NULL,breaks=NULL,verbose=FALSE,...) {
  if (verbose) print('map.ds')
  stopifnot(inherits(x,'ds'))
  x <- subset(x,it=it,is=is)

## REB 2015-03-26
  if (inherits(x,'pca')) {
    map.pca(x,pattern=pattern,verbose=verbose,new=new,
            xlim=xlim,ylim=ylim,projection=projection,
            lonR=lonR,latR=latR,axiR=axiR,gridlines=gridlines,
            col=col,breaks=breaks)
    return()
  }
  projection <- tolower(projection)
  X <- attr(x,'pattern')
  if (is.list(X)) {
    X <- X[[1]]
  }
  
  # Check if there are several patterns: one for each month/seasons
  d <- dim(X)
  if (length(d)>2) {
    dim(X) <- c(d[1],d[2]*d[3])
    X <- colMeans(X)
    dim(X) <- c(d[2],d[3])
    attr(X,'longitude') <- lon(attr(x,'pattern'))
    attr(X,'latitude') <- lat(attr(x,'pattern'))
  }
  attr(X,'variable') <- varid(x)
  attr(X,'unit') <- unit(x)[1]
  
  unit <- attr(x,'unit')
  if ( (is.na(unit) | is.null(unit)) ) unit <- " "
  for (i in 1:length(unit)) {
    if ((unit[i]=='degree Celsius') | (unit[i]=='deg C') | (unit[i]=='degC'))
         unit[i] <- 'degree*C'
  }
  
  attr(X,'unit') <- unit
  attr(X,'source') <- attr(x,'source')
  #dim(X) <- attr(x,'dimensions')[1:2]
  #print(c(dim(X),length(attr(X,'longitude')),length(attr(X,'longitude'))))

  if (projection=="lonlat") {
    lonlatprojection(x=X,n=n,col=col,breaks=breaks,
                     what='fill',colorbar=FALSE,
                     gridlines=gridlines,new=new,...)
    if (is.list(attr(x,'pattern'))) {
      Xa <- attr(x,'pattern')
      nms <- names(Xa)
      col <- c('black','darkgreen','grey','yellow','magenta','cyan',
               'brown','white','green')
      #browser()
      for (i in (2:length(nms))) 
        contour(lon(Xa[[i]]),lat(Xa[[i]]),Xa[[i]],add=TRUE,col=col[i])
    } else if (sum(is.element(what,'contour'))>0)
       contour(lon(X),lat(X),X,add=TRUE,col=col[1])
  } else
  if (projection=="sphere") map2sphere(x=X,lonR=lonR,latR=latR,axiR=axiR,
                                       what=what,gridlines=gridlines,
                                       col=col,new=new,...) else
  if (projection=="np") map2sphere(X,lonR=lonR,latR=90,axiR=axiR,
                                       what=what,gridlines=gridlines,
                                       col=col,new=new,...) else
  if (projection=="sp") map2sphere(X,lonR=lonR,latR=-90,axiR=axiR,
                                   what=what,gridlines=gridlines,
                                   col=col,new=new,...)
  invisible(X)
}


map.field <- function(x,it=NULL,is=NULL,new=TRUE,
                      xlim=NULL,ylim=NULL,zlim=NULL,
                      what=c("fill","contour"),
                      FUN='mean',n=15,projection="lonlat",
                      lonR=NULL,latR=NULL,na.rm=TRUE,colorbar=TRUE,
                      axiR=0,gridlines=FALSE,col=NULL,breaks=NULL,
                      verbose=FALSE,...) {

  stopifnot(inherits(x,'field'))
  if (verbose) print('map.field')
  x <- subset(x,it=it,is=is)
  #print(length(x)); print(attr(x,'dimensions')[1:2])
  projection <- tolower(projection)
  if (FUN=='trend') FUN <- 'trend.coef'

  if (!is.null(xlim)) {
    if (xlim[1] < 0) x <- g2dl(x,greenwich=FALSE)
  }
  #str(X)
  X <- coredata(x)
  

  ## If one time slice, then map this time slice
  if (dim(X)[1]==1) X <- coredata(x[1,]) else
  if (is.null(X)) X <- coredata(X) else if (inherits(X,"matrix")) {
  ## If several time slices, map the required statistics
    good <- apply(coredata(x),2,nv) > 1
    X <- rep(NA,length(good))
    xx <- coredata(x[,good])
    X[good] <- apply(xx,2,FUN=FUN,na.rm=na.rm)
  }

    ## if zlim is specified, then mask data outside this range
  if (!is.null(zlim)) {
    d <- dim(X)
    mask <- (X < min(zlim)) | (X > max(zlim))
    rng <- range(X,na.rm=TRUE)
    X[mask] <- NA
    dim(X) <- d
    if (verbose) print(paste('zlim=',zlim[1],zlim[2],
                             '  sum(mask)=',print(sum(mask)),
                             '  range(X)=',rng[1],rng[2]))
  }
  
  #print(length(X))
  attr(X,'longitude') <- attr(x,'longitude')
  attr(X,'latitude') <- attr(x,'latitude')
  attr(X,'variable') <- attr(x,'variable')[1]
#  if (attr(x,'unit')=="deg C") attr(X,'unit') <- expression(degree*C) else
  unit <- attr(x,'unit')[1]
  if ( (is.na(unit) | is.null(unit)) ) unit <- " "
  if ((unit=='degree Celsius') | (unit=='deg C') | (unit=='degC'))
       unit <- 'degree*C'

  attr(X,'unit') <- unit
  attr(X,'source') <- attr(x,'source')
  attr(X,'time') <- range(index(x))
  attr(X,'method') <- FUN
  attr(X,'timescale') <- class(x)[2]
  #print(length(X)); print(attr(x,'dimensions'))
  dim(X) <- attr(x,'dimensions')[1:2]
  #class(X) <- class(x)
  #str(X)
  
  if (projection=="lonlat") lonlatprojection(x=X,xlim=xlim,ylim=ylim,n=n,
                                             col=col,breaks=breaks,
                                             what=what,new=new,
                                             gridlines=gridlines,...) else
  if (projection=="sphere") map2sphere(x=X,lonR=lonR,latR=latR,axiR=axiR,
                                       what=what,gridlines=gridlines,
                                       col=col,new=new,
                                       breaks=breaks,colorbar=colorbar,...) else
  if (projection=="np") map2sphere(X,lonR=lonR,latR=90,axiR=axiR,
                                   what=what,gridlines=gridlines,
                                   col=col,new=new,
                                   breaks=breaks,colorbar=colorbar,...) else
  if (projection=="sp") map2sphere(X,lonR=lonR,latR=-90,axiR=axiR,
                                   what=what,gridlines=gridlines,
                                   col=col,new=new,
                                   breaks=breaks,colorbar=colorbar,...)
  invisible(X)
}


map.corfield <- function(x,new=TRUE,
                         xlim=NULL,ylim=NULL,zlim=NULL,n=15,
                         projection="lonlat",
                         col=NULL,breaks=seq(-1,1,0.1),
                         lonR=NULL,latR=NULL,axiR=0,what=c("fill","contour"),
                         gridlines=TRUE,verbose=FALSE,...) {
  if (verbose) print("map.corfield")
  stopifnot(inherits(x,'corfield'))
  x <- subset(x,it=it,is=is)
  projection <- tolower(projection)
  dim(x) <- attr(x,'dimensions')[1:2]
  
  ## if zlim is specified, then mask data outside this range
  if (!is.null(zlim)) {
    d <- dim(x)
    mask <- (x < min(zlim)) | (x > max(zlim))
    x[mask] <- NA
    dim(x) <- d
    if (verbose) {print(zlim); print(dim(x)); print(sum(mask))}
  }

  if (projection=="lonlat") lonlatprojection(x=x,n=n,col=col,breaks=breaks,
                             what=what,gridlines=gridlines,new=new,...) else
  if (projection=="sphere") map2sphere(x=x,lonR=lonR,latR=latR,axiR=axiR,
                                       what=what,gridlines=gridlines,
                                       col=col,new=new,...) else
  if (projection=="np") map2sphere(x,lonR=lonR,latR=90,axiR=axiR,
                                   what=what,gridlines=gridlines,
                                   col=col,new=new,...) else
  if (projection=="sp") map2sphere(x,lonR=lonR,latR=-90,axiR=axiR,
                                   what=what,gridlines=gridlines,
                                   col=col,new=new,...)
  if (!is.null(attr(x,'x.longitude')) & !is.null(attr(x,'x.latitude')))
      points(attr(x,'x.longitude'),attr(x,'x.latitude'),lwd=2,cex=1.2)
  invisible(x)
}


map.trend <- function(x,it=NULL,is=NULL,new=TRUE,
                      xlim=NULL,ylim=NULL,zlim=NULL,n=15,
                      projection="lonlat",
                      col=NULL,breaks=NULL,
                     lonR=NULL,latR=NULL,axiR=0,what=c("fill","contour"),
                     gridlines=TRUE,verbose=FALSE,...) {
  if (verbose) print('map.trend')
  stopifnot(inherits(x,'field'),inherits(x,'trend'))
  x <- subset(x,it=it,is=is)
  projection <- tolower(projection)
  X <- attr(x,'pattern')

  ## if zlim is specified, then mask data outside this range
  if (!is.null(zlim)) {
    d <- dim(x)
    mask <- (x < min(zlim)) | (x > max(zlim))
    x[mask] <- NA
    dim(x) <- d
    if (verbose) {print(zlim); print(dim(x)); print(sum(mask))}
  } 
  attr(X,'longitude') <- attr(x,'longitude')
  attr(X,'latitude') <- attr(x,'latitude')
  attr(X,'variable') <- paste(attr(x,'variable'),'trend')
  attr(X,'time') <- range(index(x))
  attr(X,'unit') <- paste('d',attr(x,'unit'),'/decade')
  attr(X,'source') <- attr(x,'source')
  dim(X) <- attr(x,'dimension')[1:2]
  #str(X)
  if (projection=="lonlat") lonlatprojection(x=X,n=n,col=col,breaks=breaks,
                             what=what,gridlines=gridlines,new=new,...) else
  if (projection=="sphere") map2sphere(x=X,lonR=lonR,latR=latR,axiR=axiR,
                                       what=what,gridlines=gridlines,
                                       col=col,new=new,...) else
  if (projection=="np") map2sphere(X,lonR=lonR,latR=90,axiR=axiR,
                                   what=what,gridlines=gridlines,
                                   col=col,new=new,...) else
  if (projection=="sp") map2sphere(X,lonR=lonR,latR=-90,axiR=axiR,
                                   what=what,gridlines=gridlines,
                                   col=col,new=new,...)
  invisible(X)
}


lonlatprojection <- function(x,it=NULL,is=NULL,xlim=NULL,ylim=NULL,
                             n=15,col=NULL,breaks=NULL,geography=TRUE,
                             what=c("fill","contour"),gridlines=TRUE,
                             new=TRUE,colorbar=NULL,verbose=FALSE,...) {
  if (verbose) print('lonlatprojection')
  #print(dim(x));
  #print(c(length(attr(x,'longitude')),length(attr(x,'latitude'))))
  data("geoborders",envir=environment())
  if(sum(is.finite(x))==0) stop('No valid data')
  # To deal with grid-conventions going from north-to-south or east-to-west:
  srtx <- order(attr(x,'longitude')); lon <- lon(x)[srtx]
  srty <- order(attr(x,'latitude'));  lat <- lat(x)[srty]
  #print('meta-stuff')
  unit <- unit(x); variable <- varid(x); varid <- varid(x)
  if ( (unit=="degC") | (unit=="deg C") | (unit=="degree C") )
    unit <- "degree*C"
  if (unit=="%") unit <- "'%'"
  if ( (tolower(variable)=="t(2m)") | (tolower(variable)=="t2m") |
      (tolower(variable)=="2t") )
    variable <- "T[2*m]"
#  if (inherits(x,'corfield'))
#    main=eval(parse(text=paste('expression(paste("correlation: ",',
#                               variable," *(",unit,")))",sep=""))) else
  main=eval(parse(text=paste('expression(',variable," *(",unit,"))",sep="")))
  sub <- attr(x,'source')
  colid <- 't2m'; if (is.precip(x)) colid <- 'precip'
  if (is.null(colorbar)) colorbar <- (sum(is.element(what,'fill')>0))
                                      
  #print('time')
  if (!is.null(attr(x,'timescale'))) {
    #print(attr(x,'timescale'))
    timescale <- attr(x,'timescale')
    if (timescale == 'annual') {
      t1 <- year(attr(x,'time'))[1]
      t2 <- year(attr(x,'time'))[2]
    } else
    if (sum(is.element(c('month','season'),timescale))>0) {
#      t1 <- format(attr(x,'time')[1],'%Y-%m')
#      t2 <- format(attr(x,'time')[2],'%Y-%m')
      t1 <- paste(year(attr(x,'time'))[1],month(attr(x,'time'))[1])
      t2 <- paste(year(attr(x,'time'))[2],month(attr(x,'time'))[2])
    } else {
      t1 <- attr(x,'time')[1]  
      t2 <- attr(x,'time')[2]
    }
    period <- paste('[',t1,', ',t2,']',sep='')
  } else period <- NULL
  #print(period)
  method <- attr(x,'method')
  x <- x[srtx,srty]
  #print("HERE"); print(xlim); str(x)
  if (!is.null(xlim)) {
    outside <- (lon < xlim[1]) | (lon > xlim[2])
    x[outside,] <- NA
  } else xlim <- range(lon)
  
  if (!is.null(ylim)) {
    outside <- (lat < ylim[1]) | (lat > ylim[2])
    x[,outside] <- NA
  } else ylim=range(lat)

  if (!is.null(col))
    if (length(col)>1) n <- length(col)
  #print(n); print(summary(c(x)))
  if (is.null(breaks))
    breaks <- pretty(c(x),n=n)
  #print(breaks)

  if (is.null(col)) col <- colscal(n=length(breaks)-1,col=colid) else
  if (length(col)==1) {
     palette <- col
     col <- colscal(col=palette,n=length(breaks)-1)
  }
  if (length(breaks) != length(col)+1)
    breaks <- seq(min(c(x),na.rm=TRUE),max(c(x),na.rm=TRUE),
                  length=length(col)+1)
                      
  #print(variable)
#  if ( (tolower(variable)=='precip') | (tolower(variable)=='tp') )
   if (colid=='precip') col <- rev(col)
  
  #print(c(length(breaks),length(col)))
  #if (is.Date(what))

  if ( (par()$mfcol[1]> 1) | (par()$mfcol[2]> 1) ) new <- FALSE
      
  if (new) {
    #dev.new()
    par(bty="n",xaxt="n",yaxt="n",xpd=FALSE,
        fig=c(0.05,0.95,0.13,0.95),mar=rep(1,4))
#    par(bty="n",xaxt="n",yaxt="n",xpd=FALSE,
#        fig=c(0.05,0.95,0.12,0.95))
  } else {
    par(bty="n",xaxt="n",yaxt="n",xpd=FALSE,mar=rep(1,4))
#    par(bty="n",xaxt="n",yaxt="n",xpd=FALSE)
  }

  plot(range(lon),range(lat),type="n",xlab="",ylab="", # REB 10.03
         xlim=xlim,ylim=ylim)                # to sumerimpose.
  par0 <- par()
  if (sum(is.element(tolower(what),'fill'))>0)   
    image(lon,lat,x,xlab="",ylab="",add=TRUE,
          col=col,breaks=breaks,xlim=xlim,ylim=ylim,...)
  
  if (geography) {
    lines(geoborders$x,geoborders$y,col="darkblue")
    lines(attr(geoborders,'borders')$x,attr(geoborders,'borders')$y,col="pink")
    lines(geoborders$x+360,geoborders$y,col="darkblue")
  }
  if (sum(is.element(tolower(what),'contour'))>0)
     contour(lon,lat,x,lwd=1,col="grey70",add=TRUE)
  if (gridlines) grid()
  par(xpd=TRUE)
  dlat <- diff(range(lat))/60
  #print(dlat)
  text(lon[1],lat[length(lat)] + dlat,main,pos=4,font=2)
  text(lon[1],lat[1] - dlat,sub,col="grey30",pos=4,cex=0.7)

  if (!is.null(period))
    text(lon[length(lon)],lat[length(lat)] + dlat,period,pos=2,cex=0.7,col="grey30")
  if (!is.null(method))
    text(lon[length(lon)],lat[1] - 0.5*dlat,method,col="grey30",pos=2,cex=0.7)
  if (colorbar) {
    par(xaxt="s",fig=c(0.05,0.95,0.01,1))
    breaks <- round(seq(min(x,na.rm=TRUE),max(x,na.rm=TRUE),length=length(col)),1)
    image.plot(horizontal=TRUE,legend.only=TRUE,zlim=range(x,na.rm=TRUE),
               lab.breaks=breaks,col=col,axis.args=list(cex.axis=0.8),
               border=FALSE)
#    image.plot(horizontal=TRUE,legend.only=TRUE,zlim=range(x,na.rm=TRUE),
#               lab.breaks=pretty(x),col=col,legend.width=1,axis.args=list(cex.axis=0.8),
#               border=FALSE,add=TRUE,graphics.reset=TRUE)
#    par(fig = c(0.3, 0.7, 0.05, 0.10),mar=rep(0,4),cex=0.8,
#        new = TRUE, mar=c(1,0,0,0), xaxt = "s",yaxt = "n",bty = "n")
#     print("colourbar"); print(pretty(x))
#    bar <- cbind(breaks,breaks)
#    image(breaks,c(1,2),bar,col=col,breaks=breaks)
  
#    par(bty="n",xaxt="n",yaxt="n",xpd=FALSE,
#        fig=c(0.05,0.95,0.12,0.95),new=TRUE)
#    plot(range(lon),range(lat),type="n",xlab="",ylab="")
#    par(fig=par0$fig,new = TRUE, mar=par0$mar, xaxt = "n",yaxt = "n",bty = 'n')
#    plot(range(lon),range(lat),type="n",xlab="",ylab="",
#         xlim=xlim,ylim=ylim)
    #browser()
    par(fig=par0$fig,mar=par0$mar,new=TRUE,xaxt="n")
    plot(range(lon),range(lat),type="n",xlab="",ylab="", # REB 10.03
         xlim=xlim,ylim=ylim)                # to sumerimpose.
  }
  par(xaxt="s",yaxt="s",las=1,col.axis='grey',col.lab='grey',cex.lab=0.7,cex.axis=0.7)
  axis(2,at=pretty(lat(x)),col='grey')
  axis(3,at=pretty(lon(x)),col='grey')
  grid()

  par(col.axis='black',col.lab='black',cex.lab=1,cex.axis=1)
  result <- list(x=lon,y=lat,z=x,breaks=breaks)
  invisible(result)
}


#map.pca <- function(x,it=NULL,is=NULL,ipca=1,cex=1.5,xlim=NULL,ylim=NULL,
#                    n=100,projection="lonlat",
#                    FUN=NULL,col=NULL,breaks=NULL,
#                    lonR=NULL,latR=NULL,axiR=0,what=c("fill","contour"),
#                    gridlines=TRUE) {
#  print('map.pca')
#  x <- subset(x,it=it,is=is)
#  lon <- attr(x,'longitude')
#  lat <- attr(x,'latitude')
#
#  U <- attr(x,'pattern')
#  L <- attr(x,'eigenvalues')
#  V <- coredata(x)
#  #print(dim(U)); print(length(L)); print(dim(V))
  #X <- U %*% diag(L) %*% t(V)
#  X <- U %*% diag(L)
  #print(dim(X))
#  X <- X[,ipca]
#  X <- apply(X,1,FUN,na.rm=TRUE);
#  N <- length(X)
#  col <- colscal(n=n)
#  ax <- quantile(abs(X),0.95,na.rm=TRUE)
#  scale0 <- seq(-ax,ax,length=n)
  #print(scale0); print(n)
#  a.T <- rep(NA,N)
#  for (i in 1:N) a.T[i] <-  sum(X[i] > scale0)
#  a.T[a.T < 1] <- 1; a.T[a.T > 100] <- 100
#
#  par(bty="n",xaxt="n",yaxt="n",xpd=FALSE,
#      fig=c(0.05,0.95,0.1,0.95),mar=rep(1,4))
  #print(rbind(lon,lat,a.T))
#  plot(lon,lat,col=col[a.T],pch=19,xlab="",ylab="",cex=cex)
#  data(geoborders,envir=environment())
#  lines(geoborders)
#  lines(geoborders$x - 360,geoborders$y)
#  colbar(scale0,col,fig=c(0.90,0.95,0.05,0.25))
#}

map.pca <- function(x,pattern=1,new=TRUE,FUN='mean',
                    col=NULL,cex=1.5,xlim=NULL,ylim=NULL,zlim=NULL,
                    colbar=NULL,verbose=FALSE,...) {
    ##
    #args <- list(...)
    #print(args)
    X <- rbind(attr(x,'pattern')[,pattern],attr(x,'pattern')[,pattern])
  #print(dim(X))
  #str(x)
  X <- attrcp(x,X)

  ## if zlim is specified, then mask data outside this range
  if (!is.null(zlim)) {
    d <- dim(X)
    mask <- (X < min(zlim)) | (X > max(zlim))
    X[mask] <- NA
    dim(X) <- d
    if (verbose) {print(zlim); print(dim(X)); print(sum(mask))}
  }    
  attr(X,'longitude') <- lon(x)
  attr(X,'latitude') <- lat(x)
  class(X) <- 'station'
  if (is.null(col)) {
    col <- colscal(30,col=varid(x))
  }
  
  map.station(X,new=new,FUN=FUN,col=col,bg=col,
              colbar=list(col=col,type='r',v=0),
              cex=cex,xlim=xlim,ylim=ylim,...)
}

map.mvr <- function(x,it=NULL,is=NULL,new=TRUE,xlim=NULL,ylim=NULL,
                    n=15,projection="lonlat",
                    col=NULL,breaks=NULL,
                    lonR=NULL,latR=NULL,axiR=0,what=c("fill","contour"),
                    gridlines=TRUE,verbose=FALSE,...) {
  x <- subset(x,it=it,is=is)
  
}

map.cca <- function(x,it=NULL,is=NULL,new=TRUE,icca=1,
                    xlim=NULL,ylim=NULL,zlim=NULL,
                    what=c("fill","contour"),cex=1.5,
                    n=15,projection="lonlat",
                    lonR=NULL,latR=NULL,colorbar=TRUE,
                    axiR=0,gridlines=TRUE,col=NULL,
                    breaks=NULL,verbose=FALSE,...) {
  #print('map.cca')
  #x <- subset(x,it=it,is=is)
  ## browser()
  # For plotting, keep the same kind of object, but replace the patterns in
  # the eof/pca with the CCA patterns
  Y <- x$Y
  #print(dim(attr(Y,'pattern'))); print(dim(U))
  #attr(Y,'pattern') <- U
  U <- x$B.m
  dim(U) <- c(dim(attr(Y,'pattern'))[-length(dim(attr(Y,'pattern')))],
                                     length(x$i.eofs))
  attr(Y,'pattern') <- U
  attr(Y,'eigenvalues') <- rep(1,length(x$i.eofs))
  attr(Y,'time') <- range(index(x))
  X <- x$X
  #print(dim(attr(X,'pattern'))); print(dim(V))
  #attr(X,'pattern') <- V
  V <- x$A.m
  dim(V) <- c(dim(attr(X,'pattern'))[-length(dim(attr(X,'pattern')))],
                                     length(x$i.eofs))
  attr(X,'pattern') <- V
  attr(X,'eigenvalues') <- rep(1,length(x$i.eofs))
  attr(X,'time') <- range(index(x))

  # REB removed '...' in the two following map calls.
#  map(Y,icca,xlim=xlim,ylim=ylim,what=what,
#      projection=projection,lonR=lonR,latR=latR,axiR=axiR,
#      gridlines=gridlines,col=col,breaks=breaks,FUN='mean')
#  dev.new()
#  map(X,icca,xlim=xlim,ylim=ylim,what=what,
#      projection=projection,lonR=lonR,latR=latR,axiR=axiR,
#      gridlines=gridlines,col=col,breaks=breaks,FUN='mean')
#  print('Need to fix breaks and map.station')

#  col <- rgb( c(rep(0,15),1-sqrt(seq(0,1,length=15))),
#              abs(sin(seq(0,pi,length=30))),
#              c(sqrt(seq(0,1,length=15)),rep(1,15)) )
#  col <- colscal(30,col=varid(x))
#  if (is.precip(X)) col.x <- rev(col) else
#                    col.x <- col
#  if (is.precip(Y)) col.y <- rev(col) else
#                    col.y <- col
# REB: removed col=col.y,bg=col.y

  if (sum(is.element(what,'map'))>0)
    par(fig=c(0,0.5,0.5,1)) ## mar=c(0.05,.05,0.05,0.05),
  else 
    par(fig=c(0,0.5,0.5,1),mar=c(0.2,.2,0.2,0.2))

  
  map(Y,pattern=icca,xlim=xlim,ylim=ylim,what=what,cex=cex,
      projection=projection,lonR=lonR,latR=latR,axiR=axiR,
      gridlines=gridlines,FUN='mean',colorbar=colorbar,
      colbar=list(col=col,type='r',v=0),
      showall=FALSE,new=FALSE)
  ## browser()
  if (sum(is.element(what,'ts'))>0)
      par(fig=c(0,1,0.5,1),new=TRUE) else
  par(fig=c(0.5,1,0.5,1),new=TRUE) ## mar=c(0,0,0,0),
  map(X,pattern=icca,xlim=xlim,ylim=ylim,what=what,cex=cex,
      projection=projection,lonR=lonR,latR=latR,axiR=axiR,
      gridlines=gridlines,FUN='mean',colorbar=colorbar,
      colbar=list(col=col,type='r',v=0),
      showall=FALSE,new=FALSE)
  
  invisible(list(U=U,V=V))
}


# Produce a KMZ-file to show the data in GoogleEarth.
map.googleearth <- function(x) {
}
