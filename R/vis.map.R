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
