#' Plot maps for esd objects
#' 
#' Make map of geophysical data. These plot functions are S3 methods for esd
#' objects.
#' 
#' @aliases map map.default map.matrix map.data.frame map.station
#' map.stationmeta map.stationsummary map.comb map.eof map.ds
#' map.field map.corfield map.cca map.events map.mvr map.pca map.array
#' map.trend map.mvcomb lonlatprojection rotM map2sphere
#' @seealso map.trajectory vec map.dsensemble
#'
#' @param x the object to be plotted; in \code{rotM}, x holds a vector of
#' x-coordinates.
#' @param FUN The function to be applied on x before mapping (e.g. \code{mean})
#' @param colbar The colour scales defined through \code{colscal}. Users can
#' specify the colour `pal'*ette (`pal'), the number of breaks (`n'), values of
#' `breaks', and the position of the color bar from the main plot (`pos'). The
#' `rev' argument, will produce a reversed color bar if set to TRUE. The other
#' arguments (`type',`h' and `v') are more specific to \code{col.bar} and are
#' used only if argument `fancy' is set to TRUE (not yet finished). \code{colbar=NULL}
#' is used if the colourbar is not to be shown. Also use \code{colbar=NULL} to present several 
#' maps in one figure (e.g. with \code{par(mfcol=c(2,2))}).
#' @param lab \code{'default'} to show a lable saying what variable (unit) and time period. 
#' \code{'simple'} to just use \code{varid(x)}, and \code{'unit'} to show variable and unit. 
#' Other strings will be used as the label for the plot and \code{NULL} shows no lable.
#'
#' @param it see \code{\link{subset}}
#' @param is see \code{\link{subset}}
#' @param new TRUE: create a new graphic device.
#' @param projection Projections: c("lonlat","sphere","np","sp") - the latter
#' gives stereographic views from the North and south poles.
#' @param xlim see \code{\link{plot}} - only used for 'lonlat' and 'sphere'
#' projections.
#' @param ylim see \code{\link{plot}} - only used for 'lonlat' and 'sphere'
#' projections.
#' @param n The number of colour breaks in the color bar
#' @param breaks graphics setting - see \code{\link{image}}
#' @param type graphics setting - colour shading or contour
#' @param gridlines Only for the lon-lat projection
#' @param lonR Only for the spherical projection used by \code{map2sphere} to change viewing angle
#' @param latR Only for the spherical projection used by \code{map2sphere} to change viewing angle
#' @param axiR Only for the spherical projection used by \code{map2sphere} to change viewing angle
#' @param style Only for the spherical projection used by \code{map2sphere} to apply night shade effect. c('plain','night')
#' @param y a vector of y coordinates
#' @param z a vector of z coordinates
#' @param ip Selects which pattern (see \code{\link{EOF}}, \code{\link{CCA}}) to plot
#' @param geography TRUE: plot geographical features
#' @param angle for hatching
#' @param a used in \code{\link{vec}} to scale the length of the arrows
#' @param r used in \code{\link{vec}} to make a 3D effect of plotting the arrows up in the air.
#' @param ix used to subset points for plotting errors
#' @param iy used to subset points for plotting errors
#' @param colorbar Show the color bar in the map (default TRUE). If FALSE, the
#' colorbar is not added into the map (ignored).
#' @param cex Size of symbols.
#' @param cex0 Scaling of symbols if \code{cex} is defined by a variable with
#' different size for different locations.
#' @param cex.subset ...
#' @param add.text Add abbreviated location names.
#' @param full.names Show the full name of the location.
#' @param showall Default is set to FALSE
#' @param showaxis If set to FALSE, the axis are not displayed in the plot.
#' @param fancy If set to true, will use \code{\link{col.bar}} instead of
#' \code{image} to produce the colour bar
#' @param text If TRUE, display text info on the map.The default is set to
#' FALSE
#' @param show.val Display the values of 'x' or 'FUN(x)' on top of the coloured
#' map.
#' @param legend.shrink If set, the size of the color bar is shrinked (values
#' between 0 and 1)
#' @param ... further arguments passed to or from other methods.
#' @param land if TRUE mask land, else mask ocean
#' @param what What to map: ['eof','field] for EOF pattern or the field
#' @param type - default c('fill','contour')
#' recovered from the EOFs.
#' @param fig see \code{\link{par}}
#' @param add set add=TRUE if adding the map as a subplot
#' @param nbins number of bins/colour categories
#'
#' @return A field object
#' 
#' @seealso \code{\link{plot.station}}, \code{\link{showmaps}}
#' @keywords map
#' @examples
#' 
#' # Select stations in ss and map the geographical location 
#' # of the selected stations with a zoom on Norway.
#' ss <- select.station(cntr="NORWAY",param="precip",src="GHCND")
#' map(ss, col="blue",bg="lightblue",xlim = c(-10,30) , ylim = c(50,70), new=FALSE)
#' 
#' ## Get NACD data and map the mean values
#' y <- station.nacd()
#' map(y,FUN='mean',colbar=list(pal="t2m",n=10), cex=2, new=FALSE)
#' 
#' # Examples of cyclone maps (map.events)
#' data(storms)
#' # Subset cyclones from the start of January 2016 lasting a minimum of 10 time steps 
#' x <- subset(storms,it=c("2016-01-01","2016-01-07"),ic=list(param="trackcount",pmin=10))
#' # Map with points and lines showing the cyclone centers and trajectories
#' map(x, type=c("trajectory","points"), col="blue")
#' ## Map with only the trajectory and start and end points
#' map(x, type=c("trajectory","start","end"), col="red")
#' ## Map showing the cyclone depth (slp at center) as a color scale (rd = red scale)
#' map(x, param="pcent", type=c('trajectory','start'), 
#'     colbar=list(pal="rd", rev=TRUE, breaks=seq(980,1000,2)), 
#'     alpha=0.9, new=FALSE)
#' 
#' # Example showing two maps in one figure:
#' \dontrun{
#' rr <- retrieve('~/data/data.ECAD/rr_ens_mean_0.1deg_reg.nc',it=c(2000,2001))
#' it <- 195
#' par(mfcol=c(1,2))
#' attr(rr,'variable') <- 'precip'
#' map(rr,FUN='mean',new=FALSE,type='fill',colbar=NULL,lab=paste(range(index(rr)),collapse=' - '))
#' map(subset(rr,it=it),FUN='mean',new=FALSE,type='fill',colbar=NULL,lab=as.character(index(rr)[it]))
#' }


#' @export map
map <- function(x,...) UseMethod("map")

#' @exportS3Method
#' @export
map.default <- function(x,...,FUN='mean',it=NULL,is=NULL,new=FALSE,
                        projection="lonlat",xlim=NULL,ylim=NULL,zlim=NULL,lab='default',
                        colbar= list(pal=NULL,rev=FALSE,n=10,breaks=NULL,pos=0.05,
                                     show=TRUE,type="p",cex=2,h=0.6,v=1),
                        type=c("fill","contour"),gridlines=FALSE,cex=2,
                        lonR=NULL,latR=NULL,axiR=NULL,style='plain',
                        verbose=FALSE,plot=TRUE,add=FALSE) {
  
  
  ## default with no arguments will produce a map showing available station
  ## data in the esd package.
  
  if (verbose) print('map.default')
  #def.par <- par(no.readonly = TRUE) # save default, for resetting...
  if (is.logical(colbar)) colbar <- NULL
  ## If only a few items are provided in colbar - then set the rest to the default
  if (!is.null(colbar)) {
    colbar <- colbar.ini(x,FUN=FUN,colbar=colbar,verbose=FALSE)
  }
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
  if (attr(X,'unit') =='%') attr(X,'unit') <- "'%'"
  attr(X,'source') <- attr(x,'source')
  attr(X,'variable') <- varid(x)
  if (inherits(X,'zoo')) {
    attr(X,'time') <- range(index(x))
  } else if (!is.null(attr(x,'time'))) {
    attr(X,'time') <- attr(x,'time')
  }
  if (plot) {
    if (projection=="lonlat") {
      z <- lonlatprojection(x=X,xlim=xlim,ylim=ylim,colbar=colbar,verbose=verbose,
                            type=type,new=new,gridlines=gridlines,...)
    } else if (projection=="sphere") {
      z <- map2sphere(x=X,lonR=lonR,latR=latR,axiR=axiR,xlim=xlim,ylim=ylim,
                      type=type,gridlines=gridlines,colbar=colbar,new=new,...)
    } else if (projection=="np") {
      z <- map2sphere(X,lonR=lonR,latR=90,axiR=axiR,xlim=xlim,ylim=ylim,
                      type=type,gridlines=gridlines,colbar=colbar,new=new,...)
    } else if (projection=="sp") {
      z <- map2sphere(X,lonR=lonR,latR=-90,axiR=axiR,new=new,xlim=xlim,ylim=ylim,
                      type=type,gridlines=gridlines,colbar=colbar,...)
    } 
  } else z <- X
  invisible(z)
}

#' @exportS3Method
#' @export
map.matrix <- function(x,...,it=NULL,is=NULL,new=FALSE,projection="lonlat",
                       xlim=NULL,ylim=NULL,zlim=NULL,n=15,
                       colbar= list(pal=NULL,rev=FALSE,n=10,breaks=NULL,
                                    pos=0.05,show=TRUE,type="p",cex=2,h=0.6,v=1),
                       type=c("fill","contour"),gridlines=FALSE,
                       lonR=NULL,latR=NULL,axiR=NULL,verbose=FALSE,
                       ip=1,plot=TRUE) {
  
  ## If x is provided, map only x...
  ## default with no arguments will produce a map showing the station data in the esd package.
  ##  image(lon(x),lat(x),x)
  
  if (verbose) print('map.matrix')
  #def.par <- par(no.readonly = TRUE) # save default, for resetting...
  if (!is.null(is)) x <- subset(x,is=is)  # if is is set, then call subset
  if (inherits(x,'zoo')) attr(x,'time') <- range(index(x))
  if (verbose) str(x)
  if (plot) {
    if (projection=="lonlat") {
      z <- lonlatprojection(x=x,new=new,xlim=xlim,ylim=ylim,zlim=zlim,colbar=colbar,
                            type=type,gridlines=gridlines,verbose=verbose,...)
    } else if (projection=="sphere") {
      z <- map2sphere(x=x,new=new,xlim=xlim,ylim=ylim,zlim=zlim,colbar=colbar,
                      lonR=lonR,latR=latR,axiR=axiR,verbose=verbose,...)
    } else if (projection=="np") {
      z <- map2sphere(x,new=new,xlim=xlim,ylim=ylim,zlim=zlim,lonR=lonR,latR=90,
                      colbar=colbar,verbose=verbose,...)
    } else if (projection=="sp") {
      z <- map2sphere(x,new=new,xlim=xlim,ylim=ylim,zlim=zlim,lonR=lonR,latR=-90,
                      colbar=colbar,verbose=verbose,...)
    }
  }
  invisible(z)
  #map.station(NULL,...)
}

#' @exportS3Method
#' @export
map.data.frame <- function(x,...,it=NULL,is=NULL,new=FALSE,projection="lonlat",
                           xlim=NULL,ylim=NULL,zlim=NULL,n=15,
                           colbar= list(pal=NULL,rev=FALSE,n=10,breaks=NULL,
                                        pos=0.05,show=TRUE,type="p",cex=2,h=0.6,v=1),
                           type=c("fill","contour"),gridlines=FALSE,
                           lonR=NULL,latR=NULL,axiR=NULL,verbose=FALSE,
                           ip=1,plot=TRUE) {
  if (verbose) print('map.data.frame')
  #def.par <- par(no.readonly = TRUE) # save default, for resetting...
  attr(x,'location') <- x$location; x$location <- NULL
  attr(x,'longitude') <- x$longitude; x$longitude <- NULL
  attr(x,'latitude') <- x$latitude; x$latitude <- NULL
  attr(x,'country') <- x$country; x$country <- NULL
  x <- as.matrix(x)
  z <- map(x,it=it,is=is,new=new,projection=projection,
           xlim=xlim,ylim=ylim,zlim=zlim,n=15,
           colbar= colbar,type=type,gridlines=gridlines,
           lonR=lonR,latR=latR,axiR=axiR,verbose=verbose,
           ip=ip,plot=plot,...)
  invisible(z)
}

#' @exportS3Method
#' @export
map.array <- function(x,...,FUN='mean',ip=NULL,is=NULL,new=FALSE,
                      projection="lonlat",na.rm=TRUE,
                      xlim=NULL,ylim=NULL,zlim=NULL,##n=15,
                      colbar=list(col=NULL,rev=FALSE,breaks=NULL,pos=0.05,
                                  show=TRUE,type="r",cex=2,h=0.6,v=1),
                      type=c("fill","contour"),gridlines=FALSE,
                      lonR=NULL,latR=NULL,axiR=NULL,verbose=FALSE,plot=TRUE) {
  if (verbose) print('map.array')
  #def.par <- par(no.readonly = TRUE) # save default, for resetting...
  if (!is.null(is)) x <- subset(x,is=is)  # if is is set, then call subset
  if (is.null(ip)) {
    ## If it is NULL, then aggregate all of 3rd dimension
    D <- dim(x)
    x2d <- x
    dim(x2d) <- c(D[1]*D[2],D[3])
    z <- apply(x2d,1,FUN,na.rm=na.rm)
    z <- as.matrix(z)
    dim(z) <- c(D[1],D[2])
    str(z)
  } else  z <- x[,,ip]
  d <- dim(z)
  
  ## if it is a vector of indices aggregate the selected indices
  if (length(d)==3) {
    dim(z) <- c(d[1]*d[2],d[3])
    z <- apply(z,2,FUN)
    dim(z) <- c(d[1],d[2])
  }
  attr(z,'longitude') <- lon(x)
  attr(z,'latitude') <- lat(x)
  attr(z,'variable') <- varid(x)
  attr(z,'unit') <- attr(x,'unit')[1]
  attr(z,'colbar') <- colbar
  
  if (plot) map(z,new=new,xlim=xlim,ylim=ylim,zlim=zlim,colbar=colbar,
                lonR=lonR,latR=latR,axiR=axiR,
                type=type,gridlines=gridlines,projection=projection,
                verbose=verbose,...)
  invisible(z)
}

#' @exportS3Method
#' @export
map.comb <- function(x,...,it=NULL,is=NULL,new=FALSE,projection="lonlat",
                     xlim=NULL,ylim=NULL,zlim=NULL,#n=15,
                     colbar=list(pal=NULL,rev=FALSE,n=10,breaks=NULL,
                                 pos=0.05,show=TRUE,type="p",cex=2,h=0.6,v=1),
                     type=c("fill","contour"),gridlines=FALSE,
                     lonR=NULL,latR=NULL,axiR=NULL,verbose=FALSE,
                     ip=1,plot=TRUE) {
  if (verbose) print('map.comb')
  #def.par <- par(no.readonly = TRUE) # save default, for resetting...
  stopifnot(inherits(x,'eof'))
  x <- subset(x,it=it,is=is)
  projection <- tolower(projection)
  ## if (is.null(col)) col <- colscal(pal=colbar$pal,n=n-1,rev=colbar$rev) else
  ## if (length(col)==1) {
  ##   pal <- col
  ##   col <- colscal(pal=pal,n=n-1,rev=colbar$rev)
  ## }
  if (is.null(varid(x))) attr(x,'variable') <- 'NA'
  ## if (tolower(varid(x))=='precip') col <- rev(col) 
  
  z <- map.eof(x=x,xlim=xlim,ylim=ylim,zlim=zlim,ip=ip,
               projection=projection,colbar=colbar,new=new,
               lonR=lonR,latR=latR,axiR=axiR,type=type,
               gridlines=gridlines,verbose=verbose,plot=plot,...) -> result
  invisible(z)
}


#' @exportS3Method
#' @export
map.eof <- function(x,...,it=NULL,is=NULL,new=FALSE,projection="lonlat",what="eof",
                    xlim=NULL,ylim=NULL,zlim=NULL,lab="default",
                    colbar=list(pal=NULL,rev=FALSE,n=10,breaks=NULL,
                                pos=0.05,show=TRUE,type="p",cex=2,h=0.6,v=1),
                    type=c("fill","contour"),gridlines=FALSE,
                    lonR=NULL,latR=NULL,axiR=NULL,verbose=FALSE,
                    ip=1,cex=1,plot=TRUE) {
  
  if (verbose) print('map.eof')
  #def.par <- par(no.readonly = TRUE) # save default, for resetting...
  stopifnot(inherits(x,'eof'))
  ##x <- subset(x,it=it,is=is)
  projection <- tolower(projection)
  
  ## REB 2016-10-19: one option is to recover the field and then maps the field
  if (what=="field") {
    if (verbose) print('what=field: recover the field before mapping')
    x <- subset(x,it=it,is=is)
    if (verbose) {print(dim(x)); range(index(x)); range(lon(x)); range(lat(x))}
    y <- as.field(x)
    z <- map(y,new=new,projection=projection,xlim=xlim,ylim=ylim,
             zlim=zlim,colbar=colbar,type=type,gridlines=gridlines,
             lonR=lonR,latR=latR,axiR=axiR,verbose=verbose,cex=cex,plot=plot)
    invisible(z)
  } else {
    
    tot.var <- attr(x,'tot.var')
    D <- attr(x,'eigenvalues')
    var.eof <- 100* D^2/tot.var
    d.eof <- dim(attr(x,'pattern'))
    if (verbose) {print('Dimensions of pattern'); print(d.eof)}
    if (length(d.eof)==2) dim(attr(x,'pattern')) <- c(d.eof,1)
    X <- attr(x,'pattern')[,,ip]
    
    ## if zlim is specified, then mask data outside this range
    if (!is.null(zlim)) {
      d <- dim(X)
      mask <- (X < min(zlim)) | (X > max(zlim))
      X[mask] <- NA
      dim(X) <- d
      if (verbose) {print(zlim); print(dim(X)); print(sum(mask))}
    }
    ##str(x)
    attr(X,'longitude') <- attr(x,'longitude')
    attr(X,'latitude') <- attr(x,'latitude')
    attr(X,'variable') <- attr(x,'variable')
    attr(X,'unit') <- attr(x,'unit')[1]
    if (attr(X,'unit') =='%') attr(X,'unit') <- "'%'"
    attr(X,'source') <- attr(x,'source')
    attr(X,'time') <- range(index(x))
    attr(X,'greenwich') <- attr(x,"greenwich")
    attr(X,'colbar') <- colbar
    
    if ( (ip==1) & !is.null(attr(x, "area.mean.expl")) )
      if (attr(x, "area.mean.expl"))
        type <- "fill"
    if (plot) {
      if (projection=="lonlat") {
        z <- lonlatprojection(x=X,it=it,xlim=xlim,ylim=ylim,lab=lab,
                              colbar=colbar,new=new,type=type,
                              gridlines=gridlines,verbose=verbose,...)
      } else if (projection=="sphere") {
        z <- map2sphere(x=X,it=it,lonR=lonR,latR=latR,axiR=axiR,lab=lab,
                        xlim=xlim,ylim=ylim,type=type,gridlines=gridlines,
                        colbar=colbar,new=new,verbose=verbose,...)
      } else if (projection=="np") {
        z <- map2sphere(X,it=it,lonR=lonR,latR=90,axiR=axiR,lab=lab,
                        xlim=xlim,ylim=ylim,type=type,gridlines=gridlines,
                        colbar=colbar,new=new,verbose=verbose,...)
      } else if (projection=="sp") {
        z <- map2sphere(X,it=it,lonR=lonR,latR=-90,axiR=axiR,lab=lab,
                        xlim=xlim,ylim=ylim,type=type,gridlines=gridlines,
                        colbar=colbar,new=new,verbose=verbose,...)
      }
    } else z <- X
  }
  invisible(z)
}

#' @exportS3Method
#' @export
map.ds <- function(x,...,it=NULL,is=NULL,new=FALSE,projection="lonlat",
                   xlim=NULL,ylim=NULL,zlim=NULL,
                   colbar=list(pal=NULL,rev=FALSE,n=10,breaks=NULL,
                               pos=0.05,show=TRUE,type="p",cex=2,h=0.6,v=1),
                   type=c("fill","contour"),gridlines=FALSE,
                   lonR=NULL,latR=NULL,axiR=NULL,verbose=FALSE,plot=TRUE) {
  if (verbose) print('map.ds')
  #def.par <- par(no.readonly = TRUE) # save default, for resetting...
  stopifnot(inherits(x,'ds'))
  x <- subset(x,is=is)
  
  ## REB 2015-03-26
  if (inherits(x,'mvcomb')) {
    z <- map.mvcomb(x,it=it,verbose=verbose,new=new,
                    xlim=xlim,ylim=ylim,projection=projection,
                    lonR=lonR,latR=latR,axiR=axiR,gridlines=gridlines,
                    colbar=colbar,...) ##col=col,breaks=breaks)
    invisible(z)
  } else if (inherits(x,'pca')) {
    z <- map.pca(x,it=it,verbose=verbose,new=new,
                 xlim=xlim,ylim=ylim,projection=projection,
                 lonR=lonR,latR=latR,axiR=axiR,gridlines=gridlines,
                 colbar=colbar,...) ##col=col,breaks=breaks)
    invisible(z)
  } else if (inherits(x,'eof')) {
    z <- map.eof(x,it=it,verbose=verbose,new=new,
                 xlim=xlim,ylim=ylim,projection=projection,
                 lonR=lonR,latR=latR,axiR=axiR,gridlines=gridlines,
                 colbar=colbar,...) ##col=col,breaks=breaks)
    invisible(z)
  }
  projection <- tolower(projection)
  if (!is.null(attr(x,'pattern')))  { 
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
    attr(X,'variable') <- varid(attr(x,'eof'))
    attr(X,'unit') <- 'weight'
  } else X <- NULL
  
  unit <- attr(x,'unit')
  if ( (is.na(unit) | is.null(unit)) ) unit <- " "
  for (i in 1:length(unit)) {
    if ((unit[i]=='degree Celsius') | (unit[i]=='deg C') | (unit[i]=='degC'))
      unit[i] <- 'degree*C'
  }
  
  if (is.null(X)) { 
    attr(X,'unit') <- unit
    attr(X,'source') <- attr(x,'source')
  }
  if ((plot) & (!is.null(X))) {
    if (projection=="lonlat") {
      z <- lonlatprojection(x=X,colbar=colbar,verbose=verbose,xlim=xlim,ylim=ylim,
                            type='fill',gridlines=gridlines,new=new,...)
      if (is.list(attr(x,'pattern'))) {
        Xa <- attr(x,'pattern')
        nms <- names(Xa)
        col <- c('black','darkgreen','grey','yellow','magenta','cyan',
                        'brown','white','green')
                        
        for (i in (2:length(nms))) 
          contour(lon(Xa[[i]]),lat(Xa[[i]]),Xa[[i]],add=TRUE,col=col[i])
      } else if (sum(is.element(type,'contour'))>0)
        contour(lon(X),lat(X),X,add=TRUE,col="grey50")
    } else if (projection=="sphere") {
      z <- map2sphere(x=X,lonR=lonR,latR=latR,axiR=axiR,
                      xlim=xlim,ylim=ylim,type=type,
                      gridlines=gridlines,colbar=colbar,
                      new=new,verbose=verbose,...)
    } else if (projection=="np") {
      z <- map2sphere(X,lonR=lonR,latR=90,axiR=axiR,
                      xlim=xlim,ylim=ylim,type=type,
                      gridlines=gridlines,colbar=colbar,
                      new=new,verbose=verbose,...)
    } else if (projection=="sp") {
      z <- map2sphere(X,lonR=lonR,latR=-90,axiR=axiR,
                      xlim=xlim,ylim=ylim,type=type,
                      gridlines=gridlines,colbar=colbar,
                      new=new,verbose=verbose,...)
    }
  } else z <- X
  invisible(z)
}

#' @exportS3Method
#' @export
map.field <- function(x,...,FUN='mean',it=NULL,is=NULL,new=FALSE,
                      projection="lonlat",
                      xlim=NULL,ylim=NULL,zlim=NULL,n=15,
                      colbar= list(pal=NULL,rev=FALSE,n=10,breaks=NULL,
                                   pos=0.05,show=TRUE,type="p",cex=2,h=0.6,v=1),
                      type=c("fill","contour"),gridlines=FALSE,
                      lonR=NULL,latR=NULL,axiR=NULL,verbose=FALSE,
                      na.rm=TRUE,plot=TRUE,add=TRUE) {
  
  stopifnot(inherits(x,'field'))
  if (verbose) print('map.field')
  #def.par <- par(no.readonly = TRUE)
  
  x <- subset(x,it=it,is=is)
  #print(length(x)); print(attr(x,'dimensions')[1:2])
  projection <- tolower(projection)
  if (FUN=='trend') FUN <- 'trend.coef'
  if (!is.null(xlim)) {
    if (xlim[1] < 0) x <- g2dl(x,greenwich=FALSE)
  }
  #str(X)
  X <- coredata(x)
  ## Fudge, since these don't work with apply
  if (FUN=='firstyear') attr(x,FUN) <- firstyear(x)
  if (FUN=='lastyear') attr(x,FUN) <- lastyear(x)
  
  natts <- names(attributes(x))
  ## REB 2019-01-30
  if (sum(is.element(natts,FUN))) {
    if (verbose) print(paste('Use attribute',FUN))
    X <- attr(x,FUN)
  } else
    ## If one time slice, then map this time slice
    if (dim(X)[1]==1) {
      if (verbose) print('One point in time')
      X <- coredata(x[1,])
    } else if (is.null(X)) {
      if (verbose) print('Data is a vector')
      X <- coredata(X)
    } else if (inherits(X,"matrix")) {
      if (verbose) {print(paste('Data is a matrix. FUN=',FUN)); print(dim(x))}
      ## If several time slices, map the required statistics
      good <- apply(coredata(x),2,nv) > 1
      X <- rep(NA,length(good))
      xx <- x[,good]
      X[good] <- apply(coredata(xx),2,FUN=FUN,na.rm=na.rm)
    } else {print('Do not know what to do')}
  
  ## if zlim is specified, then mask data outside this range
  if (!is.null(zlim)) {
    d <- dim(X)
    mask <- (X < min(zlim)) | (X > max(zlim))
    rng <- range(X,na.rm=TRUE)
    X[mask] <- NA
    dim(X) <- d
    if (verbose) print(paste('zlim=',zlim[1],zlim[2],
                             '  sum(mask)=',sum(mask),
                             '  range(X)=',rng[1],rng[2]))
  }
  
  #print(length(X))
  attr(X,'longitude') <- attr(x,'longitude')
  attr(X,'latitude') <- attr(x,'latitude')
  attr(X,'variable') <- attr(x,'variable')[1]
  #  if (attr(x,'unit')=="deg C") attr(X,'unit') <- expression(degree*C) else
  unit <- attr(x,'unit')[1]
  
  if (unit =='%') unit <- "'%'"
  if ( (is.na(unit) | is.null(unit)) ) unit <- " "
  if ((unit=='degree Celsius') | (unit=='deg C') | (unit=='degC'))
    unit <- 'degree*C'
  
  unit <- as.character(unit)
  attr(X,'unit') <- unit
  attr(X,'source') <- attr(x,'source')
  attr(X,'time') <- range(index(x))
  attr(X,'method') <- FUN
  attr(X,'timescale') <- class(x)[2]
  if (verbose) {print(length(X)); print(attr(x,'dimensions'))}
  dim(X) <- attr(x,'dimensions')[1:2]
  
  if (verbose) {print(str(X)); print(summary(c(X)))}
  if (plot) {
    if (projection=="lonlat") {
      z <- lonlatprojection(x=X,xlim=xlim,ylim=ylim,zlim=zlim,n=n,
                            colbar=colbar,type=type,new=new,
                            gridlines=gridlines,verbose=verbose,...)
    } else if (projection=="sphere") {
      z <- map2sphere(x=X,xlim=xlim,ylim=ylim,zlim=zlim,n=n,
                      lonR=lonR,latR=latR,axiR=axiR,
                      type=type,gridlines=gridlines,
                      colbar=colbar,new=new,verbose=verbose,...)
    } else if (projection=="np") {
      z <- map2sphere(X,xlim=xlim,ylim=ylim,zlim=zlim,n=n,
                      lonR=lonR,latR=90,axiR=axiR,
                      type=type,gridlines=gridlines,
                      colbar=colbar,new=new,verbose=verbose,...)
    } else if (projection=="sp") {
      z <- map2sphere(X,xlim=xlim,ylim=ylim,zlim=zlim,n=n,
                      lonR=lonR,latR=-90,axiR=axiR,
                      type=type,gridlines=gridlines,
                      colbar=colbar,new=new,verbose=verbose,...)
    }
  } else z <- X
  invisible(z)
}

#' @exportS3Method
#' @export
map.corfield <- function(x,...,it=NULL,is=NULL,new=FALSE,projection="lonlat",
                         xlim=NULL,ylim=NULL,zlim=NULL,n=15,
                         colbar= list(pal=NULL,rev=FALSE,n=NULL,
                                      breaks=seq(-1,1,by=0.05),pos=0.05,show=TRUE,
                                      type="p",cex=2,h=0.6,v=1),
                         type=c("fill","contour"),gridlines=FALSE,
                         lonR=NULL,latR=NULL,axiR=NULL,verbose=FALSE,plot=TRUE) {
  
  if (verbose) print("map.corfield")
  stopifnot(inherits(x,'corfield'))
  def.par <- par(no.readonly = TRUE) # save default, for resetting...
  x <- subset(x,it=it,is=is,verbose=verbose)
  projection <- tolower(projection)
  dim(x) <- attr(x,'dimensions')[1:2]
  ## if (!is.null(colbar)) colbar$pal <- varid(x)[2] ## AM 08-07-2015 comment 
  attr(x,'variable') <- paste(varid(x),collapse='/')
  
  #if (length(attr(x,'unit'))>1) attr(x,'unit') <- paste(attr(x,'unit'),collapse='/')
  attr(x,'unit') <- attr(x,'unit')[1]
  
  ## if zlim is specified, then mask data outside this range
  if (!is.null(zlim)) {
    if (verbose) print(zlim)
    d <- dim(x)
    mask <- (x < min(zlim,na.rm=TRUE)) | (x > max(zlim,na.rm=TRUE))
    x[mask] <- NA
    dim(x) <- d
    if (verbose) {print(zlim); print(dim(x)); print(sum(mask))}
  }
  
  if (verbose) {print(projection); print(dim(x))}
  
  if (plot) {
    if (projection=="lonlat") {
      z <- lonlatprojection(x,xlim=xlim,ylim=ylim,zlim=zlim,n=n,
                            colbar=colbar,type=type,new=new,verbose=verbose,
                            gridlines=gridlines,...)
    } else if (projection=="sphere") {
      z <- map2sphere(x,xlim=xlim,ylim=ylim,zlim=zlim,n=n,
                      lonR=lonR,latR=latR,axiR=axiR,type=type,gridlines=gridlines,
                      colbar=colbar,new=new,verbose=verbose,...)
    } else if (projection=="np") {
      z <- map2sphere(x,xlim=xlim,ylim=ylim,zlim=zlim,n=n,
                      lonR=lonR,latR=90,axiR=axiR,type=type,gridlines=gridlines,
                      colbar=colbar,new=new,verbose=verbose,...) 
    } else if (projection=="sp") {
      z <- map2sphere(x,xlim=xlim,ylim=ylim,zlim=zlim,n=n,
                      lonR=lonR,latR=-90,axiR=axiR,type=type,gridlines=gridlines,
                      colbar=colbar,new=new,verbose=verbose,...)
    }
  } else z <- x
  if (!is.null(attr(x,'x.longitude')) & !is.null(attr(x,'x.latitude')))
    points(attr(x,'x.longitude'),attr(x,'x.latitude'),lwd=2,cex=1.2)
  invisible(z)
}

#' @exportS3Method
#' @export
map.trend <- function(x,...,it=NULL,is=NULL,new=FALSE,projection="lonlat",
                      xlim=NULL,ylim=NULL,zlim=NULL,n=15,
                      colbar= list(pal=NULL,rev=FALSE,n=10,breaks=NULL,
                                   pos=0.05,show=TRUE,type="p",cex=2,h=0.6,v=1),
                      type=c("fill","contour"),gridlines=FALSE,
                      lonR=NULL,latR=NULL,axiR=NULL,verbose=FALSE,plot=TRUE) {
  if (verbose) print('map.trend')
  stopifnot(inherits(x,'field'),inherits(x,'trend'))
  def.par <- par(no.readonly = TRUE) # save default, for resetting...
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
  if (plot) {
    if (projection=="lonlat") {
      z <- lonlatprojection(x=x,xlim=xlim,ylim=ylim,zlim=zlim,n=n,
                            colbar=colbar,type=type,new=new,
                            verbose=verbose,gridlines=gridlines,...)
    } else if (projection=="sphere") {
      z <- map2sphere(x=x,xlim=xlim,ylim=ylim,zlim=zlim,n=n,
                      lonR=lonR,latR=latR,axiR=axiR,
                      type=type,gridlines=gridlines,
                      colbar=colbar,new=new,verbose=verbose,...)
    } else if (projection=="np") {
      z <- map2sphere(x,xlim=xlim,ylim=ylim,zlim=zlim,n=n,
                      lonR=lonR,latR=90,axiR=axiR,
                      type=type,gridlines=gridlines,
                      colbar=colbar,new=new,verbose=verbose,...)
    } else if (projection=="sp") {
      z <- map2sphere(x,xlim=xlim,ylim=ylim,zlim=zlim,n=n,
                      lonR=lonR,latR=-90,axiR=axiR,
                      type=type,gridlines=gridlines,
                      colbar=colbar,new=new,verbose=verbose,...)
    }
  } else z <- X
  invisible(z)
}


#' @exportS3Method
#' @export map.pca
map.pca <- function(x,...,it=NULL,is=NULL,ip=1,new=FALSE,add=FALSE,projection="lonlat",
                    xlim=NULL,ylim=NULL,zlim=NULL,FUN='mean',##n=15,
                    colbar=list(pal=NULL,rev=FALSE,n=10,breaks=NULL,
                                pos=0.05,show=TRUE,type="p",cex=1,h=0.6,v=1),
                    #cex.axis=1,cex.main=1,cex.lab=1,
                    type=c("fill","contour"),gridlines=FALSE,fig=c(0,1,0.05,0.95),
                    lonR=NULL,latR=NULL,axiR=NULL,verbose=FALSE,plot=TRUE) {
  ##
  if (verbose) print(paste('map.pca',FUN))
  def.par <- par(no.readonly = TRUE) # save default, for resetting...
  if(inherits(x,"trajectory")) {
    z <- map.pca.trajectory(x,projection=projection,lonR=lonR,latR=latR,
                            xlim=xlim,ylim=ylim,...)
  } else {
    args <- list(...)
    #print(args)
    ## REB 2016-11-02 fix
    if (is.null(dim(attr(x,'pattern')))) {
      dim(attr(x,'pattern')) <- c(1,length(attr(x,'pattern')))
    }
    X <- rbind(attr(x,'pattern')[,ip],attr(x,'pattern')[,ip])
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
    attr(X,'mean') <- NULL
    class(X) <- 'station'
    ##if (is.null(colbar$col) | is.null(colbar)) {
    ##  colbar$col <- colscal(30,pal=varid(x))
    ##}
    if (verbose) str(X)
    
    if (is.element(FUN,args)) {
      z <- map.station(X,new=new,colbar=colbar,
                       xlim=xlim,ylim=ylim,zlim=zlim,
                       plot=TRUE,add=add,fig=fig,verbose=verbose,...) -> z
    } else {
      z <- map.station(X,new=new,colbar=colbar,FUN=FUN,
                       xlim=xlim,ylim=ylim,zlim=zlim,
                       plot=TRUE,add=add,fig=fig,verbose=verbose,...) -> z
    } 
  } 
  invisible(z)
}

#' @exportS3Method
#' @export
map.mvr <- function(x,...,it=NULL,is=NULL,new=FALSE,projection="lonlat",
                    xlim=NULL,ylim=NULL,zlim=NULL,
                    colbar= list(pal=NULL,rev=FALSE,n=10,breaks=NULL,
                                 pos=0.05,show=TRUE,type="p",cex=2,h=0.6,v=1),
                    type=c("fill","contour"),gridlines=FALSE,
                    verbose=FALSE,plot=TRUE) {
  if(verbose) print("map.mvr")
  def.par <- par(no.readonly = TRUE) # save default, for resetting...
  x <- subset(x,it=it,is=is)
  z <- map.field(x,new=new,FUN="mean",colbar=colbar,xlim=xlim,ylim=ylim,
                 verbose=verbose,plot=TRUE,...) -> z
  invisible(z)
  
}

#' @exportS3Method
#' @export
map.cca <- function(x,...,ip=1,it=NULL,is=NULL,new=FALSE,projection="lonlat",
                    xlim=NULL,ylim=NULL,zlim=NULL,
                    colbar1=list(pal=NULL,rev=FALSE,n=10,
		                 breaks=seq(-1,1,by=0.1),type="p",
		                 cex=2,show=FALSE,h=0.6,v=1,pos=0.05),
                    colbar2= NULL,
                    type=c("fill","contour"),gridlines=FALSE,
                    lonR=NULL,latR=NULL,axiR=NULL,cex=2,
		    plot=TRUE,verbose=FALSE) {
  if (verbose) print('map.cca')

  ## KMP 2023-02-07: It's not a good idea to set the layout here
  ## because map.cca is used in plot.cca and this will interfere with
  ## the layout specified there
  #layout(matrix(c(1,2),1,2,byrow = TRUE), rep(1,2), rep(1,2), TRUE)
  
  if (is.null(colbar2)) colbar2 <- colbar1
  ##x <- subset(x,it=it,is=is)
  
  ## For plotting, keep the same kind of object, but replace the patterns in
  ## the eof/pca with the CCA patterns
  Y <- x$Y
  ##print(dim(attr(Y,'pattern'))); print(dim(U))
  ##attr(Y,'pattern') <- U
  U <- x$B.m
  dim(U) <- c(dim(attr(Y,'pattern'))[-length(dim(attr(Y,'pattern')))],
              length(x$ip))
  attr(Y,'pattern') <- U
  attr(Y,'eigenvalues') <- rep(1,length(x$ip))
  attr(Y,'time') <- range(index(x))
  X <- x$X
  ##print(dim(attr(X,'pattern'))); print(dim(V))
  ##attr(X,'pattern') <- V
  V <- x$A.m
  dim(V) <- c(dim(attr(X,'pattern'))[-length(dim(attr(X,'pattern')))],
              length(x$ip))
  attr(X,'pattern') <- V
  attr(X,'eigenvalues') <- rep(1,length(x$ip))
  attr(X,'time') <- range(index(x))
  z1 <- map(Y,ip=ip,xlim=xlim,ylim=ylim,type=type,cex=cex,
            projection=projection,lonR=lonR,latR=latR,axiR=axiR,
            gridlines=gridlines,FUN='mean',verbose=verbose,
            colbar=colbar1,showall=FALSE,new=FALSE,plot=TRUE)

  z2 <- map(X,ip=ip,xlim=xlim,ylim=ylim,type=type,cex=cex,
            projection=projection,lonR=lonR,latR=latR,axiR=axiR,
            gridlines=gridlines,FUN='mean',verbose=verbose,
            colbar=colbar2,showall=FALSE,new=FALSE,plot=TRUE)
  
  invisible(list(z1 = z1, z2 = z2))
}

#' @exportS3Method
#' @export
map.events <- function(x,Y=NULL,...,it=NULL,is=NULL,xlim=NULL,ylim=NULL,main=NULL,
                       param=NA,alpha=0.3,lwd=3,col="black",bg="white",pch=21,cex=1,
                       colbar=list(pal="budrd",rev=FALSE,n=10,breaks=NULL,
                                   pos=0.05,show=TRUE,type="p",cex=2,h=0.6,v=1),
                       showaxis=TRUE,fig=c(0,1,0.05,0.95),mgp=c(2,0.5,0),mar=rep(2,4),
                       lty=1,type=c("points","trajectory","start","end"),
                       border=FALSE,
                       projection="lonlat",latR=NULL,lonR=NULL,new=TRUE,
                       verbose=FALSE) {
  if(verbose) print("map.events")
  x0 <- x
  x <- subset(x,it=it,is=is,verbose=verbose)
  if(is.null(attr(x,"calendar"))) calendar <- "gregorian" else calendar <- attr(x,"calendar")
  if(is.null(it) & dim(x)[1]>0) {
    if (requireNamespace("PCICt", quietly = TRUE)) {
      it <- range(PCICt::as.PCICt(as.character(x$date),cal=calendar,format="%Y%m%d"))
    } else {
      it <- range(as.Date(as.character(x$date),cal=calendar,format="%Y%m%d"))
    }
  }
  
  if (is.null(is$lon) & !is.null(xlim)) {
    is$lon <- xlim
  } else if (is.null(is$lon) & is.null(xlim)) {
    if(length(Y)>0) {
      is$lon <- range(lon(Y))
    } else if(dim(x)[1]>0) {
      is$lon <- range(x[,"lon"])+c(-5,5)
    }
  }
  if (is.null(xlim) & projection=="lonlat") xlim <- is$lon
  
  if(projection=="lonlat" & !any(xlim<0) & any(xlim>180)) x <- g2dl(x,greenwich=TRUE)
  
  if (is.null(is$lat) & !is.null(ylim)) {
    is$lat <- ylim
  } else if (is.null(is$lat) & is.null(ylim)) {
    if(length(Y)>0) {
      is$lat <- range(lat(Y))
    } else if(dim(x)[1]>0) {
      is$lat <- range(x[,"lat"])+c(-2,2)
    }
  }
  if (is.null(ylim) & projection=="lonlat") ylim <- is$lat
  
  if (!is.null(Y)) {
    Y <- subset(Y,is=is)
    if(length(Y)>0) {
      if(dim(x)[1]==0) {
        Y <- subset(Y,it=it)
      } else {
        ty <- index(Y)
        if (inherits(Y,"month")) {
          tx <- round(x[,"date"]*1E-2)*1E2+1
          ty <- as.numeric(format(ty,"%Y%m%d"))
        } else if (inherits(ty,"Date")) {
          tx <- x[,"date"]
          ty <- as.numeric(format(ty,"%Y%m%d"))
        } else if (inherits(ty,c("POSIXt","PCICt"))) {
          tx <- x[,"date"]*1E2 + x[,"time"]
          ty <- as.numeric(format(ty,"%Y%m%d%H"))
        }
        ii <- is.element(ty,tx)
        Y <- subset(Y,it=ii)
      }
    }
  }
  if(length(Y)!=0) {
    if (is.null(lonR)) lonR <- mean(lon(Y))
    if (is.null(latR)) latR <- max(lat(Y))
    z <- map(Y,colbar=colbar,new=new,projection=projection,main="",
             fig=fig,mar=mar,mgp=mgp,showaxis=showaxis,
             add=add,xlim=xlim,ylim=ylim,latR=latR,lonR=lonR,verbose=verbose)
  } else {
    if(!is.null(xlim)) {
      lonR <- mean(xlim)
    } else if (is.null(lonR)) {
      if (dim(x)[1]>0) {
        lonR <- mean(x[,"lon"])
      } else {
        lonR <- 0
      }
    }
    if(!is.null(ylim)) {
      latR <- mean(ylim)
      #latR <- sign(ylim[ylim==max(abs(ylim))])*max(abs(ylim))
    } else if (is.null(latR)) {
      if(dim(x)[1]>0) {
        latR <- mean(x[,"lat"])
        #latR <- sign(x[,"lat"][x[,"lat"]==max(abs(x[,"lat"]))])*max(abs(x[,"lat"]))
      } else {
        latR <- 90
      }
    }
    xs <- events2station(x, FUN="location", param="pcent", verbose=verbose)
    map(xs, FUN="mean", col="grey", cex=0.1, pch='.',
        new=new,projection=projection,main="",xlab="",ylab="",
        fig=fig,mar=mar,mgp=mgp,showaxis=showaxis,
        border=border,
        xlim=xlim,ylim=ylim,latR=latR,lonR=lonR,
        verbose=verbose)
  }
  if(dim(x)[1]>0) {
    cols <- adjustcolor(col,alpha.f=alpha)
    if("points" %in% type) {
      if(verbose) print("plot points")
      if(projection=="lonlat") {
        points(x[,"lon"],x[,"lat"],col=cols,bg=bg,cex=cex,pch=pch,lwd=lwd)
      } else {
        theta <- pi*x[,"lon"]/180
        phi <- pi*x[,"lat"]/180
        ax <- sin(theta)*cos(phi)
        ay <- cos(theta)*cos(phi)
        az <- sin(phi)
        if (verbose) {print('Rotation:');print(dim(rotM(x=0,y=0,z=lonR))); print(dim(rbind(ax,ay,az)))}
        a <- rotM(x=0,y=0,z=lonR) %*% rbind(ax,ay,az)
        a <- rotM(x=latR,y=0,z=0) %*% a
        ax <- a[1,]; ay <- a[2,]; az <- a[3,]
        points(ax[ay>0],az[ay>0],col=cols,bg=bg,cex=cex,pch=pch,lwd=lwd)    
      }
    }
    if("trajectory" %in% colnames(x0) &
       any(c("trajectory","start","end") %in% type)) {
      xt <- subset(x0,it=(x0$trajectory %in% x$trajectory))
      if(!("trackcount" %in% names(x)) & dim(xt)[1]>1) {
        xt <- trackstats(xt)
        xt <- subset(xt,it=xt$trackcount>1)
      }
      if(dim(xt)[1]>1) {
        xall <- as.trajectory(xt,nmin=2,n=45,verbose=verbose)
        map(xall,lty=lty,lwd=lwd,alpha=alpha,new=FALSE,
            col=col,lonR=lonR,latR=latR,
            projection=projection,type=type,param=param,
            showaxis=FALSE,
            colbar=colbar,verbose=verbose,...)
      }
    }
  }
  period <- unique(c(min(it),max(it)))
  if (!is.null(period) & length(Y)==0) {
    title(sub = paste(period,collapse=" - "))
    # text(par("usr")[1] + 0.05*diff(range(par("usr")[3:4])),
    #      par("usr")[4] - 0.05*diff(range(par("usr")[3:4])),
    #      paste(period,collapse=" - "),pos=4,cex=0.75,col="grey30")
  }
  
  if (!is.null(main)) {
    text(par("usr")[1] + 0.05*diff(range(par("usr")[3:4])),
         par("usr")[4] - 0.10*diff(range(par("usr")[3:4])),
         main,pos=4,cex=1,col="black")
  }
  #par(def.par) # reset to default
}

#' Function that masks either ocean or land
#'
#' Uses topography from \code{\link{etopo5}} to mask either land or ocean.
#'
#' @param x a \code{field} object 
#' @param land a boolean; if TRUE mask land, if FALSE mask ocean
#'
#' @export
mask <- function(x,land=FALSE) {
  data(etopo5, envir = environment())
  h <- regrid(etopo5,is=x)
  if (!land) {
    h[h < -5] <- NA
  } else {
    h[h > 5] <- NA
  }
  X <- coredata(x)
  X[,is.na(h)] <- NA
  X -> coredata(x)
  return(x)
}
