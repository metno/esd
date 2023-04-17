# Plot esd objects
# 
# These plot functions are S3 methods for esd objects, based on \code{plot}.
# 
# @importFrom graphics par grid segments text axis legend polygon mtext abline layout rect boxplot
# @importFrom grDevices dev.new dev.off dev.copy2eps dev.copy2pdf dev.list
# @importFrom stats pbinom
# 
# @param x the object to be plotted
# @param plot.type "single"
# @param new if TRUE plot in new window
# @param lwd width of line
# @param type type of plot: 'l' = line, 'p' = point, 'b' = both
# @param pch type of marker
# @param main main title
# @param xlab label of x-axis
# @param ylab label of y-axis
# @param errorbar if TRUE show errorbar
# @param legend.show if TRUE show legend
# @param map.show show map of stations
# @param map.type 'points' to show stations on map, 'rectangle' to show area
# @param map.insert if TRUE show map as insert, else show map in new window
# @param it For subsetting in time - See \code{\link{subset}}.
# @param is For subsetting in space - See \code{\link{subset}}. Can also be
# a station value and if provided, the plotting will involve an interpolation
# to the same coordinates as defined by \code{is}.
# @param ip Which EOF/CCA pattern (mode) to plot
# @param cex magnification factor, see \code{\link[graphics]{par}}
# @param cex.axis see \code{\link[graphics]{par}}
# @param cex.lab see \code{\link[graphics]{par}}
# @param cex.main see \code{\link[graphics]{par}}
# @param mar see \code{\link[graphics]{par}}
# @param fig coordinates of figure region, see \code{\link[graphics]{par}}
# @param alpha transparency factor for main plot
# @param alpha.map transparency factor for map
# @param verbose a boolean; if TRUE print information about progress
# @param col Colour see \code{\link[graphics]{par}}
# @param lwd width of line
# @param xlim range of x-axis
# @param ylim range of y-axis
# @param what Indicate what to plot. For \code{plot.eof}, c('pc', 'eof', 'var') is the default setting which means that
# the plot will include the principle components, EOF patterns and explained variance. 'field' expands eof to field before
# plotting
# @param colbar a list, see \code{\link{colbar}}
# @param \dots additional arguments
# 
# @return None
# 
# @keywords plot graphics
# 
# @examples
# 
# # Example: use aggregate to compute annual mean temperature for Svalbard:
# data(Svalbard)
# y <- aggregate(Svalbard, by=year(Svalbard), FUN='mean', na.rm = FALSE)
# plot(y, new=FALSE)
# 
# # Example with downscaling:
# lon <- c(-12,37)
# lat <- c(52,72)
# t2m <- t2m.DNMI(lon=lon,lat=lat)
# data(Oslo)
# ds <- DS(Oslo,t2m)
# 
# # Plot the results for January month
# # plot(subset(ds,it='Jan'))
# 
# # Plot the residuals:
# residual <- as.residual(ds)
# obs <- as.anomaly(as.calibrationdata(ds))
# 
# plot.zoo(obs[,1], lwd=2, new=FALSE)
# lines(residual, col="red")
# 
# print("Global climate model simulation NorESM")
# T2m <- t2m.NorESM.M(lon=lon,lat=lat)
# 
# # Plot the global mean of the field:
# plot(T2m, new=FALSE)
# # Plot area mean of a sub region
# plot(T2m,is=list(lon=c(0,10),lat=c(60,70)), new=FALSE)
# 
# # Plot interpolated results corresponding to ferder
# data(ferder)
# plot(T2m,ferder, new=FALSE)
# 
# # Plot Hovmuller diagram: Not working ...
# ## plot(T2m,is=list(lon=0))
# 
# print("Extract a subset - the January month")
# x <- subset(t2m,it="jan")
# X <- subset(T2m,it="jan")
# 
# print("Combine the fields for computing common EOFs:")
# XX <- combine(x,X)
# 
# print("Compute common EOFs")
# eofxx <- EOF(XX)
# plot(eofxx, new=FALSE)
# 
# print("Downscale the January mean temperature")
# ds.jan <- DS(Oslo,eofxx)
# plot(ds.jan, new=FALSE)
# 
# @export 
# plot <- function(x,...)  UseMethod("plot")

#' Plot esd objects
#' 
#' These plot functions are S3 methods for esd objects, based on \code{plot}.
#'
#' @importFrom graphics par grid segments text axis legend polygon mtext abline layout rect boxplot
#' @importFrom grDevices dev.new dev.off dev.copy2eps dev.copy2pdf dev.list
#' @importFrom stats pbinom
#'
#' @seealso plot.list plot.numeric plot.xsection plot.eof plot.station plot.ds plot.eof plot.pca plot.spell plot.nevents plot.IDF
#'
#' @param x the object to be plotted
#' @param plot.type "single"
#' @param new if TRUE plot in new window
#' @param lwd width of line
#' @param type type of plot: 'l' = line, 'p' = point, 'b' = both
#' @param pch type of marker
#' @param main main title
#' @param xlab label of x-axis
#' @param ylab label of y-axis
#' @param errorbar if TRUE show errorbar
#' @param legend.show if TRUE show legendp
#' @param map.show show map of stations
#' @param map.type 'points' to show stations on map, 'rectangle' to show area
#' @param map.insert if TRUE show map as insert, else show map in new window
#' @param it For subsetting in time - See \code{\link{subset}}.
#' @param is For subsetting in space - See \code{\link{subset}}. Can also be
#' a station value and if provided, the plotting will involve an interpolation
#' to the same coordinates as defined by \code{is}.
#' @param ip Which EOF/CCA pattern (mode) to plot
#' @param cex magnification factor, see \code{\link[graphics]{par}}
#' @param cex.axis see \code{\link[graphics]{par}}
#' @param cex.lab see \code{\link[graphics]{par}}
#' @param cex.main see \code{\link[graphics]{par}}
#' @param mar see \code{\link[graphics]{par}}
#' @param fig coordinates of figure region, see \code{\link[graphics]{par}}
#' @param alpha transparency factor for main plot
#' @param alpha.map transparency factor for map
#' @param verbose a boolean; if TRUE print information about progress
#' @param col Colour see \code{\link[graphics]{par}}
#' @param lwd width of line
#' @param xlim range of x-axis
#' @param ylim range of y-axis
#' @param what Indicate what to plot. For \code{plot.eof}, c('pc', 'eof', 'var') is the default setting which means that
#' the plot will include the principle components, EOF patterns and explained variance. 'field' expands eof to field before
#' plotting
#' @param colbar a list, see \code{\link{colbar}} 
#' @param \dots additional arguments
#' 
#' @return None
#'
#' @keywords plot graphics
#'
#' @examples
#' 
#' # Example: use aggregate to compute annual mean temperature for Svalbard:
#' data(Svalbard)
#' y <- aggregate(Svalbard, by=year(Svalbard), FUN='mean', na.rm = FALSE) 
#' plot(y, new=FALSE)
#' 
#' # Example with downscaling:
#' lon <- c(-12,37)
#' lat <- c(52,72)
#' t2m <- t2m.DNMI(lon=lon,lat=lat)
#' data(Oslo)
#' ds <- DS(Oslo,t2m)
#' 
#' # Plot the results for January month
#' plot(subset(ds,it='Jan'), new=FALSE)
#' 
#' # Plot the residuals:
#' residual <- as.residual(ds)
#' obs <- as.anomaly(as.calibrationdata(ds)$y)
#' plot.zoo(obs[,1], lwd=2, new=FALSE)
#' lines(residual, col="red")
#' 
#' print("Global climate model simulation NorESM")
#' T2m <- t2m.NorESM.M(lon=lon,lat=lat)
#' 
#' # Plot the global mean of the field:
#' plot(T2m, new=FALSE)
#' # Plot area mean of a sub region
#' plot(T2m, is=list(lon=c(0,10),lat=c(60,70)), new=FALSE)
#' 
#' # Plot interpolated results corresponding to ferder
#' data(ferder)
#' plot(T2m, ferder, new=FALSE)
#' 
#' # Plot Hovmuller diagram: Not working ...
#' ## plot(T2m,is=list(lon=0)) 
#' 
#' print("Extract a subset - the January month")
#' x <- subset(t2m,it="jan")
#' X <- subset(T2m,it="jan")
#' 
#' print("Combine the fields for computing common EOFs:")
#' XX <- combine(x,X)
#' 
#' print("Compute common EOFs")
#' eofxx <- EOF(XX)
#' plot(eofxx, new=FALSE)
#' 
#' print("Downscale the January mean temperature") 
#' ds.jan <- DS(Oslo,eofxx)
#' plot(ds.jan, new=FALSE)
#'
#' @export plot.double
plot.double <- plot.default

#' @exportS3Method
#' @export
plot.numeric <- plot.default

#' @exportS3Method
#' @export
plot.list <- function(x,...,is=NULL,
                      col=c(rgb(1,1,0.5,0.05),rgb(1,0.5,0.5,0.05),rgb(0.5,1,0.5,0.05)),
                      lwd=3,xlim=NULL,ylim=NULL) {
  if (!is.null(is)) y <- subset(x,it=is) else y <- x[[1]]
  plot(y,img=img,col=col[1],lwd=lwd,xlim=xlim,ylim=ylim)
  for (j in c(2:length(x),1)) {
    if (!is.null(it)) y <- subset(x[[j]],it=it) else y <- x[[j]]
    for (i in 1:dim(y)[2]) lines(y[,i],lwd=7,col=col[j])
    lines(attr(y,'station'),lwd=3,col=rgb(0.5,0.5,0.5,0.25))
  }
}

#' Plot esd objects
#' 
#' These plot functions are S3 methods for esd objects, based on \code{plot}.
#'
#' @importFrom graphics par grid segments text axis legend polygon mtext abline layout rect boxplot
#' @importFrom grDevices dev.new dev.off dev.copy2eps dev.copy2pdf dev.list
#' @importFrom stats pbinom
#'
#' @param x the object to be plotted
#' @param plot.type "single"
#' @param new if TRUE plot in new window
#' @param lwd width of line
#' @param type type of plot: 'l' = line, 'p' = point, 'b' = both
#' @param pch type of marker
#' @param main main title
#' @param xlab label of x-axis
#' @param ylab label of y-axis
#' @param errorbar if TRUE show errorbar
#' @param legend.show if TRUE show legendp
#' @param map.show show map of stations
#' @param map.type 'points' to show stations on map, 'rectangle' to show area
#' @param map.insert if TRUE show map as insert, else show map in new window
#' @param it For subsetting in time - See \code{\link{subset}}.
#' @param is For subsetting in space - See \code{\link{subset}}. Can also be
#' a station value and if provided, the plotting will involve an interpolation
#' to the same coordinates as defined by \code{is}.
#' @param ip Which EOF/CCA pattern (mode) to plot
#' @param cex magnification factor, see \code{\link[graphics]{par}}
#' @param cex.axis see \code{\link[graphics]{par}}
#' @param cex.lab see \code{\link[graphics]{par}}
#' @param cex.main see \code{\link[graphics]{par}}
#' @param mar see \code{\link[graphics]{par}}
#' @param fig coordinates of figure region, see \code{\link[graphics]{par}}
#' @param alpha transparency factor for main plot
#' @param alpha.map transparency factor for map
#' @param verbose a boolean; if TRUE print information about progress
#' @param col Colour see \code{\link[graphics]{par}}
#' @param lwd width of line
#' @param xlim range of x-axis
#' @param ylim range of y-axis
#' @param what Indicate what to plot. For \code{plot.eof}, c('pc', 'eof', 'var') is the default setting which means that
#' the plot will include the principle components, EOF patterns and explained variance. 'field' expands eof to field before
#' plotting
#' @param colbar a list, see \code{\link{colbar}} 
#' @param \dots additional arguments
#' 
#' @return None
#'
#' @keywords plot graphics
#'
#' @examples
#' 
#' # Example: use aggregate to compute annual mean temperature for Svalbard:
#' data(Svalbard)
#' y <- aggregate(Svalbard, by=year(Svalbard), FUN='mean', na.rm = FALSE) 
#' plot(y, new=FALSE)
#' 
#' # Example with downscaling:
#' lon <- c(-12,37)
#' lat <- c(52,72)
#' t2m <- t2m.DNMI(lon=lon,lat=lat)
#' data(Oslo)
#' ds <- DS(Oslo,t2m)
#' 
#' # Plot the results for January month
#' # plot(subset(ds,it='Jan'))
#' 
#' # Plot the residuals:
#' residual <- as.residual(ds)
#' obs <- as.anomaly(as.calibrationdata(ds))
#' 
#' plot.zoo(obs[,1], lwd=2, new=FALSE)
#' lines(residual, col="red")
#' 
#' print("Global climate model simulation NorESM")
#' T2m <- t2m.NorESM.M(lon=lon,lat=lat)
#' 
#' # Plot the global mean of the field:
#' plot(T2m, new=FALSE)
#' # Plot area mean of a sub region
#' plot(T2m,is=list(lon=c(0,10),lat=c(60,70)), new=FALSE)
#' 
#' # Plot interpolated results corresponding to ferder
#' data(ferder)
#' plot(T2m,ferder, new=FALSE)
#' 
#' print("Extract a subset - the January month")
#' x <- subset(t2m,it="jan")
#' X <- subset(T2m,it="jan")
#' 
#' print("Combine the fields for computing common EOFs:")
#' XX <- combine(x,X)
#' 
#' print("Compute common EOFs")
#' eofxx <- EOF(XX)
#' plot(eofxx, new=FALSE)
#' 
#' print("Downscale the January mean temperature") 
#' ds.jan <- DS(Oslo,eofxx)
#' plot(ds.jan, new=FALSE)
#'
#' @exportS3Method
#' @export plot.station
plot.station <- function(x,...,plot.type="single",new=TRUE,
                         lwd=3,type='l',pch=0,main=NULL,col=NULL,
                         xlim=NULL,ylim=NULL,xlab="",ylab=NULL,
                         errorbar=TRUE,legend.show=FALSE,
                         map.show=TRUE,map.type=NULL,map.insert=TRUE,
                         xrange=NULL,yrange=NULL,
                         cex.axis=1.2,cex.lab=1.2,cex.main=1.2,
                         mar=NULL,
                         alpha=0.5,alpha.map=0.7,#add=FALSE,
                         verbose=FALSE) {
  
  if (verbose) print('plot.station')
  par0 <- par()
  par(las=1)
  
  if (!is.numeric(lon(x)) | !is.numeric(lat(x))) {
    map.show <- FALSE
  }
  if(map.show) {
    if (verbose) print('show map')
    if (is.null(map.type)) {
      if( inherits(x,"field") | length(lon(x))!=length(lat(x)) |
          (length(lon(x))==2 & length(lat(x))==2) ) {
        map.type <- "rectangle"
      } else {
        map.type <- "points"
      }
    }
    if (verbose) print(map.type)
  }
  
  if(is.null(mar)) {
    if(map.show & map.insert) {
      mar <- c(4.5, 4.5, 0.75, 0.5)
    } else {
      mar <- c(4.5, 4.5, 4.5, 1.0)
    }
  }

  if (is.null(ylim)) {
    ylim <- range(pretty(range(x,na.rm=TRUE)))
    ylim <- ylim + diff(ylim)*0.1*c(-0.05,0.05)
    if(legend.show) ylim[1] <- ylim[1]-diff(ylim)*0.05
    if(map.show & map.insert) ylim[2] <- ylim[2]+diff(ylim)*0.05
  }
  if (is.null(xlim)) xlim <- range(index(x))
  if (verbose) {print(xlim); print(ylim)}
  
  if (plot.type=="single") {
    if (is.null(ylab)) {
      ylab <- esd::ylab(x) # ggplot2 ylab can interfere with esd
    }
    if (inherits(ylab,"try-error")) ylab <- attr(x,'unit')
  } else {
    if (is.null(ylab)) { 
      if ((length(levels(factor(stid(x))))>1) & (length(levels(factor(varid(x))))<=1)) {
        ylab <- stid(x)
      } else 
        ylab <- varid(x)
    } else {
      if (is.null(main)) {
        if ((length(levels(factor(stid(x))))>1) & (length(levels(factor(varid(x))))<=1)) {
          main <- levels(factor((attr(x,'longname'))))[1]
        } else {
          main <- levels(factor(loc(x)))[1]
        }
      }
    }  
  }
  #if (is.null(main)) main <- attr(x,'longname')[1] 
  if (is.null(col)) {
    if (is.null(dim(x))) {
      col <- "blue"
    } else if (!is.null(lon(x)) & !is.null(lat(x)) &
               length(lon(x))==dim(x)[2] &
               length(lat(x))==dim(x)[2]) {
      nx <- (lon(x)-min(lon(x)))/diff(range(lon(x)))
      ny <- (lat(x)-min(lat(x)))/diff(range(lat(x)))
      if ( all(is.finite(nx) & is.finite(ny)) ) {
        col <- rgb(1-ny,nx,ny,1)
      } else {
        col <- rainbow(dim(x)[2])
      }
    } else {
      col <- rainbow(length(x[1,]))  
    }
  }
  if(is.null(alpha.map)) alpha.map <- alpha
  col.map <- adjustcolor(col,alpha.f=alpha.map)
  col <- adjustcolor(col,alpha.f=alpha)
  
  ns <- length(stid(x))

  errorbar <- errorbar & !is.null(err(x))
  if(new) dev.new()
  
  if(map.show & !map.insert) {
    vis.map(x,col=col.map,map.type,add.text=FALSE,map.insert=FALSE,
            cex.axis=cex.axis,xrange=xrange,yrange=yrange,
            cex=1.8,verbose=verbose)
  }
  
  cls <- class(x)
  if("seasonalcycle" %in% cls) xaxt <- "n" else  xaxt <- NULL
  class(x) <- "zoo"
  
  par(cex.axis=1,mar=mar, bty="n",xaxt="s",yaxt="s",xpd=FALSE)
  plot.zoo(x,...,plot.type=plot.type,xlab=xlab,ylab=ylab,
           col=col,xlim=xlim,ylim=ylim,lwd=lwd,type=type,pch=pch,
           cex.axis=cex.axis,cex.lab=cex.lab,cex.main=cex.main,
           xaxt=xaxt,main=main)
  setfig <- FALSE
  fig <- par()$fig; usr <- par()$usr; xaxp <- par()$xaxp; 
  yaxp <- par()$yaxp; plt <- par()$plt
  par1 <- par()
  if("seasonalcycle" %in% cls) {
    axis(1,at=seq(1,12),labels=month.abb,cex.axis=cex.axis,las=2)
  }
  
  if (plot.type=="single") {
    if (errorbar) {
      segments(index(x),x-err(x),index(x),x+err(x),
               lwd=3,col=rgb(0.5,0.5,0.5,0.25))
    }
    
    if(legend.show) {
      legend("bottomleft",legend=paste(attr(x,'location'),": ",
                                 round(attr(x,'longitude'),2),"E/",
                                 round(attr(x,'latitude'),2),"N (",
                                 attr(x,'altitude')," masl)",sep=""),
             bty="n",cex=0.65,ncol=2,
             text.col="grey40",lty=1,col=col)
    }
    if (map.show & map.insert) {
      vis.map(x,col=col.map,map.type=map.type,cex=1,cex.axis=0.65,
              add.text=FALSE,map.insert=map.insert,
              xrange=xrange,yrange=yrange,verbose=verbose)
      setfig <- TRUE
    }
  }
  ## Don't reset fig unless it has been set within the function. This will ruin the use of layout and mfrow.
  if(setfig) par(fig=fig) 
  
  #dontset <- c("cin","cra","csi","cxy","din","page","fig") 
  #for(p in names(par1)[!names(par1) %in% dontset]) eval(parse(text=paste0("par(",p,"=par1$",p,")")))
  par(cex.axis=par0$cex.axis, mar=par0$mar, bty=par0$bty,
      xaxt=par0$xaxt, yaxt=par0$yaxt, xpd=par0$xpd)
}

#' Plot esd objects
#' 
#' These plot functions are S3 methods for esd objects, based on \code{plot}.
#'
#' @importFrom graphics par grid segments text axis legend polygon mtext abline layout rect boxplot
#' @importFrom grDevices dev.new dev.off dev.copy2eps dev.copy2pdf dev.list
#' @importFrom stats pbinom
#'
#' @aliases plot.eof.comb plot.eof.field plot.eof.var
#'
#' @param x the object to be plotted
#' @param plot.type "single"
#' @param new if TRUE plot in new window
#' @param lwd width of line
#' @param type type of plot: 'l' = line, 'p' = point, 'b' = both
#' @param pch type of marker
#' @param main main title
#' @param xlab label of x-axis
#' @param ylab label of y-axis
#' @param errorbar if TRUE show errorbar
#' @param legend.show if TRUE show legendp
#' @param map.show show map of stations
#' @param map.type 'points' to show stations on map, 'rectangle' to show area
#' @param map.insert if TRUE show map as insert, else show map in new window
#' @param it For subsetting in time - See \code{\link{subset}}.
#' @param is For subsetting in space - See \code{\link{subset}}. Can also be
#' a station value and if provided, the plotting will involve an interpolation
#' to the same coordinates as defined by \code{is}.
#' @param ip Which EOF/CCA pattern (mode) to plot
#' @param cex magnification factor, see \code{\link[graphics]{par}}
#' @param cex.axis see \code{\link[graphics]{par}}
#' @param cex.lab see \code{\link[graphics]{par}}
#' @param cex.main see \code{\link[graphics]{par}}
#' @param mar see \code{\link[graphics]{par}}
#' @param fig coordinates of figure region, see \code{\link[graphics]{par}}
#' @param add if TRUE add plot to existing figure
#' @param alpha transparency factor for main plot
#' @param alpha.map transparency factor for map
#' @param verbose a boolean; if TRUE print information about progress
#' @param col Colour see \code{\link[graphics]{par}}
#' @param lwd width of line
#' @param xlim range of x-axis
#' @param ylim range of y-axis
#' @param what Indicate what to plot. For \code{plot.eof}, c('pc', 'eof', 'var') is the default setting which means that
#' the plot will include the principle components, EOF patterns and explained variance. 'field' expands eof to field before
#' plotting
#' @param colbar a list, see \code{\link{colbar}} 
#' @param \dots additional arguments
#' 
#' @return None
#'
#' @keywords plot graphics
#'
#' @examples
#' 
#' # Example: use aggregate to compute annual mean temperature for Svalbard:
#' data(Svalbard)
#' y <- aggregate(Svalbard, by=year(Svalbard), FUN='mean', na.rm = FALSE)
#' plot(y, new=FALSE)
#' 
#' # Example with downscaling:
#' lon <- c(-12,37)
#' lat <- c(52,72)
#' t2m <- t2m.DNMI(lon=lon,lat=lat)
#' data(Oslo)
#' ds <- DS(Oslo,t2m)
#' 
#' # Plot the results for January month
#' plot(subset(ds,it='Jan'), new=FALSE)
#' 
#' # Plot the residuals:
#' residual <- as.residual(ds)
#' obs <- as.anomaly(as.calibrationdata(ds))
#' 
#' plot.zoo(obs[,1], lwd=2, new=FALSE)
#' lines(residual, col="red")
#' 
#' print("Global climate model simulation NorESM")
#' T2m <- t2m.NorESM.M(lon=lon,lat=lat)
#' 
#' # Plot the global mean of the field:
#' plot(T2m, new=FALSE)
#' # Plot area mean of a sub region
#' plot(T2m,is=list(lon=c(0,10),lat=c(60,70)), new=FALSE)
#' 
#' # Plot interpolated results corresponding to ferder
#' data(ferder)
#' plot(T2m,ferder, new=FALSE)
#' 
#' # Plot Hovmuller diagram: Not working ...
#' ## plot(T2m,is=list(lon=0))
#' 
#' print("Extract a subset - the January month")
#' x <- subset(t2m,it="jan")
#' X <- subset(T2m,it="jan")
#' 
#' print("Combine the fields for computing common EOFs:")
#' XX <- combine(x,X)
#' 
#' print("Compute common EOFs")
#' eofxx <- EOF(XX)
#' plot(eofxx, new=FALSE)
#' 
#' print("Downscale the January mean temperature") 
#' ds.jan <- DS(Oslo,eofxx)
#' plot(ds.jan, new=FALSE)
#'
#' @exportS3Method
#' @export plot.eof
plot.eof <- function(x,...,new=FALSE,xlim=NULL,ylim=NULL,
                     ip=1,what=c("pc","eof","var"),
                     colbar=list(pal=NULL,rev=FALSE,n=10,alpha=0.8,
                                 breaks=NULL,type="p",cex=2,show=TRUE,
                                 h=0.6,v=1,pos=0.05),
                     verbose=FALSE,is=NULL,it=NULL) {
  if (verbose) print(paste('plot.eof',paste(what,collapse=',')))
  if (inherits(x,"comb"))
    plot.eof.comb(x,new=new,xlim=xlim,ylim=ylim,
                  ip=ip,what=what,colbar=colbar,verbose=verbose,...) else
                    if (inherits(x,c("field","station")))
                      plot.eof.field(x,new=new,xlim=xlim,ylim=ylim,
                                     ip=ip,what=what,colbar=colbar,verbose=verbose,
                                     it=it,is=is,...) else
                                       print("x does not have 'comb' or 'field' aspects...")
}

#' @exportS3Method
#' @export plot.eof.field
plot.eof.field <- function(x,...,new=FALSE,xlim=NULL,ylim=NULL,ip=1,
                           what=c("pc","eof","var"), colbar=NULL,
                           cex.axis=0.9,cex.main=0.9,cex.lab=0.9,
                           verbose=FALSE,it=NULL,is=NULL,cex=1) {
  ##layout(matrix(c(1,2,3,3),nrow = 2,ncol = 2,byrow = TRUE)) # REB: this does not work well at the moment
  if (verbose) print(paste('plot.eof.field',paste(what,collapse=',')))
  ## Save the original graphics settings
  par0 <- par()
  n <- ip
  what <- tolower(what)
  if ('field' %in% what) {
    ## Expand EOF to original field before plotting
    if (verbose) print('Transform eof to field before plot')
    x <- subset(x,it=it,is=is)
    y <- as.field(x)
    z <- plot(y,xlim=xlim,ylim=ylim,new=new,...)
    invisible(z)
  }
  #str(ip); stop("HERE")
  D <- attr(x,'eigenvalues')
  tot.var <- attr(x,'tot.var')
  var.eof <- 100* D^2/tot.var
  if (length(what)==3) {
    mfrow <- c(2,2)
  } else
    if (length(what)==2) mfrow <- c(2,1) else
      if (length(what)==1) mfrow <- c(1,1)
  if (new) dev.new()
  ## par(cex.axis=0.75,cex.lab=0.7,cex.main=0.8)
  #par(mfrow=mfrow)##,mar=c(1,1,1,2)) ##,bty="n",xaxt="n",yaxt="n")
  if (length(grep('eof',what))>0) {
    if (verbose) {print('Show map'); print(class(x))}
    if(length(dim(attr(x,"pattern")))==3) {
      zlim <- range(attr(x,"pattern")[,,ip], na.rm=TRUE)
    } else if(length(dim(attr(x,"pattern")))==2) {
      zlim <- range(attr(x,"pattern")[,ip], na.rm=TRUE)
    } else {
      zlim <- range(attr(x,"pattern"), na.rm=TRUE)
    }
    digs <- ceiling(max(c(1,1-log10(max(abs(zlim))))))
    zlim <- round(zlim, digits=digs)
    if(all(zlim>=0) | all(zlim<=0)) {
      if(is.null(colbar$breaks)) colbar$breaks <- pretty(zlim, n=10)
    } else {
      if(is.null(colbar$breaks)) colbar$breaks <- pretty(c(-1,1)*max(abs(zlim)), n=10)
      if(is.null(colbar$pal)) colbar$pal <- "t2m"
    }
    if (inherits(x,'eof')) {  ## inherits(x,'pca') |
      par(fig=c(0,0.5,0.5,1),mar=c(3,3,2,2))
      ## par(fig=c(0.025,0.5,0.5,0.975)) ## c(0,0.45,0.5,0.975) c(0.05,0.5,0.55,0.95)
      map(x,ip=ip,verbose=verbose,
          cex.main=cex.main,cex.axis=cex.axis,
          cex.lab=cex.lab,cex=cex,new=FALSE,colbar=colbar,...) 
    } else if (inherits(x,'pca')) {
      #par(fig=c(0.5,1,0.5,1),mar=c(3,3,2,2))
      fig <-c(0,0.5,0.5,1)
      main1 <- paste0('Leading EOF#',ip, ' (',
                      round(var.eof[ip],digits=2),"%)")
      if (inherits(x,'station')) 
        map(x,ip=ip,verbose=verbose,
            cex.main=cex.main,cex.axis=cex.axis,
            cex.lab=cex.lab,cex=cex, fig=fig,
            new=FALSE,colbar=colbar,...)  else
              if (inherits(x,'radiosonde')) {
                if (verbose) print('plot PCA of radiosonde data')
                par(fig=fig)
                plot(attr(x,'pattern')[,1],alt(x),type='l',col='grey',xlab='')
                np <- dim(attr(x,'pattern'))[1]
                for (i in 2:np) lines(attr(x,'pattern')[,i],alt(x),col='grey')
                lines(attr(x,'pattern')[,ip],alt(x))
              }
      title(main=src(x)[1],cex.main=cex.main*0.8,
            col.main="grey40",adj=0,line=0)
      title(main=main1,cex.main=cex.main)
      par(xaxt='s',yaxt='s',mar=c(3,3,2,2))
    }
  }
  ##  if (length(grep('pc',what))>0) result <- as.station(x) else
  #  if (length(grep('var',what))>0) result <- attr(x,'tot.var')
  
  ylab <- paste("PC",n)
  main <- paste('First',n,"leading EOFs: ", ## attr(x,'longname')
                round(sum(var.eof[1:n]),1),"% of variance")
  
  if (length(grep('var',what))>0) {
    par(new=TRUE,fig=c(0.5,1,0.5,1),mar=c(3,3,2,2))##,xaxt="s",yaxt="s")fig=c(0.5,0.95,0.5,0.975) 
    plot.eof.var(x,ip=ip,new=FALSE,cex.main=cex.main,
                 cex.axis=cex.axis,bty="n",cex=cex)
  }
  
  #print(main)
  if (length(grep('pc',what))>0) {
    ##par(bty="n", ##,xaxt="s",yaxt="s",xpd=FALSE,
    par(fig=c(0.05,1,0.025,0.475),mar=c(3,3,2,2),new=TRUE) ##,cex.axis=0.9,cex.lab=1) ##(0.05,0.95,0.02,0.45)
    main <- paste0('Leading PC#',ip,' of ',attr(x,'longname'),
                   " - Explained variance = ",round(var.eof[ip],digits=2),"%")
    if(inherits(x,"seasonalcycle")) xaxt <- "n" else  xaxt <- NULL
    xn <- x[,n]
    if(inherits(index(xn),"PCICt")) {
      # KMP 2019-05-25: To handle data with PCICt format time index (special calendar data)
      # works but the date format on the x-axis sometimes looks weird...
      caldays <- as.numeric(substr(attr(x,"calendar"),1,3))
      index(xn) <- as.numeric(format(index(x),"%Y")) + 
        (as.numeric(format(index(x),"%j"))+as.numeric(format(index(x),"%H"))/24)/caldays
    }
    plot.zoo(xn,#x[,n],
             lwd=2,ylab=ylab,main=main,xlim=xlim,ylim=ylim,
             cex.main=cex.main,bty="n",cex.axis=cex.axis,
             cex.lab=cex.lab,xaxt=xaxt)
    if(inherits(x,"seasonalcycle")) axis(1,at=seq(1,12),labels=month.abb,
                                         cex.axis=cex.axis,las=2)
    grid()
  }
  
  # par(fig=c(0,1,0,0.55),new=TRUE, mar=c(1,1,1,1),xaxt="n",yaxt="n",bty="n")
  # plot(c(0,1),c(0,1),type="n",xlab="",ylab="")
  # 
  # varnm <- varid(x)[1]
  # legend(0,0.83,varnm,bty="n",cex=0.8,ncol=2,text.col="grey40")
  ## Reset the graphics settings that have been changed to original
  par(fig=par0$fig, mar=par0$mar, mgp=par0$mgp, xaxt=par0$xaxt , yaxt=par0$yaxt)
  
  #par(bty="n",xaxt="n",yaxt="n",xpd=FALSE,
  #    fig=c(0,1,0.1,1),new=FALSE)
  #par(fig=c(0,1,0,0.1),new=NEW, mar=c(0,0,0,0))  
}

#' @exportS3Method
#' @export  plot.eof.comb
plot.eof.comb <- function(x,...,new=FALSE,xlim=NULL,ylim=NULL,
                          ip=1,col=c("red"),alpha=1,
                          what=c("pc","eof","var"),colbar=NULL,verbose=FALSE) {
  if (verbose) print("plot.eof.comb")
  par0 <- par()
  n <- ip
  D <- attr(x,'eigenvalues')
  tot.var <- attr(x,'tot.var')
  var.eof <- 100* D^2/tot.var
  
  if (new) dev.new()
  if (length(what)==4) {
    par(mfrow=c(2,2))
  } else if(length(what)==3) {
    layout( matrix(c(1,2,3,3), nrow=2, byrow=TRUE) )
  } else if (length(what)==2) {
    par(mfrow=c(2,1))
  }
  
  #par(cex.axis=0.75,cex.lab=0.7,cex.main=0.8)

  if (length(grep('eof',what))>0) {
    #if (!is.null(mfrow)) par(fig=c(0,0.5,0.5,1))
    map(x,ip=ip,verbose=verbose,colbar=colbar,...)
  }
  if (verbose) print('...')
  n.app <- attr(x,'n.apps')
  col <- rep(col,n.app)
  src <- rep("",n.app+1)
  src[1] <- attr(x,'source')
  ylab <- paste("PC",n)
  main <- paste("EOF: ",n,"accounts for",
                round(var.eof[n],1),"% of variance")
  if (verbose) {print(main); print(what)}
  
  if (length(grep('var',what))>0)  {
    #    par(xaxt="s",yaxt="s")
    #    plot.eof.var(x,new=FALSE,cex.main=0.7)
    #if (!is.null(mfrow)) par(new=TRUE,fig=c(0.5,1,0.5,1))##,xaxt="s",yaxt="s")fig=c(0.5,0.95,0.5,0.975) 
    plot.eof.var(x,ip=ip,new=FALSE,cex.main=0.8,cex.axis=0.9,bty="n",verbose=verbose)
  }
  if(verbose) {print(ylim); print(names(attributes(x))); print(n.app)}
  anms <- names(attributes(x))
  apps <- anms[grep('appendix',anms)]
  n.app <- length(apps)
  if (is.null(ylim)) {
     ylim <- range(coredata(x[,n]))
     for (i in 1:n.app) {
       if(verbose) print(apps[i])
       z <- attr(x,apps[i])
       zz <- try(coredata(z[,n]))
       if (!inherits(zz,'try-error')) ylim <- range(c(ylim,zz),na.rm=TRUE)
     }  
  }
  
  if(verbose) print(xlim)
  if (is.null(xlim)) {
    xlim <- range(index(x))
    for (i in 1:n.app) {
      z <- attr(x,apps[i])
      xlim <- range(xlim,index(z))
    }
  }
  if(is.character(xlim)) xlim <- as.Date(xlim)
  
  if (length(grep('pc',what))>0) {
    if (verbose) {print('time series');print(index(x)); print(index(attr(x,'appendix.1')))}
    if ( (sd(coredata(x)[,n])/sd(coredata(attr(x,'appendix.1'))[,n]) > 100) |
         (sd(coredata(attr(x,'appendix.1'))[,n])/sd(coredata(x)[,n]) > 100) )
      warning('plot.comb.eof: PCs have very different scales')
    #    par(bty="n",xaxt="s",yaxt="s",xpd=FALSE,
    #      fig=c(0.1,0.9,0.1,0.5),new=TRUE,cex.axis=0.6,cex.lab=0.6)
    #    plot.zoo(x[,n],lwd=2,ylab=ylab,main=main,sub=attr(x,'longname'),
    #                                          xlim=xlim,ylim=ylim)
    #if (!is.null(mfrow)) par(fig=c(0.025,1,0.025,0.475),new=TRUE) ##,cex.axis=0.9,cex.lab=1) ##(0.05,0.95,0.02,0.45)
    main <- paste0('Leading PC#',ip,' of ',attr(x,'longname'),
                   " - Explained variance = ",round(var.eof[ip],digits=2),"%")
    
    plot.zoo(x[,ip],lwd=2,ylab=ylab,main=main,xlim=xlim,ylim=ylim,
             cex.main=0.8,bty="n",cex.axis=0.9,cex.lab=1,xaxt="n")
    taxis <- range(index(x))

    ## Plot the common PCs
    for (i in 1:n.app) {
      z <- attr(x,apps[i])
      zz <- try(z[,ip])
      if (!inherits(zz,'try-error')) {
        if (verbose) {print(apps[i]); print(c(dim(z),ip))}
        lines(zz,col=adjustcolor(col[i],alpha.f=alpha),lwd=2)
        taxis <- range(c(taxis, index(zz)))
      }
      if (verbose) print(attr(z,'source'))
      if (!is.null(attr(z,'source'))) src[i+1] <- attr(z,'source') else
        src[i+1] <- paste('x',i,sep='.')
    }
    
    taxis <- pretty(taxis, n=10)
    if (min(diff(taxis))> 360) taxisl <- year(taxis)  else
      taxisl <- taxis      # REB 2016-03-03
    if (verbose) print(taxisl)
    axis(1,at=taxis,labels=taxisl,cex.axis=0.9)      # REB 2016-03-03
    grid()
    
    
    lines(x[,ip],lwd=2,col="black")
  }
  #    par(xaxt="n",yaxt="n",bty="n",fig=c(0,1,0,0.1),
  #        mar=rep(0,4),new=TRUE)
  #    plot(c(0,1),c(0,1),type="n",xlab="",ylab="")
  #    legend(0,1,src,col=c("black",col),lwd=2,ncol=4,bty="n",cex=0.7)
  #    par(xaxt="n",yaxt="n",bty="n",fig=par0$fig,mar=par0$mar,new=TRUE)
  #    plot.zoo(x[,n],type="n",xlab="",ylab="")
  
  #par(fig=c(0,1,0,0.55),new=TRUE, mar=c(0,0,0,0),xaxt="n",yaxt="n",bty="n")
  #plot(c(0,1),c(0,1),type="n",xlab="",ylab="")
  varnm <- varid(x)
  #legend(0,0.83,varnm,bty="n",cex=0.8,ncol=2,text.col="grey40")
  legend("bottomright",varnm,bty="n",cex=0.8,ncol=2,text.col="grey40")
  
  #par(bty="n",xaxt="n",yaxt="n",xpd=FALSE,
  #    fig=c(0,1,0.1,1),new=TRUE)
  #par(fig=c(0,1,0,0.1),new=TRUE, mar=c(0,0,0,0))  
  
  ## Reset the graphics settings that have been changed to original
  par(fig=par0$fig, mar=par0$mar, mgp=par0$mgp, xaxt=par0$xaxt , yaxt=par0$yaxt)
}

#' Plot esd objects
#' 
#' These plot functions are S3 methods for esd objects, based on \code{plot}.
#'
#' @importFrom graphics par grid segments text axis legend polygon mtext abline layout rect boxplot
#' @importFrom grDevices dev.new dev.off dev.copy2eps dev.copy2pdf dev.list
#' @importFrom stats pbinom
#'
#' @aliases plot.ds.eof plot.ds.pca
#'
#' @param x the object to be plotted
#' @param plot.type "single"
#' @param new if TRUE plot in new window
#' @param lwd width of line
#' @param type type of plot: 'l' = line, 'p' = point, 'b' = both
#' @param pch type of marker
#' @param main main title
#' @param xlab label of x-axis
#' @param ylab label of y-axis
#' @param errorbar if TRUE show errorbar
#' @param legend.show if TRUE show legendp
#' @param map.show show map of stations
#' @param map.type 'points' to show stations on map, 'rectangle' to show area
#' @param map.insert if TRUE show map as insert, else show map in new window
#' @param it For subsetting in time - See \code{\link{subset}}.
#' @param is For subsetting in space - See \code{\link{subset}}. Can also be
#' a station value and if provided, the plotting will involve an interpolation
#' to the same coordinates as defined by \code{is}.
#' @param ip Which EOF/CCA pattern (mode) to plot
#' @param cex magnification factor, see \code{\link[graphics]{par}}
#' @param cex.axis see \code{\link[graphics]{par}}
#' @param cex.lab see \code{\link[graphics]{par}}
#' @param cex.main see \code{\link[graphics]{par}}
#' @param mar see \code{\link[graphics]{par}}
#' @param fig coordinates of figure region, see \code{\link[graphics]{par}}
#' @param alpha transparency factor for main plot
#' @param alpha.map transparency factor for map
#' @param verbose a boolean; if TRUE print information about progress
#' @param col Colour see \code{\link[graphics]{par}}
#' @param lwd width of line
#' @param xlim range of x-axis
#' @param ylim range of y-axis
#' @param what Indicate what to plot. For \code{plot.eof}, c('pc', 'eof', 'var') is the default setting which means that
#' the plot will include the principle components, EOF patterns and explained variance. 'field' expands eof to field before
#' plotting
#' @param colbar a list, see \code{\link{colbar}} 
#' @param \dots additional arguments
#' 
#' @return None
#'
#' @keywords plot graphics
#'
#' @examples
#' 
#' # Example: use aggregate to compute annual mean temperature for Svalbard:
#' data(Svalbard)
#' y <- aggregate(Svalbard, by=year(Svalbard), FUN='mean', na.rm = FALSE) 
#' plot(y, new=FALSE)
#' 
#' # Example with downscaling:
#' lon <- c(-12,37)
#' lat <- c(52,72)
#' t2m <- t2m.DNMI(lon=lon,lat=lat)
#' data(Oslo)
#' ds <- DS(Oslo,t2m)
#' 
#' # Plot the results for January month
#' # plot(subset(ds,it='Jan'))
#' 
#' # Plot the residuals:
#' residual <- as.residual(ds)
#' obs <- as.anomaly(as.calibrationdata(ds))
#' 
#' plot.zoo(obs[,1], lwd=2, new=FALSE)
#' lines(residual, col="red")
#' 
#' print("Global climate model simulation NorESM")
#' T2m <- t2m.NorESM.M(lon=lon,lat=lat)
#' 
#' # Plot the global mean of the field:
#' plot(T2m, new=FALSE)
#' # Plot area mean of a sub region
#' plot(T2m,is=list(lon=c(0,10),lat=c(60,70)), new=FALSE)
#'
#' # Plot interpolated results corresponding to ferder
#' data(ferder)
#' plot(T2m,ferder, new=FALSE)
#'
#' print("Extract a subset - the January month")
#' x <- subset(t2m,it="jan")
#' X <- subset(T2m,it="jan")
#' 
#' print("Combine the fields for computing common EOFs:")
#' XX <- combine(x,X)
#' 
#' print("Compute common EOFs")
#' eofxx <- EOF(XX)
#' plot(eofxx, new=FALSE)
#' 
#' print("Downscale the January mean temperature") 
#' ds.jan <- DS(Oslo,eofxx)
#' plot(ds.jan, new=FALSE)
#'
#' @exportS3Method
#' @export plot.ds
plot.ds <- function(x,...,plot.type="multiple",what=NULL,new=TRUE,
                    lwd=1,type='l',pch=0,main=NULL,col=NULL,
                    colbar=list(pal=NULL,rev=FALSE,n=10,
                                breaks=NULL,type="p",cex=2,show=TRUE,
                                h=0.6, v=1,pos=0.05),
                    xlim=NULL,ylim=NULL,xlab="",ylab=NULL,verbose=FALSE) {
  if (verbose) {print(paste('plot.ds')); print(names(attributes(ds)))}
  
  if (inherits(x,'pca')) {
    plot.ds.pca(x,what=what,verbose=verbose,new=new,...)
    return()
  } else if (inherits(x,'eof')) {
    plot.ds.eof(x,verbose=verbose,...)
    return()
    #}
    ## KMP 2021-11-29: What's the intention here? This redirects plot.ds
    ## to plot.ds.station.pca even for single stations (not PCA based downscaling). 
  } else if ( (inherits(x,'station')) & (inherits(attr(x,'eof'),'pca')) ) {
    plot.ds.station.pca(x,verbose=verbose,...)
    return()
  } 
  
  if(is.null(what)) what <- c("map","ts",'xval')
  if (verbose) print(paste('plot.ds',paste(what,collapse=',')))
  par0 <- par()
  
  unit <- attr(x,'unit')
  if ( (is.na(unit) | is.null(unit)) ) unit <- " "
  for (i in 1:length(unit)) {
    if ((unit[i]=='degree Celsius') | (unit[i]=='deg C') | (unit[i]=='degC'))
      unit[i] <- 'degree*C'
  }
  
  if (is.null(ylab))
    ylab <- esd::ylab(x)
  
  if (verbose)  print(ylab)
  if (is.null(main)) main <- attr(x,'longname')[1]               
  if (is.null(col)) col <- rainbow(length(x[1,]))  
  
  cols <- rep("blue",100)
  model <- attr(x,'model')
  
  if (length(what)==2) mfrow <- c(2,1) else
    if (length(what)==1) mfrow <- c(1,1)
  
  if (new) dev.new()
  if (plot.type=="single") new <- TRUE
  par(cex.axis=0.75,cex.lab=0.7,cex.main=0.8)
  
  if (sum(is.element(what,'map'))>0) {
    if (verbose) print('Show map...')
    par(fig=c(0,0.5,0.5,1))
    map(x,new=FALSE,colbar=list(show=FALSE),verbose=verbose,...)
    points(lon(x),lat(x),lwd=3,cex=1.5)
  }
  
  if ( (sum(is.element(what,'xval'))>0)  & (!is.null(attr(x,'evaluation'))) ) {
    par(new=TRUE,fig=c(0.5,1,0.5,1)) 
    plot(attr(x,'evaluation')[,1],attr(x,'evaluation')[,2],
         main='Cross-validation',xlab='original data',
         ylab='prediction',pch=19,col="grey")
    lines(range(c(attr(x,'evaluation')),na.rm=TRUE),
          range(c(attr(x,'evaluation')),na.rm=TRUE),lty=2)
    cal <- data.frame(y=coredata(attr(x,'evaluation')[,1]),
                      x=coredata(attr(x,'evaluation')[,2]))
    xvalfit <- lm(y ~ x, data = cal)
    abline(xvalfit,col=rgb(1,0,0,0.3),lwd=2)
    #par(bty="n",fig=c(0.6,0.95,0.48,0.52),mar=c(0,0,0,0),new=TRUE,
    #    xaxt='n',yaxt='n',cex.sub=0.7)
    #plot(c(0,1),c(0,1),type='n',xlab='',ylab='')
    ok <- is.finite(attr(x,'evaluation')[,1]) &
      is.finite(attr(x,'evaluation')[,2])
    text(par()$usr[1] + diff(range(par()$usr[1:2]))/24,
         par()$usr[4] - diff(range(par()$usr[3:4]))/12,
         paste('correlation=',
               round(cor(attr(x,'evaluation')[ok,1],attr(x,'evaluation')[ok,2]),2)),
         pos=4,cex=0.8,col='grey')
  }  else {
    xvalfit <- NULL
  }
  
  Y0 <- as.original.data(x)
  #print(index(Y0)); print(index(x))
  yX <- merge.zoo(Y0,x,all=FALSE)
  #print(summary(yX))
  y0 <- yX$Y0
  ## KMP 2021-02-26: Plot attribute 'fitted value' instead of coredata
  if (!is.null(attr(x,'n.apps'))) ns <- attr(x,'n.apps') else
    ns <- 0
  y.rng <- NA; x.rng <- NA
  if (ns > 0) {
    #print("Add other DS results")
    for (i in 1:ns) {
      eval(parse(text=paste("y <- attr(x,'appendix.",i,"')",sep="")))
      #print(summary(y))
      y.rng <- range(y.rng,y,na.rm=TRUE)
      x.rng <- range(x.rng,index(y),na.rm=TRUE)
    }
    ## KMP 2021-04-26: added following to solve problem when index is Date
    if(is.numeric(x.rng) & is.dates(index(x))) x.rng <- as.Date(x.rng)
  }
  
  if (sum(!is.finite(x.rng))>0) x.rng <- NULL
  if (sum(!is.finite(y.rng))>0) y.rng <- NULL
  
  if (is.null(ylim)) {
    #ylim <- range(coredata(x),coredata(y0),y.rng,na.rm=TRUE)
    ylim <- range(attr(x,"fitted_values"),coredata(y0),y.rng,na.rm=TRUE)
  }
  if (is.null(xlim)) {
    xlim <- range(index(x),index(y0),x.rng,na.rm=TRUE)
  }
  
  #par(fig=c(0.025,1,0.025,0.475),new=TRUE)
  par(bty="n",fig=c(0,1,0.1,0.5),mar=c(1,4.5,1,1),new=TRUE, xaxt='s',yaxt='s')
  ds <- list(obs=y0)
  ## REB 2022-08-10 testing for sensible ranges
  if (sum(!is.finite(xlim))>0) xlim <- NULL
  if (sum(!is.finite(ylim))>0) ylim <- NULL
  plot.zoo(y0,plot.type=plot.type,ylab=ylab,xlab=xlab,
           main=main,xlim=xlim,ylim=ylim,lwd=1,type='b',pch=19)
  grid()
  if (verbose) print(c(class(index(x)),class(index(y0))))
  if ( (class(index(x))=='Date') & (class(index(y0))=='numeric') & inherits(x,'annual') ) 
    index(x) <- year(index(x))
  if ( (class(index(x))=='numeric') & (class(index(y0))=='Date') & inherits(x,'annual') ) 
    index(x) <- as.Date(paste(index(x),'01-01',sep='-'))
  #lines(x,col="red",type="l",lwd=lwd)
  lines(attr(x,"fitted_values"),col="red",type="l",lwd=lwd)
  
  cal0 <- data.frame(y=coredata(y0),t=year(y0))
  #cal1 <- data.frame(y=coredata(x),t=year(x))
  cal1 <- data.frame(y=attr(x,"fitted_values"),t=year(x))
  
  trend0 <- lm(y ~ t, data=cal0)
  trend1 <- lm(y ~ t, data=cal1)
  lines(zoo(predict(trend0),order.by=index(y0)),lty=2)
  lines(zoo(predict(trend1),order.by=index(x)),lty=2,col='red')
  
  st0 <- summary(trend0); st1 <- summary(trend1)
  obstrend <- paste('obs. trend: ', round(st0$coefficients[2],2),' (',
                    round(st0$coefficients[2]-2*st0$coefficients[4],2),',',
                    round(st0$coefficients[2]+2*st0$coefficients[4],2),')',
                    attr(x,'unit'),'/decade',sep='')
  dstrend <- paste('obs. trend: ', round(st1$coefficients[2],2),' (',
                   round(st1$coefficients[2]-2*st1$coefficients[4],2),',',
                   round(st1$coefficients[2]+2*st1$coefficients[4],2),')',
                   attr(x,'unit'),'/decade',sep='')
  
  if (is.null(attr(x,'source'))) attr(x,'source') <- 'ESD'
  if (is.na(attr(x,'source'))) attr(x,'source') <- 'ESD'
  legtext <- c("Observations",attr(x,'source')) 
  legcol <- c("black","red")
  
  if (ns > 0) {
    #print("Add other DS results")
    for (i in 1:ns) {
      eval(parse(text=paste("y <- attr(x,'appendix.",i,"')",sep="")))
      lines(zoo(coredata(y),order.by=index(y)),col=cols[i],lwd=lwd)
      legcol <- c(legcol,cols[i])
      legtext <- c(legtext,attr(y,'source'))
      eval(parse(text=paste("ds$result.",i,
                            " <- attr(x,'appendix.",i,"')",sep="")))
    }
    
  } else  legcol <- c("black","red")
  
  ## Replot observations and prediction for calibration period
  
  lines(y0,lwd=1,type='b',pch=19)
  #lines(x,col="red",type="l",lwd=lwd)
  lines(attr(x,"fitted_values"),col="red",type="l",lwd=lwd)
  #print(legcol)
  if (!is.null(attr(x,'appendix.1'))) legend <- c("Obs.","Cal.","Proj") else
    legend <- c("Obs.","Cal.")
  legend(x="topleft",legend=legend,bty="n",horiz=TRUE,
         col=c("black","red","blue"),lwd=c(1,1,1),pch=c(19,1,1))
  
  if (plot.type=="single") {
    par(fig=c(0,1,0,0.1),new=TRUE, mar=c(0,0,0,0),xaxt="s",yaxt="s",bty="n")
    plot(c(0,1),c(0,1),type="n",xlab="",ylab="")
    legend(0.01,0.95,c(paste(attr(x,'location'),": ",
                             #attr(x,'aspect'),
                             #attr(x,'longname')," - ",
                             round(attr(x,'longitude'),2),"E/",
                             round(attr(x,'latitude'),2),"N (",
                             attr(x,'altitude')," masl)",sep=""),
                       obstrend,dstrend),ncol=3,
           bty="n",cex=0.6,text.col="grey40",
           lwd=c(1,rep(2,ns),1,1),lty=c(1,rep(1,ns),2,2),
           col=c(col,'black',col[1]))
    
    par(fig=c(0,1,0.05,0.95),new=TRUE,mar=par0$mar,xaxt="n",yaxt="n",bty="n")
    #plot.zoo(x,plot.type=plot.type,type="n",ylab="",xlab="",xlim=xlim,ylim=ylim)
    plot.zoo(attr(x,"fitted_values"),plot.type=plot.type,type="n",
             ylab="",xlab="",xlim=xlim,ylim=ylim)
  }
  par(fig=par0$fig, mar=par0$mar, mgp=par0$mgp, xaxt=par0$xaxt, yaxt=par0$yaxt)
  invisible(list(trend0=trend0,trend1=trend1,xvalfit=xvalfit))
}

#' @exportS3Method
#' @export plot.eof.var
plot.eof.var <- function(x,...,ip=1,new=TRUE,xlim=NULL,ylim=NULL,n=20,verbose=FALSE) {
  if(verbose) print("plot.eof.var")
  n <- min(c(n,length(attr(x,'eigenvalues'))))
  D <- attr(x,'eigenvalues')
  tot.var <- attr(x,'tot.var')
  var.eof <- 100* D^2/tot.var
  nt <- length(index(x)) 
  n.eff <- round(nt * (1.0-attr(x,'max.autocor'))/
                   (1.0+attr(x,'max.autocor')))  
  dD <- D*sqrt(2.0/n.eff)
  main <- paste(n,"leading EOFs: ", ## attr(x,'longname')
                round(sum(var.eof[1:n]),digits=1),"% of variance")
  ##main <- paste(attr(x,'longname'),n,"leading EOFs: ",
  ##               round(sum(var.eof[1:n]),1),"% of variance")
  if (is.null(xlim)) xlim <- c(0,n) ##c(0.7,n+0.3)
  if (is.null(ylim)) ylim <- c(0,100)
  if (new) dev.new()
  ##par(bty="n") ##,xaxt="n")
  plot(var.eof,type="n",
       main=main,ylab="Variance (%)",xlab="EOF order",
       xlim=xlim,ylim=ylim,...)
  lines(cumsum(var.eof),lwd=3,col=rgb(1,0.8,0.8))
  grid() ## nx=21,ny=12)
  lines(var.eof,type="b",pch=19)
  for (i in 1:length(var.eof)) {
    lines(rep(i,2),100*c((D[i]+dD[i])^2/tot.var,(D[i]-dD[i])^2/tot.var),
          lty=1,col="darkgrey")
    lines(c(i-0.25,i+0.25),100*rep((D[i]+dD[i])^2/tot.var,2),
          lwd=1,col="darkgrey")
    lines(c(i-0.25,i+0.25),100*rep((D[i]-dD[i])^2/tot.var,2),
          lwd=1,col="darkgrey")
  }
  points(var.eof,cex=1.5)
  points(var.eof,pch=20,cex=1.2,col="darkgrey")
  points(ip,var.eof[ip],pch=20,col="red",cex=1.2)
  attr(var.eof,'errorbar') <- cbind(100*(D-dD)^2/tot.var,100*(D+dD)^2/tot.var)
  invisible(var.eof)
}

#' Plot esd objects
#' 
#' These plot functions are S3 methods for esd objects, based on \code{plot}.
#'
#' @importFrom graphics par grid segments text axis legend polygon mtext abline layout rect boxplot
#' @importFrom grDevices dev.new dev.off dev.copy2eps dev.copy2pdf dev.list
#' @importFrom stats pbinom
#'
#' @param x the object to be plotted
#' @param plot.type "single"
#' @param new if TRUE plot in new window
#' @param lwd width of line
#' @param type type of plot: 'l' = line, 'p' = point, 'b' = both
#' @param pch type of marker
#' @param main main title
#' @param xlab label of x-axis
#' @param ylab label of y-axis
#' @param errorbar if TRUE show errorbar
#' @param legend.show if TRUE show legendp
#' @param map.show show map of stations
#' @param map.type 'points' to show stations on map, 'rectangle' to show area
#' @param map.insert if TRUE show map as insert, else show map in new window
#' @param it For subsetting in time - See \code{\link{subset}}.
#' @param is For subsetting in space - See \code{\link{subset}}. Can also be
#' a station value and if provided, the plotting will involve an interpolation
#' to the same coordinates as defined by \code{is}.
#' @param ip Which EOF/CCA pattern (mode) to plot
#' @param cex magnification factor, see \code{\link[graphics]{par}}
#' @param cex.axis see \code{\link[graphics]{par}}
#' @param cex.lab see \code{\link[graphics]{par}}
#' @param cex.main see \code{\link[graphics]{par}}
#' @param mar see \code{\link[graphics]{par}}
#' @param fig coordinates of figure region, see \code{\link[graphics]{par}}
#' @param alpha transparency factor for main plot
#' @param alpha.map transparency factor for map
#' @param verbose a boolean; if TRUE print information about progress
#' @param col Colour see \code{\link[graphics]{par}}
#' @param lwd width of line
#' @param xlim range of x-axis
#' @param ylim range of y-axis
#' @param what Indicate what to plot. For \code{plot.eof}, c('pc', 'eof', 'var') is the default setting which means that
#' the plot will include the principle components, EOF patterns and explained variance. 'field' expands eof to field before
#' plotting
#' @param colbar a list, see \code{\link{colbar}} 
#' @param \dots additional arguments
#' 
#' @return None
#'
#' @keywords plot graphics
#'
#' @examples
#' 
#' # Example: use aggregate to compute annual mean temperature for Svalbard:
#' data(Svalbard)
#' y <- aggregate(Svalbard, by=year(Svalbard), FUN='mean', na.rm = FALSE) 
#' plot(y, new=FALSE)
#' 
#' # Example with downscaling:
#' lon <- c(-12,37)
#' lat <- c(52,72)
#' t2m <- t2m.DNMI(lon=lon,lat=lat)
#' data(Oslo)
#' ds <- DS(Oslo,t2m)
#' 
#' # Plot the results for January month
#' # plot(subset(ds,it='Jan'))
#' 
#' # Plot the residuals:
#' residual <- as.residual(ds)
#' obs <- as.anomaly(as.calibrationdata(ds))
#' 
#' plot.zoo(obs[,1], lwd=2, new=FALSE)
#' lines(residual, col="red")
#' 
#' print("Global climate model simulation NorESM")
#' T2m <- t2m.NorESM.M(lon=lon,lat=lat)
#' 
#' # Plot the global mean of the field:
#' plot(T2m, new=FALSE)
#' # Plot area mean of a sub region
#' plot(T2m,is=list(lon=c(0,10),lat=c(60,70)), new=FALSE)
#' 
#' # Plot interpolated results corresponding to ferder
#' data(ferder)
#' plot(T2m, is=ferder, new=FALSE)
#' 
#' # Plot Hovmuller diagram: Not working ...
#' ## plot(T2m,is=list(lon=0)) 
#' 
#' print("Extract a subset - the January month")
#' x <- subset(t2m,it="jan")
#' X <- subset(T2m,it="jan")
#' 
#' print("Combine the fields for computing common EOFs:")
#' XX <- combine(x,X)
#' 
#' print("Compute common EOFs")
#' eofxx <- EOF(XX)
#' plot(eofxx, new=FALSE)
#' 
#' print("Downscale the January mean temperature") 
#' ds.jan <- DS(Oslo,eofxx)
#' plot(ds.jan, new=FALSE)
#'
#' @exportS3Method
#' @export plot.field
plot.field <- function(x,...,is=NULL,it=NULL,FUN="mean",map.type='rectangle',verbose=FALSE) {
  if (verbose) print("plot.field")
  stopifnot(!missing(x),inherits(x,'field'))
  
  d <- dim(x)
  if (d[2]==1) {
    if (verbose) print('one grid point')
    plot.station(x,verbose=verbose,...)
    return()
  }
  
  # To distinguish between types of plots - movmuller or are mean
  twocoord <- (length(is)==2)
  if (twocoord) {
    l1 <- length(is[[1]]); l2 <- length(is[[2]])
  } else {
    l1 <- 0; l2 <- 0
  }
  #print(c(l1,l2))
  if (inherits(is,'station')) {
    # If a station object is provided - extract a time series for a
    # the corresponding location
    lon <- attr(is,'longitude')
    lat <- attr(is,'latitude')
  } else if ((inherits(x,'list')) & twocoord & (l1==1) & (l2==1) ) {
    # Is spatial coordinates  
    lon <- is[[1]]
    lat <- is[[2]]
    z <- NULL
  } else if (!is.null(is)) {
    nms <- names(is)
    lon <- attr(x,'longitude')
    lat <- attr(x,'latitude')
    #print(nms)
    if (length(nms)==2) {
      lon <- is[[1]]
      lat <- is[[2]]
      y <- subset(x,is=is,it=it)
      z <- aggregate.area(y,FUN=FUN)
    } else if ( (length(nms)==1) & (tolower(nms)=="lon") ) {
      # Hovmuller diagram along latitude
      #print(is[[1]]); print(lon); print(d)
      if (length(is[[1]])== 1) {
        picklon <- lon[max( (1:length(lon))[lon <= is[[1]]] )]
        #print(picklon)
        xy <- rep(lon,length(lat))
        yx <- sort(rep(lat,length(lon)))
        ix <- is.element(xy,picklon)
        z <- x[,ix]
        z <- attrcp(x,z)
        #image(z)
        attr(z,'longitude') <- picklon
        attr(z,'latitude') <- attr(x,'latitude')
        #print(attr(x,'dimensions'))
        attr(z,'dimensions') <- c(1,attr(x,'dimensions')[2:3])
        class(z) <- c('xsection',class(x)[-1])
      } else if (length(is[[1]])== 2) {
        y <- subset(x,is=list(lon=is[[1]],lat=range(lat)))
        xy <- rep(attr(x,'longitude'),d[2])
        yx <- sort(rep(attr(x,'latitude'),d[1]))
        X <- as.data.frame(coredata(x)); colnames(yx)
        Z <- aggregate(x,by=yx,FUN=FUN)
        z <- zoo(Z,order.by=index(x))
        z <- attrcp(x,z)
        class(z) <- class(x)
      }
    } else if ( (length(nms)==1) & (tolower(nms)=="lat") ) {
      # Hovmuller diagram along longitude
      if (length(is[[1]])== 1) {
        picklat <- lat[max( (1:length(lat))[lat <= is[[1]]] )]
        xy <- rep(attr(x,'longitude'),length(lat))
        yx <- sort(rep(attr(x,'latitude'),length(lon)))
        iy <- is.element(yx,picklat)
        z <- x[,iy]
        z <- attrcp(x,z)
        attr(z,'latitude') <- picklat
        attr(z,'longitude') <- attr(x,'longitude')
        attr(z,'dimensions') <- c(attr(x,'dimensions')[1],1,attr(x,'dimensions')[3])
        class(z) <- c('xsection',class(x)[-1])
      } else if (length(is[[1]])== 2) {
        y <- subset(x,is=list(lon=range(lon),lat=is[[1]]))
        xy <- rep(attr(x,'longitude'),d[2])
        yx <- sort(rep(attr(x,'latitude'),d[1]))
        X <- as.data.frame(coredata(x)); colnames(xy)
        Z <- aggregate(x,by=xy,FUN=FUN)
        z <- zoo(Z,order.by=index(x))
        class(z) <- class(x)
      }
    }
    
  } else {lon <- NULL; lat <- NULL}
  
  if ( is.null(lon) & is.null(lat) ) {
    #print("aggregate")
    z <- aggregate.area(x,is=is,FUN=FUN)
    class(z) <- c('station',class(x)[-1])
  } else if ( (length(lon)==1) & (length(lat)==1) ) {
    #print("regrid")
    z <- regrid(x,list(x=lon,y=lat))
    class(z) <- c('station',class(x)[-1])
  }
  #print("plot")
  plot(z,map.type=map.type,...)
  z <- attrcp(x,z,ignore=c("longitude","latitude"))
  attr(z,'history') <- history.stamp(x)
  if (inherits(x,'station')) lines(y,col="red",lwd=2)
  invisible(z)
}

#' Plot esd objects
#' 
#' These plot functions are S3 methods for esd objects, based on \code{plot}.
#'
#' @importFrom graphics par grid segments text axis legend polygon mtext abline layout rect boxplot
#' @importFrom grDevices dev.new dev.off dev.copy2eps dev.copy2pdf dev.list
#' @importFrom stats pbinom
#'
#' @param x the object to be plotted
#' @param plot.type "single"
#' @param new if TRUE plot in new window
#' @param lwd width of line
#' @param type type of plot: 'l' = line, 'p' = point, 'b' = both
#' @param pch type of marker
#' @param main main title
#' @param xlab label of x-axis
#' @param ylab label of y-axis
#' @param errorbar if TRUE show errorbar
#' @param legend.show if TRUE show legendp
#' @param map.show show map of stations
#' @param map.type 'points' to show stations on map, 'rectangle' to show area
#' @param map.insert if TRUE show map as insert, else show map in new window
#' @param it For subsetting in time - See \code{\link{subset}}.
#' @param is For subsetting in space - See \code{\link{subset}}. Can also be
#' a station value and if provided, the plotting will involve an interpolation
#' to the same coordinates as defined by \code{is}.
#' @param ip Which EOF/CCA pattern (mode) to plot
#' @param cex magnification factor, see \code{\link[graphics]{par}}
#' @param cex.axis see \code{\link[graphics]{par}}
#' @param cex.lab see \code{\link[graphics]{par}}
#' @param cex.main see \code{\link[graphics]{par}}
#' @param mar see \code{\link[graphics]{par}}
#' @param fig coordinates of figure region, see \code{\link[graphics]{par}}
#' @param alpha transparency factor for main plot
#' @param alpha.map transparency factor for map
#' @param verbose a boolean; if TRUE print information about progress
#' @param col Colour see \code{\link[graphics]{par}}
#' @param lwd width of line
#' @param xlim range of x-axis
#' @param ylim range of y-axis
#' @param what Indicate what to plot. For \code{plot.eof}, c('pc', 'eof', 'var') is the default setting which means that
#' the plot will include the principle components, EOF patterns and explained variance. 'field' expands eof to field before
#' plotting
#' @param colbar a list, see \code{\link{colbar}} 
#' @param \dots additional arguments
#' 
#' @return None
#'
#' @keywords plot graphics
#'
#' @examples
#' 
#' # Example: use aggregate to compute annual mean temperature for Svalbard:
#' data(Svalbard)
#' y <- aggregate(Svalbard, by=year(Svalbard), FUN='mean', na.rm = FALSE) 
#' plot(y, new=FALSE)
#' 
#' # Example with downscaling:
#' lon <- c(-12,37)
#' lat <- c(52,72)
#' t2m <- t2m.DNMI(lon=lon,lat=lat)
#' data(Oslo)
#' ds <- DS(Oslo,t2m)
#' 
#' # Plot the results for January month
#' # plot(subset(ds,it='Jan'))
#' 
#' # Plot the residuals:
#' residual <- as.residual(ds)
#' obs <- as.anomaly(as.calibrationdata(ds))
#' 
#' plot.zoo(obs[,1], lwd=2, new=FALSE)
#' lines(residual, col="red")
#' 
#' print("Global climate model simulation NorESM")
#' T2m <- t2m.NorESM.M(lon=lon,lat=lat)
#' 
#' # Plot the global mean of the field:
#' plot(T2m, new=FALSE)
#' # Plot area mean of a sub region
#' plot(T2m,is=list(lon=c(0,10),lat=c(60,70)), new=FALSE)
#' 
#' # Plot interpolated results corresponding to ferder
#' data(ferder)
#' plot(T2m,ferder, new=FALSE)
#' 
#' # Plot Hovmuller diagram: Not working ...
#' ## plot(T2m,is=list(lon=0)) 
#' 
#' print("Extract a subset - the January month")
#' x <- subset(t2m,it="jan")
#' X <- subset(T2m,it="jan")
#' 
#' print("Combine the fields for computing common EOFs:")
#' XX <- combine(x,X)
#' 
#' print("Compute common EOFs")
#' eofxx <- EOF(XX)
#' plot(eofxx, new=FALSE)
#' 
#' print("Downscale the January mean temperature") 
#' ds.jan <- DS(Oslo,eofxx)
#' plot(ds.jan, new=FALSE)
#'
#' @exportS3Method
#' @export plot.pca
plot.pca <- function(x,...,ip=1,cex=1,verbose=FALSE,new=TRUE) {
  if (verbose) print('plot.pca')
  if(inherits(x,"trajectory")) {
    plot.pca.trajectory(x,cex=cex,new=new,verbose=verbose,...)
  } else {
    attr(x,'longname') <- attr(x,'longname')[1]
    if(is.null(ip)) ip <- 1
    if(length(ip)>1) {
      plot.pca.multiple(x,ip=ip,verbose=verbose,new=new,cex=cex,...)
    } else {
      plot.eof.field(x,ip=ip,verbose=verbose,new=new,cex=cex,...)
    }
  }
}

#' @export plot.pca.multiple
plot.pca.multiple <- function(x,...,new=FALSE,xlim=NULL,ylim=NULL,ip=1:3,
                              colbar=NULL,cex.axis=0.9,cex.main=0.9,cex.lab=0.9,
                              mar=c(3,3,4,2),mgp1=c(0.5,0.1,0),mgp2=c(2,0.5,0),
                              verbose=FALSE,it=NULL,is=NULL,cex=1) {
  if (verbose) print(paste('plot.pca.multiple'))
  ## Save the original graphics settings
  par0 <- par()
  n <- ip
  
  D <- attr(x,'eigenvalues')
  tot.var <- attr(x,'tot.var')
  var.eof <- 100* D^2/tot.var
  
  if (new) dev.new()
  
  for(i in ip) {
    j <- which(ip==i)
    if (verbose) {print('Show map'); print(class(x))}
    zlim <- range(attr(x,"pattern")[,i])
    digs <- ceiling(max(c(1,1-log10(max(abs(zlim))))))
    zlim <- round(zlim, digits=digs)
    if(all(zlim>=0) | all(zlim<=0)) {
      colbar$breaks <- pretty(zlim, n=10)
    } else {
      colbar$breaks <- pretty(c(-1,1)*max(abs(zlim)), n=10)
    }
    fig=c(0,1/3,1.02-j/length(ip),1-(j-1)/length(ip))
    if (inherits(x,'eof')) {
      par(fig=fig,mar=mar,mgp=c(0,0,0),new=j!=1)#TRUE)
      map(x,ip=i,verbose=verbose,
          cex.main=cex.main,cex.axis=cex.axis,
          cex.lab=cex.lab,cex=cex,new=FALSE,colbar=colbar,...)
    } else if (inherits(x,'pca')) {
      main1 <- paste0('EOF ',i)
      if(which(i==ip)==1) xlab1 <- "lon" else xlab1 <- " "
      if(which(i==ip)==1) ylab1 <- "lat" else ylab1 <- " "
      map(x,ip=i,verbose=verbose,mar=mar,mgp=mgp1,
          cex.main=cex.main,cex.axis=cex.axis,
          cex.lab=cex.lab,cex=cex,fig=fig,add=j!=1,
          xlab=xlab1,ylab=ylab1,new=FALSE, colbar=colbar,
          main=" ",...)
      #title(main=main1,cex.main=cex.main)
      mtext(main1, side=3, adj=0.1, line=-1, cex=cex.main)
      par(xaxt='s',yaxt='s')
    }
    
    ylab <- paste("PC",i)
    if(which(ip==i)==1) {
      main <- paste0('PC ',i,' of ',attr(x,'longname')[1],
                     " - Explained variance = ",round(var.eof[i],digits=2),"%")
    } else {
      main <- paste0('PC ',i,
                     " - Explained variance = ",round(var.eof[i],digits=2),"%")
    }
    if(inherits(x,"seasonalcycle")) xaxt <- "n" else  xaxt <- NULL
    xi <- x[,i]
    if(inherits(index(xi),"PCICt")) {
      # KMP 2019-05-25: To handle data with PCICt format time index (special calendar data)
      # works but the date format on the x-axis sometimes looks weird...
      caldays <- as.numeric(substr(attr(x,"calendar"),1,3))
      index(xi) <- as.numeric(format(index(x),"%Y")) +
        (as.numeric(format(index(x),"%j"))+as.numeric(format(index(x),"%H"))/24)/caldays
    }
    par(fig=c(1/3,1,1.02-j/length(ip),1-(j-1)/length(ip)),mar=mar,mgp=mgp2,new=TRUE)
    plot.zoo(xi,lwd=2,ylab=ylab,main=main,xlim=xlim,ylim=ylim,
             cex.main=cex.main,bty="n",cex.axis=cex.axis,
             cex.lab=cex.lab,xaxt=xaxt)
    if(inherits(x,"seasonalcycle")) axis(1,at=seq(1,12),labels=month.abb,
                                         cex.axis=cex.axis,las=2)
    grid()
  }
  ## Reset the graphics settings that have been changed to original
  par(fig=par0$fig, mar=par0$mar, mgp=par0$mgp, xaxt=par0$xaxt , yaxt=par0$yaxt)
}


#' @exportS3Method
#' @export plot.ds.pca
plot.ds.pca <- function(x,...,ip=1,
                        colbar1=list(pal=NULL,rev=FALSE,n=10,breaks=NULL,
                                     type="p",cex=1,show=TRUE,h=0.6, v=1,pos=0.05),
                        colbar2=NULL,mar=c(3,2.5,2,0.5),mgp=c(1.5,0.5,0),
                        what=NULL, new=FALSE, add=FALSE, verbose=FALSE) {
  y <- x # quick fix
  par0 <- par()
  if (verbose) print('plot.ds.pca')
  if (is.null(colbar2)) colbar2 <- colbar1
  attr(y,'longname') <- attr(y,'longname')[1]
  
  if(is.null(what)) what <- c("predictor.pattern","pca.pattern",
                              "xval","timeseries")
  if(is.null(attr(y,'evaluation'))) what <- what[!what %in% c("xval")]
  
  i <- 1
  if(length(what)==4) {
    figlist <- list(fig1=c(0.05,0.45,0.5,0.975), 
                    fig2=c(0.55,0.975,0.5,0.975),
                    fig3=c(0.05,0.45,0.05,0.475),
                    fig4=c(0.55,0.975,0.05,0.475))
  } else if(length(what)==3) {
    figlist <- list(fig1=c(0.05,0.45,0.5,0.975), 
                    fig2=c(0.5,0.975,0.5,0.975),
                    fig3=c(0.05,0.975,0.05,0.475))
  } else if(length(what)==2) {
    figlist <- list(fig1=c(0.05,0.45,0.05,0.975), 
                    fig2=c(0.5,0.975,0.05,0.975))
  } else figlist <- NULL

  if(new) dev.new()
  par(mar=mar, mgp=mgp)
  if('pca.pattern' %in% what) {
    if (verbose) print('PCA ip')
    map.pca(y,ip=ip,new=FALSE,colbar=colbar1,
            fig=figlist[[i]], add=(add | i>1), #fig=c(0,0.5,0.5,0.975),
            main=paste("PCA Pattern #",ip," (unitless)",sep=""),
            verbose=verbose,...)
    #title(paste("PCA Pattern # ",ip,sep=""))
    i <- i+1
  }
  
  if('predictor.pattern' %in% what) {
    if (verbose) print('Predictor pattern')
    par(new=(add | i>1))
    if(!is.null(figlist)) par(fig=figlist[[i]])
    #fig=c(0.55,0.975,0.5,0.975),new=TRUE)
    i <- i+1
    pp <- attr(y,'predictor.pattern')
    attr(pp,"variable") <- paste("EOF.Pattern.",ip,sep="")
    attr(pp,"unit") <- "unitless"
    colbar2$show <- FALSE
    map(pp,ip=ip,new=FALSE,
        colbar=colbar2,verbose=verbose,
        lab=paste("Predictor pattern # ",ip,sep=""))
    par(cex=1)  ## REB 2019-08-06
    #title(paste("Predictor pattern #",ip,sep=""), cex=0.8)
  }
  if('xval' %in% what) {
    if (verbose) print('Evaluation results')
    par(new=(add | i>1))
    if(!is.null(figlist)) par(fig=figlist[[i]])
    i <- i+1
    #par(fig=c(0.05,0.45,0.05,0.475),new=TRUE)
    ## Get the right pattern
    xvp <- (ip-1)*2 +1
    xok <- is.finite(attr(y,'evaluation')[,xvp]) & is.finite(attr(y,'evaluation')[,xvp+1])
    xcor <- cor(attr(y,'evaluation')[xok,xvp],attr(y,'evaluation')[xok,xvp+1])
    
    xlim <- range(attr(y,'evaluation')[,xvp:(xvp+1)], na.rm=TRUE)
    xlim <- xlim + c(-1,1)*diff(xlim)*0.1
    
    par(mgp=c(2,0.5,0), mar=c(3,3.5,1.5,1.5))
    par(xaxt="s", yaxt="s")
    plot(attr(y,'evaluation')[,xvp],attr(y,'evaluation')[,xvp+1],
         main=paste('Cross-validation: r=',round(xcor,2)),
         xlab='original data',ylab='prediction',pch=19,col="grey",
         xlim=xlim, ylim=xlim)
    lines(range(c(attr(y,'evaluation')),na.rm=TRUE),
          range(c(attr(y,'evaluation')),na.rm=TRUE),lty=2)
    cal <- data.frame(y=coredata(attr(y,'evaluation')[,1]),
                      x=coredata(attr(y,'evaluation')[,2]))
    xvalfit <- lm(y ~ x, data = cal)
    abline(xvalfit,col=rgb(1,0,0,0.3),lwd=2)
    #legend("bottomleft", )
  }
  if('timeseries' %in% what) {
    par(new=(add | i>1), xaxt="s", yaxt="s")
    if(!is.null(figlist)) par(fig=figlist[[i]])
    #par(fig=c(0.55,0.975,0.05,0.475),new=TRUE)
    if(is.null(attr(y, 'evaluation'))) {
      plot(attr(y,'original_data')[,ip],lwd=2,type='b',pch=19)
      lines(zoo(y[,ip]),lwd=2,col='red',type='b')
      xvalfit <- NULL
    } else {
      xlim <- range(index(attr(y,'original_data')),index(y))
      ylim <- range(attr(y,'original_data')[,ip],y[,ip],na.rm=TRUE) +
        diff(range(attr(y,'original_data')[,ip],y[,ip],na.rm=TRUE))*c(-0.1,0.2)
      y0 <- attr(y,'original_data')
      plot(y0[,ip], lwd=2, type='b', pch=19, xlim=xlim, ylim=ylim,
           xlab="Date",ylab=paste("PC #",ip,sep=""))
      grid()
      if ( (class(index(y))=='Date') & (class(index(y0))=='numeric') & inherits(y,'annual') )
        index(y) <- year(index(y))
      if ( (class(index(y))=='numeric') & (class(index(y0))=='Date') & inherits(y,'annual') )
        index(y) <- as.Date(paste(index(y),'01-01',sep='-'))
      
      lines(zoo(y[,ip]),lwd=2,col='red',type='b')
      cal0 <- data.frame(y=coredata(y0[,ip]),t=year(y0))
      #cal1 <- data.frame(y=coredata(x),t=year(x))
      cal1 <- data.frame(y=coredata(y[,ip]),t=year(y[,ip]))
      
      trend0 <- lm(y ~ t, data=cal0)
      trend1 <- lm(y ~ t, data=cal1)
      lines(zoo(predict(trend0),order.by=index(y0)),lty=2)
      lines(zoo(predict(trend1),order.by=index(x)),lty=2,col='red')
      
      legend(x=index(attr(y,'original_data')[,ip])[1],
             y=max(attr(y,'original_data')[,ip],na.rm=TRUE)+
               0.2*diff(range(attr(y,'original_data')[,ip])),
             legend=c("estimated","original"),col=c("red","black"),lty=c(1,1),
             lwd=c(2,2),pch=c(21,19),bty="n")
    }
  }
  ## Reset the graphics settings that have been changed to original
  par(fig=par0$fig, mar=par0$mar, mgp=par0$mgp, xaxt=par0$xaxt , yaxt=par0$yaxt)
}

#' @exportS3Method
#' @export plot.ds.eof
plot.ds.eof <- function(x,...,ip=1,
                        colbar1=list(pal=NULL,rev=FALSE,n=10,breaks=NULL,type="p",cex=2,show=TRUE,
                                     h=0.6, v=1,pos=0.05),colbar2=NULL,type1="fill",type2=c("fill","contour"),
                        verbose=FALSE) {
  y <- x # quick fix
  par0 <- par()
  if (verbose) print('plot.ds.eof')
  if (is.null(colbar2)) colbar2 <- colbar1 
  attr(y,'longname') <- attr(y,'longname')[1]
  par(fig=c(0,0.5,0.5,1),mar=c(3,5,4.2,1),mgp=c(3,0.5,0.5))
  map.eof(y,ip=ip,verbose=verbose,new=FALSE,colbar=colbar1,
          main=paste("(a) Predictand EOF pattern # ",ip,sep=""),type=type1,...)
  par(fig=c(0.5,1,0.5,1),mar=c(3,4,4.2,1),new=TRUE)
  map(attr(y,'predictor.pattern'),ip=ip,new=FALSE,
      colbar=colbar2,verbose=verbose,
      main=paste("(b) Predictor EOF pattern # ",ip,sep=""),type=type2)
  #title(paste("EOF Pattern # ",ip,sep=""))
  par(cex=1)  ## REB 2019-08-06
  if (!is.null(attr(y,'evaluation'))) {
    par(fig=c(0,0.5,0,0.48),mar=c(3,4.5,3,1),new=TRUE)
    pc.obs <- attr(y,'evaluation')[,1+2*(ip-1)]
    pc.ds <- attr(y,'evaluation')[,2+2*(ip-1)]
    plot(pc.obs,pc.ds,main='(c) Cross-validation',xlab='original data',
         ylab='prediction',pch=19,col="grey",
         xlim=range(pc.obs,pc.ds),ylim=range(pc.obs,pc.ds))
    lines(range(c(attr(y,'evaluation')),na.rm=TRUE),
          range(c(attr(y,'evaluation')),na.rm=TRUE),lty=2)
    cal <- data.frame(y=coredata(pc.obs),x=coredata(pc.ds))
    xvalfit <- lm(y ~ x, data = cal)
    r.xval <- round(cor(pc.obs,pc.ds),2)
    abline(xvalfit,col=rgb(1,0,0,0.3),lwd=2)
    legend("topleft", paste("r =",r.xval),bty='n')
    par(fig=c(0.5,1,0,0.48),mar=c(3,4.5,3,1),new=TRUE)
    xlim <- range(index(attr(y,'original_data')),index(y))
    ylim <- range(attr(y,'original_data')[,ip],y[,ip],na.rm=TRUE)
    y0 <- attr(y,'original_data')
    plot(y0[,ip],
         main=paste("(d) PC",ip),ylab="",
         ylim=ylim*c(1.2,1.2),xlim=xlim,
         lwd=2,type='b',pch=19)
    if ( (class(index(y))=='Date') & (class(index(y0))=='numeric') & inherits(y,'annual') ) 
      index(y) <- year(index(y))
    if ( (class(index(y))=='numeric') & (class(index(y0))=='Date') & inherits(y,'annual') ) 
      index(y) <- as.Date(paste(index(y),'01-01',sep='-'))
    lines(zoo(y[,ip]),lwd=2,col='red',type='b')
    legend("topleft",
           legend=c("estimated","original"),col=c("red","black"),lty=c(1,1),
           lwd=c(2,2),pch=c(21,19),bty="n")
  } else {
    par(fig=c(0,1,0,0.48),mar=c(3,4.5,3,1),new=TRUE)
    plot(attr(y,'original_data')[,ip],
         main="PC1",ylab="",
         ylim=range(attr(y,'original_data')[,ip])*c(1.6,1.6),
         lwd=2,type='b',pch=19)
    lines(zoo(y[,ip]),lwd=2,col='red',type='b')
    legend("topleft",
           legend=c("estimated","original"),col=c("red","black"),lty=c(1,1),
           lwd=c(2,2),pch=c(21,19),bty="n")
    xvalfit <- NULL
  }  
  ## Reset the graphics settings that have been changed to original
  par(fig=par0$fig, mar=par0$mar, mgp=par0$mgp, xaxt=par0$xaxt , yaxt=par0$yaxt)
}

#' Plot esd objects
#' 
#' The plot functions are S3 methods for esd objects, based on \code{plot}. 
#' The function \code{plot.mvr} produces a plot for an \code{mvr} object which
#' is the output of the multi-variate regression function \code{MVR},
#'
#' @seealso plot.station plot.eof plot.ds 
#'
#' @param x the object to be plotted
#' @param verbose a boolean; if TRUE print information about progress
#' @param \dots additional arguments
#' 
#' @return None
#'
#' @keywords plot graphics
#'
#' @exportS3Method
#' @export  plot.mvr
plot.mvr <- function(x,verbose=FALSE,...) {
  if(verbose) print("plot.mvr")
  plot(x$fitted.values, verbose=verbose,...)
}

#' Plot esd objects
#' 
#' These plot functions are S3 methods for esd objects, based on \code{plot}.
#'
#' @importFrom graphics par grid segments text axis legend polygon mtext abline layout rect boxplot
#' @importFrom grDevices dev.new dev.off dev.copy2eps dev.copy2pdf dev.list
#' @importFrom stats pbinom
#'
#' @param x the object to be plotted
#' @param plot.type "single"
#' @param new if TRUE plot in new window
#' @param lwd width of line
#' @param type type of plot: 'l' = line, 'p' = point, 'b' = both
#' @param pch type of marker
#' @param main main title
#' @param xlab label of x-axis
#' @param ylab label of y-axis
#' @param errorbar if TRUE show errorbar
#' @param legend.show if TRUE show legendp
#' @param map.show show map of stations
#' @param map.type 'points' to show stations on map, 'rectangle' to show area
#' @param map.insert if TRUE show map as insert, else show map in new window
#' @param it For subsetting in time - See \code{\link{subset}}.
#' @param is For subsetting in space - See \code{\link{subset}}. Can also be
#' a station value and if provided, the plotting will involve an interpolation
#' to the same coordinates as defined by \code{is}.
#' @param ip Which EOF/CCA pattern (mode) to plot
#' @param cex magnification factor, see \code{\link[graphics]{par}}
#' @param cex.axis see \code{\link[graphics]{par}}
#' @param cex.lab see \code{\link[graphics]{par}}
#' @param cex.main see \code{\link[graphics]{par}}
#' @param mar see \code{\link[graphics]{par}}
#' @param fig coordinates of figure region, see \code{\link[graphics]{par}}
#' @param alpha transparency factor for main plot
#' @param alpha.map transparency factor for map
#' @param verbose a boolean; if TRUE print information about progress
#' @param col Colour see \code{\link[graphics]{par}}
#' @param lwd width of line
#' @param xlim range of x-axis
#' @param ylim range of y-axis
#' @param what Indicate what to plot. For \code{plot.eof}, c('pc', 'eof', 'var') is the default setting which means that
#' the plot will include the principle components, EOF patterns and explained variance. 'field' expands eof to field before
#' plotting
#' @param colbar a list, see \code{\link{colbar}} 
#' @param \dots additional arguments
#' 
#' @return None
#'
#' @keywords plot graphics
#'
#' @examples
#' 
#' # Example: use aggregate to compute annual mean temperature for Svalbard:
#' data(Svalbard)
#' y <- aggregate(Svalbard, by=year(Svalbard), FUN='mean', na.rm = FALSE) 
#' plot(y, new=FALSE)
#' 
#' # Example with downscaling:
#' lon <- c(-12,37)
#' lat <- c(52,72)
#' t2m <- t2m.DNMI(lon=lon,lat=lat)
#' data(Oslo)
#' ds <- DS(Oslo,t2m)
#' 
#' # Plot the results for January month
#' # plot(subset(ds,it='Jan'))
#' 
#' # Plot the residuals:
#' residual <- as.residual(ds)
#' obs <- as.anomaly(as.calibrationdata(ds))
#' 
#' plot.zoo(obs[,1],lwd=2, new=FALSE)
#' lines(residual,col="red")
#' 
#' print("Global climate model simulation NorESM")
#' T2m <- t2m.NorESM.M(lon=lon,lat=lat)
#' 
#' # Plot the global mean of the field:
#' plot(T2m, new=FALSE)
#' # Plot area mean of a sub region
#' plot(T2m,is=list(lon=c(0,10),lat=c(60,70)), new=FALSE)
#' 
#' # Plot interpolated results corresponding to ferder
#' data(ferder)
#' plot(T2m,ferder, new=FALSE)
#' 
#' # Plot Hovmuller diagram: Not working ...
#' ## plot(T2m,is=list(lon=0)) 
#' 
#' print("Extract a subset - the January month")
#' x <- subset(t2m,it="jan")
#' X <- subset(T2m,it="jan")
#' 
#' print("Combine the fields for computing common EOFs:")
#' XX <- combine(x,X)
#' 
#' print("Compute common EOFs")
#' eofxx <- EOF(XX)
#' plot(eofxx, new=FALSE)
#' 
#' print("Downscale the January mean temperature") 
#' ds.jan <- DS(Oslo,eofxx)
#' plot(ds.jan, new=FALSE)
#'
#' @exportS3Method 
#' @export plot.cca
plot.cca <- function(x,...,ip=1,
                     colbar1=list(pal=NULL,rev=FALSE,n=10,breaks=NULL,type="p",cex=2,show=TRUE,
                                  h=0.6, v=1,pos=0.05),
                     what=c("maps","timeseries"),
                     colbar2=NULL,new=TRUE,verbose=FALSE) {
  if (verbose) print("plot.cca")
  if (new) dev.new()
  if (is.null(colbar2)) colbar2 <- colbar1
  par0 <- par()
  
  if("maps" %in% what) panels <- length(what)+1 else panels <- length(what)
  if(panels==3) {
    nf <- layout( matrix(c(1,2,3,3), nrow=2, byrow=TRUE))
  } else if(panels>1) {
    m <- seq(1,panels)
    if(panels %% 2) m <- c(m, max(m)+1)
    nf <- layout(matrix(m, nrow=round(max(m)/2), byrow=TRUE))
  }
  
  if("maps" %in% what) map(x,ip=ip,colbar1=colbar1,colbar2=colbar2,
                           verbose=verbose,new=FALSE)
  
  if("timeseries" %in% what) {
    w.m <- zoo((x$w.m[,ip]-mean(x$w.m[,ip],na.rm=TRUE))/
                 sd(x$w.m[,ip],na.rm=TRUE),order.by=x$index)
    v.m <- zoo((x$v.m[,ip]-mean(x$v.m[,ip],na.rm=TRUE))/
                 sd(x$v.m[,ip],na.rm=TRUE),order.by=x$index)
    r <- cor(x$w.m[,ip],x$v.m[,ip])
    plot(w.m,col="blue",lwd=2,new=FALSE,
         main=paste("CCA pattern ",ip," for ",varid(x),
                    "; r= ",round(r,2),sep=""),xlab="",ylab="")
    lines(v.m,col="red",lwd=2)
    legend("topleft",#0.01,0.90,
           c(paste(attr(x$X,'source')[1],attr(x$X,'variable')[1]),
             paste(attr(x$Y,'source')[1],attr(x$Y,'variable')[1])),
           col=c("red","blue"),lwd=2,lty=1,
           bty="n",cex=0.75,ncol=1,text.col="grey40")
  }
  
  ## KMP 2023-02-09: What is this fourth panel in plot.cca? 
  ## All I see here is a legend so I moved it into the time series above
  #par(xaxt="n",yaxt="n",bty="n",mar=c(0,0,0,0))
  #plot(c(0,1),c(0,1),type="n",xlab="",ylab="",new=FALSE)
  #legend(0.01,0.90,c(paste(attr(x$X,'source')[1],attr(x$X,'variable')[1]),
  #                   paste(attr(x$Y,'source')[1],attr(x$Y,'variable')[1])),
  #       col=c("red","blue"),lwd=2,lty=1,
  #       bty="n",cex=0.5,ncol=2,text.col="grey40")

  ## KMP 2023-02-22: reset if layout has been used
  if(panels>1) par(mfrow=c(1,1))
  par(mar=par0$mar, mgp=par0$mgp, xaxt=par0$xaxt , yaxt=par0$yaxt)
}

# Plot esd objects
# 
# These plot functions are S3 methods for esd objects, based on \code{plot}.
#
# @seealso plot plot.station plot.ds plot.eof plot.field plot.dsensemble
#
# @param x the object to be plotted
# @param \dots additional arguments
# 
# @return None
#
# @keywords plot graphics
#
# AM function plot.list is defined twice
# @exportS3Method plot list
# plot.list <- function(x,...) {
#  plot(combine.ds(x),...)
#}

#' plot cross-validation
#'
#' @param x input object to be shown in plot
#' @param new a boolean; if TRUE plot in new window
#' @param verbose a boolean; if TRUE print information about progress
#'
#' @export plot.xval
plot.xval <- function(x,...,ip=1,new=TRUE,verbose=FALSE) {
  if(verbose) print("plot.xval")
  if (new) dev.new()
  par0 <- par()
  par(bty="n")
  unit <- attr(x,'unit')
  cols <- rgb(seq(0,1,length=20),rep(0,20),rep(0,20))
  cexs <- seq(1.5,0.5,length=20)^2
  eofindex <- as.integer(gsub('X.','',
                              names(attr(x,'original_model')$coefficients[-1])))
  print(eofindex)
  if (unit=='deg C') unit <- expression(degree*C)
  
  class(x) <- "zoo"
  plot(x[,1],type="b",pch=19,lwd=2,
       main=paste("Cross-validation:",attr(x,'location'),attr(x,'variable')),
       xlab="",ylab= unit,
       sub=paste("r=",attr(x,'correlation'),"rmse=", attr(x,'rmse'),
                 unit))
  lines(x[,2],lwd=2,col="red")
  lines(attr(x,'fitted_values_all'),col="red",lty=3,pch="x")
  
  ## par(new=TRUE,fig=c(0.1,1,0.1,0.5),bty="n")
  ## plot(c(0,1),c(0,1),col="white",xlab="",ylab="",axes=F)
  legend(rep(range(index(x))[1],2), rep(range(x)[2],2),
         c("obs","x-valid","fit to all"), col=c("black","red","red"),
         lwd=c(2,2,1), lty=c(1,1,3), bty="n")
  
  dev.new()
  par(bty="n")
  boxplot(t(attr(x,'beta')[-1,]),col="grey90",border="grey40",
          xlim=range(eofindex),main="Model coefficients",
          ylab=expression(beta),xlab="EOF number")
  for (i in 1:20) 
    points(eofindex,attr(x,'original_model')$coefficients[-1],
           pch=19,col=cols[i],cex=cexs[i])
  grid()
}

#' plot dsensemble pca results
#'
#' @param x input object to be plotted
#' @param pts a boolean; if TRUE plot points?
#' @param target.show a boolean; if TRUE show diagnostics as a target (see \code{\link{diagnose}})
#' @param map.show a boolean; if TRUE show map of stations
#' @param legend.show a boolean; if TRUE show legend
#' @param verbose a boolean; if TRUE print information about progress
#' @param \dots additional arguments
#'
#' @exportS3Method 
#' @export plot.dsensemble.pca
plot.dsensemble.pca <- function(x,...,pts=FALSE,target.show=TRUE,map.show=TRUE,
                                it=0,ip=1,envcol=rgb(1,0,0,0.2),
                                legend.show=TRUE,verbose=FALSE) {
  if (verbose) print("plot.dsensemble.pca")
  stopifnot(inherits(x,'dsensemble') & inherits(x,'pca'))

  d <- index(x[[3]])
  pc <- x[3:length(x)]
  pc <- array(unlist(pc), dim = c(dim(pc[[1]]), length(pc)))
  pc <- lapply(seq(dim(pc)[2]), function(x) pc[ , x, ])
  fn <- function(x) {
    x <- zoo(x,order.by=d)
    class(x) <- c("dsensemble","station","zoo")
    invisible(x)
  }
  pc <- lapply(pc,fn)
  for (i in 1:length(pc)) {
    attr(pc[[i]],"station") <- as.station(x[[2]][,i],param=attr(x,"variable"),
                                          longname=attr(x,"longname"),unit=attr(x,"unit"))
  }
  plot(pc[[ip]],ylab=paste("PC",ip,sep=""))
}

#' Plot esd objects
#' 
#' These plot functions are S3 methods for esd objects, based on \code{plot}.
#'
#' \code{plot.dsensemble} plots a downscaled GCM ensemble.
#' 
#' @seealso \code{\link{plot}} plot.station plot.ds plot.eof
#'
#' @param x the object to be plotted
#' @param verbose a boolean; if TRUE print information about progress
#' @param plot a boolean; if TRUE show plot
#' @param \dots additional arguments
#' 
#' @return None
#'
#' @exportS3Method 
#' @export plot.dsensemble
plot.dsensemble <- function(x,verbose=FALSE,plot = TRUE, ...) {
  if(verbose) print("plot.dsensemble")
  if (inherits(x,c('pca','eof'))) {
    y <- plot.dsensemble.multi(x,verbose=verbose, plot = plot, ...) 
  } else if (inherits(x,'zoo')) {
    y <- plot.dsensemble.one(x,verbose=verbose,...) 
  } else if (inherits(x,'station')) {
    x <- as.station(x,verbose=verbose)
    y <- plot(x,verbose=verbose,...)
  } else {
    print(paste('Unknown class - do not know how to plot',
                paste(class(x),collapse=", ")))
    y <- x
  }
  if(verbose) print("exit plot.dsensemble")
  invisible(y)
}

#' Plot multiple stations/spatially aggregated field.
#'
#' @seealso \code{\link{plot}} plot.station plot.ds plot.eof
#'
#' @param x input object of class 'dsensemble'
#' @param it a time index, see \code{\link{subset}}
#' @param FUNX a function
#' @param verbose a boolean; if TRUE print information about progress
#' @param anomaly a boolean; if TRUE show anomalies
#' @param test a boolean; if TRUE perform some test?
#' @param plot a boolean; if TRUE show plot
#' @param \dots additional arguments
#'
#' @export plot.dsensemble.multi
plot.dsensemble.multi <- function(x,it=c(2000,2099),FUNX='mean',verbose=FALSE,
                                  anomaly=FALSE,test=FALSE, plot = TRUE, ...) {
  if (verbose) print('plot.dsensemble.multi')
  if (inherits(x,c('pca','eof'))) {
    Y <- expandpca(x,it=it,FUNX=FUNX,verbose=verbose,anomaly=anomaly,test=test)
    if (plot) plot(Y,verbose=verbose,...)
    #invisible(Y)
  } else {
    #return(NULL)
    Y <- NULL
  }
  if (verbose) print('exit plot.dsensemble.multi')
  invisible(Y)
}

#' Plot one station
#'
#' @seealso plot.dsensemble
#' 
#' @param x input object of class 'dsensemble'
#' @param pts a boolean; if TRUE show points
#' @param envcol color of envelope
#' @param obs.show if TRUE show observations
#' @param target.show if TRUE show diagnosis as target plot
#' @param map.type type of map, 'points' or 'rectangle'
#' @param map.insert if TRUE show map as insert, else show in new window
#' @param alpha transparancy factor of main plot
#' @param alpha.map transparancy factor of map
#' @param mar see \code{\link[graphics]{par}}
#' @param cex.axis see \code{\link[graphics]{par}}
#' @param cex.lab see \code{\link[graphics]{par}}
#' @param cex.main see \code{\link[graphics]{par}}
#' @param verbose a boolean; if TRUE print information about progress
#' @param anomaly a boolean; if TRUE show anomalies
#' @param test a boolean; if TRUE perform some test?
#' @param plot a boolean; if TRUE show plot
#' @param xval a boolean; if TRUE show cross-validation statistics
#' @param \dots additional arguments
#'
#' @export plot.dsensemble.one
plot.dsensemble.one <-  function(x,pts=FALSE,it=0,
                                 envcol=rgb(1,0,0,0.2),legend.show=FALSE,ylab=NULL,
                                 obs.show=TRUE,target.show=TRUE,map.show=TRUE,map.type=NULL,map.insert=TRUE,
                                 new=FALSE,xrange=NULL,yrange=NULL,
                                 alpha=0.5,alpha.map=0.7,mar=c(5.1,4.5,4.1,2.1),
                                 cex.axis=1, cex.lab=1.2, cex.main=1.2, xval.show=FALSE,
                                 verbose=FALSE,...) {
  if(verbose) print("plot.dsensemble.one")
  stopifnot(inherits(x,'dsensemble'))
  
  if (is.null(map.type)) {
    if (verbose) print(class(x))
    if( inherits(x,"field") | length(lon(x))!=length(lat(x)) |
        (length(lon(x))==2 & length(lat(x))==2) ) {
      map.type <- "rectangle"
    } else {
      map.type <- "points"
    }
  }
  
  if (verbose) {print(map.type); print(attr(x,'station'))}
  if (!is.null(attr(x,'station')) & !inherits(attr(x,'station'),c('annual','season'))) {
    z <- subset(x,it=it,verbose=verbose) 
  } else {
    z <- x
  }
  
  if (verbose) {
    print("diagnose")
    class(z)
  }
  diag <- diagnose(z,plot=FALSE,verbose=verbose)
  
  y <- attr(z,'station')
  attr(y,'standard.error') <- NULL
  if (verbose) print(paste('lon=',lon(y),'lat=',lat(y))) 
  
  d <- dim(z)
  index(y) <- year(y)
  
  if(map.show | target.show) {
    pscl <- c(0.9,1.3)
  } else {
    pscl <- c(0.9,1.1)
  }
  
  if (max(coredata(z),na.rm=TRUE) < 0) pscl <- rev(pscl)
  args <- list(...)
  if (verbose) print(names(args))
  ixl <- grep('xlim',names(args))
  if (length(ixl)==0) xlim <- range(year(z)) else
    xlim <- args[[ixl]]
  iyl <- grep('ylim',names(args))
  if (length(iyl)==0) ylim <- pscl*range(coredata(z),na.rm=TRUE) else
    ylim <- args[[iyl]]  
  index(y) <- year(y)
  if (obs.show) obscol <- 'black' else obscol='white'
  plot(y,type="b",pch=19,xlim=xlim,ylim=ylim,col=obscol,main='',
       cex.axis=cex.axis,cex.lab=cex.lab,mar=mar,
       ylab=ylab,map.show=FALSE,new=new, verbose=verbose)
  grid()
  par0 <- par()
  usr <- par()$usr; mar <- par()$mar; fig <- par()$fig
  setfig <- FALSE
  t <- index(z)
  
  if (pts) for (i in 1:d[2]) {
    points(year(t),coredata(z[,i]),pch=19,col="red",cex=0.3)
  }

  # Produce a transparent envelope
  nt <- length(index(z))
  t2 <- c(year(t),rev(year(t)))
  
  col <- rgb(rep(1,49),seq(0.95,0.1,length=49),seq(0.95,0.1,length=49),0.1)
  ## REB 2016-11-25
  if(is.null(alpha.map)) alpha.map <- alpha
  col.map <- adjustcolor(col,alpha.f=alpha.map)
  col <- adjustcolor(col,alpha.f=alpha)
  
  mu <- apply(coredata(z),1,mean,na.rm=TRUE)
  si <- apply(coredata(z),1,sd,na.rm=TRUE)
  for (ii in 1:49) {
    qp1 <- qnorm(1-ii/50,mean=coredata(mu),sd=coredata(si))
    qp2 <- qnorm(ii/50,mean=coredata(mu),sd=coredata(si))
    ci <- c(qp1,rev(qp2))
    polygon(t2[!is.na(ci)],ci[!is.na(ci)], col= envcol, border=NA)
  }
  q05 <- qnorm(0.05,mean=mu,sd=si)
  q95 <- qnorm(0.95,mean=mu,sd=si)
  
  lcol <- adjustcolor(envcol,offset=c(0.5,0.5,0.5,0.2))
  lines(zoo(mu,order.by=year(z)),lwd=3,col=lcol)
  lines(zoo(q05,order.by=year(z)),lty=2,col=lcol)  
  lines(zoo(q95,order.by=year(z)),lty=2,col=lcol)
  if (obs.show) lines(y,type="b",pch=19)
  if (!is.null(diag)) {
    index(diag$y) <- year(diag$y)
    outside <- diag$above | diag$below
    points(zoo(coredata(diag$y)[which(outside)],
               order.by=year(diag$y)[which(outside)]),col="grey")
  }
  title(main=toupper(loc(x)),cex.main=cex.main)
  if ((target.show) & (!is.null(diag))) {
    if (verbose) print('add target diagnostic')
    dx0 <- fig[2]-fig[1]
    dy0 <- fig[4]-fig[3]
    fig.target <- c(fig[1]+dx0*0.25, fig[1]+dy0*0.43,
                    fig[3]+dx0*0.72, fig[3]+dy0*0.9)
    #fig=c(0.12,0.30,0.72,0.90)
    par(fig=fig.target,
        new=TRUE, mar=c(0,0,0,0),xaxt="s",yaxt="n",bty="n",
        cex.main=0.75,xpd=NA,col.main="grey30")
    plot(diag,map.show=FALSE,new=FALSE,cex=0.75)
    setfig <- TRUE
  } 
  
  if(map.show & !map.insert) {
    vis.map(x,col="red",map.type,add.text=FALSE,map.insert=map.insert,
            cex.axis=cex.axis,cex=1.5,xrange=xrange,yrange=yrange,
            verbose=verbose,...)
    new <- TRUE
  }
  
  if (legend.show & !is.null(diag)) {
    dx0 <- fig[2]-fig[1]
    dy0 <- fig[4]-fig[3]
    fig.legend <- c(fig[1]+dx0*0.1, fig[1]+dy0*0.5,
                    fig[3]+dx0*0.2, fig[3]+dy0*0.25)
    par(fig=fig.legend,#c(0.1,0.5,0.2,0.25),
        new=TRUE,mar=c(0,0,0,0),xaxt="n",yaxt="n",bty="n",xpd=NA)
    plot(c(0,1),c(0,1),type="n",xlab="",ylab="")
    legend(0.05,0.90,c(paste("Past trend:",round(diag$deltaobs,2)),
                       paste(diag$robs,'%'),
                       paste(diag$outside,"observations"),
                       "p-value: "),
           bty="n",cex=0.7,text.col="grey40")
    legend(0.5,0.90,c(paste(levels(factor(attr(y,'unit')))[1],"/decade",sep=""),
                      "ensemble trends > obs.",
                      "outside ensemble 90% conf.int.",
                      paste(round(100*pbinom(diag$outside,size=diag$N,prob=0.1)),"%")),
           bty="n",cex=0.7,text.col="grey40")
    setfig <- TRUE
  }
  
  if (map.show & map.insert) {
    par(fig=par0$fig)
    vis.map(x,col="red",map.type=map.type,cex=1.5,
            cex.axis=cex.axis*0.65,add.text=FALSE,
            map.insert=map.insert,
            xrange=xrange,yrange=yrange,
            verbose=verbose,...)
  }
  if(xval.show) {
    ## REB 2021-05-15 - this part is unfinished and disabled.
    ## crossval()
    dx0 <- fig[2]-fig[1]
    dy0 <- fig[4]-fig[3]
    fig.xval <- c(fig[1]+dx*0.35, fig[2]+dx*0.53,
                  fig[3]+dx*0.85, fig[4]+dx*0.9)
    x$info <- NULL; x$eof <- NULL; x$pca <- NULL
    xval <- lapply(x,function(x) diag(cor(attr(x,'evaluation'))[seq(2,2*n,by=2),seq(1,2*n-1,by=2)]))
    for (i in 1:n) { 
      iy <- (i-1)*0.5
      par(fig=fig.xval + c(0,0,-iy,-iy),#c(0.32,0.50,0.85-iy,0.90-iy),
          new=TRUE, mar=c(0,0,1,0),xaxt="s",yaxt="n",bty="n",
          cex=0.5,xpd=NA,col.main="grey30")
      hist(unlist(lapply(xval,function(x) x[i])),col='grey',lwd=2,xlim=c(-1,1),
           main=paste('X-validation correlation for PCA',i),xlab='correlation')
    }
    setfig <- TRUE
  }
  ## Reset to the settings of the main plot so that additional things may be 
  ## added after the function has finished (except for some things that cannot be set)
  dontset <- c("cin","cra","csi","cxy","din","page")
  for(p in names(par0)[!names(par0) %in% dontset]) {
    if(p!="fig" | (p=="fig" & setfig)) {
      eval(parse(text=paste0("par(",p,"=par0$",p,")")))
    }
  }
  if(verbose) print("exit plot.dsensemble.one")
  invisible(z)
}

#' @export plot.xsection
plot.xsection <- function(x,...) {
  #print("plot.xsection")
  d <- attr(x,'dimensions')
  #print(d)
  X <- coredata(x)
  
  if (d[1]==1) {
    attr(X,'longitude') <- index(x)
    attr(X,'latitude') <- attr(x,'latitude')
    attr(X,'dimensions') <- attr(x,'dimensions')[c(3,2)]
    
  } else {
    attr(X,'longitude') <- attr(x,'longitude')
    attr(X,'latitude') <- index(x)
    attr(X,'dimensions') <- attr(x,'dimensions')[c(1,3)]
    X <- t(X)
  }
  attr(X,'variable') <- attr(x,'variable')
  
  attr(X,'unit') <- attr(x,'unit')
  attr(X,'source') <- attr(x,'source')
  
  # print(dim(X)); print(c(length(lon(X)),length(lat(X))))
  lonlatprojection(x=X,what="fill",geography=FALSE,...)
}



#' Plot esd objects
#' 
#' These plot functions are S3 methods for esd objects, based on \code{plot}.
#'
#' \code{plot.spell} plots wet/dry or cold/warm spells.
#' 
#' @seealso \code{\link{plot}} plot.station plot.ds plot.eof plot.field
#'
#' @param x the object to be plotted
#' @param xlim range of x-axis
#' @param ylim range of y-axis
#' @param \dots additional arguments
#' 
#' @return None
#'
#' @exportS3Method
#' @export plot.spell
plot.spell <- function(x,...,xlim=NULL,ylim=NULL,verbose=FALSE) {
  if(verbose) print("plot.spell")
  bar <- function(x,col) {
    rect(x[1],x[2],x[3],x[4],col=col,border=col)
  }
  t <- index(x)
  h <- coredata(x[,1]); l <- coredata(x[,2])
  ih <- is.finite(h); il <- is.finite(l)
  h <- h[ih]; th1 <- t[ih]
  l <- l[il]; tl1 <- t[il]
  th2 <- th1 + h; tl2 <- tl1 + l
  par(bty="n")
  
  col <- c('red','blue')
  runs <- c('hot','cold')
  spelltype <- 'hot and cold'
  if (sum(is.element(attr(x,'variable'),c('wet','dry')))>0) {
    col <- c('darkblue','brown')
    runs <- c('wet','dry')
    spelltype <- 'wet and dry' 
  }
  
  tunit <- attr(x,'threshold.unit')[1]
  for (i in 1:length(tunit)) {
    if ( (is.na(tunit[i]) | is.null(tunit[i])) ) tunit[i] <- " "
    if ((tunit[i]=='degree Celsius') | (tunit[i]=='deg C') | (tunit[i]=='degC'))
      tunit[i] <- 'degree*C'
  }
  
  
  plot(range(t),c(-1,1)*max(c(h,l),na.rm=TRUE),type="n",
       xlab="",ylab="Spell length",xlim=xlim,ylim=ylim,
       main=paste(attr(x,'location')[1],": ",spelltype[1],sep=""))
  leg <- try(eval(parse(text=paste("expression(paste(X > ",
                                   attr(x,'threshold'),"*",tunit,"))"))))
  if (inherits(leg,'try-error')) leg <- ''
  text(t[1],0.75*max(c(h,l),na.rm=TRUE),leg,srt=90,cex=0.7,col="grey")
  leg <- try(eval(parse(text=paste("expression(paste(X <= ",
                                   attr(x,'threshold'),"*",tunit,"))"))))
  if (inherits(leg,'try-error')) leg <- ''
  text(t[1],-0.75*max(c(h,l),na.rm=TRUE),leg,srt=90,cex=0.7,col="grey")
  lines(range(t),rep(0,2))
  apply(cbind(th1,rep(0,length(h)),th2,h),1,bar,col[1])
  apply(cbind(tl1,rep(0,length(l)),tl2,-l),1,bar,col[2])
  
}

#' Plot esd objects
#' 
#' These plot functions are S3 methods for esd objects, based on \code{plot}.
#'
#' @seealso \code{\link{plot}} plot.station plot.ds plot.eof plot.field
#'
#' @param x the object to be plotted
#' @param main main title
#' @param sub subtitle
#' @param verbose a boolean; if TRUE print information about progress
#' @param \dots additional arguments
#' 
#' @return None
#'
#' @exportS3Method
#' @export plot.ssa
plot.ssa <- function(x,...,main="SSA analysis",sub="",verbose=FALSE)  {
  if(verbose) print("plot.ssa")
  ssa <- x
  if ( (class(ssa)[1]!="SSA") ) stop("Need an 'SSA' object")
  nt <- ssa$nt
  dev.new()
  plot(ssa$d,main=main,sub=sub,ylab="Singular value",pch=20,col="grey50")
  points(ssa$d)
  grid()
  dev.new()
  par(mfcol=c(3,1))
  plot(ssa$v[,1],type="l",main=main,sub=sub,
       xlab="Time",ylab="SSA vector: mode 1",lwd=3,col="grey70")
  grid()
  plot(ssa$v[,2],type="l",main=main,sub=sub,
       xlab="Time",ylab="SSA vector: mode 2",lwd=3,col="grey70")
  grid()
  plot(ssa$v[,3],type="l",main=main,sub=sub,
       xlab="Time",ylab="SSA vector: mode 3",lwd=3,col="grey70")
  grid()
  dev.new()
  par(mfcol=c(3,1))
  if (class(ssa)[3] == "monthly.station.record") {
    yy <- sort(rep(ssa$x$yy,12)); yy <- yy[1:ssa$Nm]
    mm <- rep(1:12,nt); mm <- mm[1:ssa$Nm]
    plot(yy + (mm-0.5)/12, ssa$u[,1],type="l",main=main,sub=sub,
         xlab="Time",ylab="SSA loadings",lwd=3,col="grey70")
    grid()
    plot(yy + (mm-0.5)/12, ssa$u[,2],type="l",main=main,sub=sub,
         xlab="Time",ylab="SSA loadings",lwd=3,col="grey70")
    grid()
    plot(yy + (mm-0.5)/12, ssa$u[,3],type="l",main=main,sub=sub,
         xlab="Time",ylab="SSA loadings",lwd=3,col="grey70")
    grid()
  } else if (class(ssa)[3] == "daily.station.record") {
    plot(ssa$x$yy[1:ssa$Nm] + ssa$x$mm[1:ssa$Nm]/12 + ssa$x$dd[1:ssa$Nm]/365,
         ssa$u[,1],type="l",main=main,sub=sub,
         xlab="Time",ylab="SSA loadings",lwd=3,col="grey70")
    grid()
    plot(ssa$x$yy[1:ssa$Nm] + ssa$x$mm[1:ssa$Nm]/12 + ssa$x$dd[1:ssa$Nm]/365,
         ssa$u[,2],type="l",main=main,sub=sub,
         xlab="Time",ylab="SSA loadings",lwd=3,col="grey70")
    grid()
    plot(ssa$x$yy[1:ssa$Nm] + ssa$x$mm[1:ssa$Nm]/12 + ssa$x$dd[1:ssa$Nm]/365,
         ssa$u[,3],type="l",main=main,sub=sub,
         xlab="Time",ylab="SSA loadings",lwd=3,col="grey70")
    grid()
  } else {
    plot(ssa$x$yy[1:ssa$Nm] + ssa$x$mm[1:ssa$Nm]/12 + ssa$x$dd[1:ssa$Nm]/365,
         ssa$u[,1],type="l",main=main,sub=sub,
         xlab="Time",ylab="SSA loadings",lwd=3,col="grey70")
    grid()
    plot(ssa$x$yy[1:ssa$Nm] + ssa$x$mm[1:ssa$Nm]/12 + ssa$x$dd[1:ssa$Nm]/365,
         ssa$u[,2],type="l",main=main,sub=sub,
         xlab="Time",ylab="SSA loadings",lwd=3,col="grey70")
    grid()
    plot(ssa$x$yy[1:ssa$Nm] + ssa$x$mm[1:ssa$Nm]/12 + ssa$x$dd[1:ssa$Nm]/365,
         ssa$u[,3],type="l",main=main,sub=sub,
         xlab="Time",ylab="SSA loadings",lwd=3,col="grey70")
    grid()
  }
}

#' Plot esd objects
#' 
#' These plot functions are S3 methods for esd objects, based on \code{plot}.
#'
#' @seealso \code{\link{plot}}
#'
#' @param x the object to be plotted
#' @param main main title
#' @param xlab label of x-axis
#' @param ylab label of y-axis
#' @param verbose a boolean; if TRUE print information about progress
#' @param col Colour see \code{\link[graphics]{par}}
#' @param \dots additional arguments
#' 
#' @return None
#'
#' @exportS3Method 
#' @export plot.nevents
plot.nevents <- function(x,verbose=FALSE,main=NULL,xlab=NULL,ylab=NULL,col=NULL,...) {
  # Plot the results from 
  if (verbose) print('plot.nevents')
  par(bty='n')
  if (is.null(main)) main <- loc(x)
  if (is.null(xlab)) xlab <- ""
  if (is.null(ylab)) ylab <- attr(x,'info')
  if (is.null(col)) {
    if (is.T(attr(x,'observation')))
      col <- c(rgb(0.5,0.5,0.7,0.5),rgb(0.8,0.5,0.5,0.5),rgb(0.8,0.5,0.8,0.5),
               rgb(0.3,0.3,0.6,0.5),rgb(0.6,0.3,0.3,0.5),rgb(0.6,0.3,0.6,0.5)) else
                 col <- c(rgb(0.3,0.3,0.6,0.5),rgb(0.8,0.5,0.5,0.5),rgb(0.5,0.5,0.7,0.5),
                          rgb(0.6,0.3,0.6,0.5),rgb(0.6,0.3,0.3,0.5),rgb(0.8,0.5,0.8,0.5))
  }
  plot.zoo(x,plot.type='single',lwd=5,main=main,
           xlab=xlab,ylab=ylab,col=col,...)
  grid()
  points(attr(x,'observation'),pch=19)
  lines(attr(x,'nwd.pre'),col=rgb(0.5,0.5,0.5,0.5))
}

#' Plot esd objects
#' 
#' These plot functions are S3 methods for esd objects, based on \code{plot}.
#'
#' \code{plot.trajectory} plots the number of events per year.
#' 
#' @seealso \code{\link{plot}}
#'
#' @param x the object to be plotted
#' @param main main title
#' @param col colour
#' @param xlab label of x-axis
#' @param ylab label of y-axis
#' @param legend.show if TRUE show legendp
#' @param is For subsetting in space - See \code{\link{subset}}. Can also be
#' a station value and if provided, the plotting will involve an interpolation
#' to the same coordinates as defined by \code{is}.
#' @param verbose a boolean; if TRUE print information about progress
#' @param col Colour see \code{\link[graphics]{par}}
#' @param lwd width of line
#' @param xlim range of x-axis
#' @param ylim range of y-axis
#' @param \dots additional arguments
#' 
#' @return None
#'
#' @examples
#' data(imilast.M03)
#' plot(imilast.M03, new=FALSE)
#'
#' @exportS3Method
#' @export plot.trajectory
plot.trajectory <- function(x,it=NULL,is=NULL,
                            main=NULL,xlim=NULL,ylim=NULL,
                            col=NULL,pch=0,type='l',lwd=3,
                            xlab="",ylab=NULL,new=TRUE,verbose=FALSE,...) {
  if(verbose) print("plot.trajectory")
  y <- subset(x,it=it,is=is)
  n <- count.trajectory(y,by='year')
  if(new) dev.new()
  plot.station(n,main=main,new=new,col=col,
               xlim=xlim,ylim=ylim,type=type,
               lwd=lwd,pch=pch,xlab=xlab,ylab=ylab,
               legend.show=FALSE,...)
  if (is.null(col)) col <- rainbow(length(y[1,]))
  if (!is.null(attr(y,"longitude")) |
      !is.null(attr(y,"latitude"))) {
    if(verbose) print("adding legend")
    par(xpd=TRUE)
    leg <- ""
    if (!any(is.na(attr(y,"longitude"))) & !is.null(attr(y,"longitude"))) {
      leg <- paste(leg,paste(round(range(attr(y,"longitude")),2),sep="-"),"E/",sep="")
    }
    if (!any(is.na(attr(y,"latitude"))) & !is.null(attr(y,"latitude"))) {
      leg <- paste(leg,paste(round(range(attr(y,"latitude")),2),sep="-"),"N",sep="")
    }
    if (is.null(col)) col <- rainbow(1)
    if(verbose) print(paste('legend:',leg,', color:',col[1]))
    legend("bottomleft",inset=c(0,-0.25),legend=leg,bty="n",cex=0.6,ncol=3,
           text.col="grey40",lty=1,col=col)
  }
  invisible(n)
}

#' name to expression - only valid for temperature
#'
#' @export
nam2expr <- function(x) {
  y <- x
  for (i in 1:length(y)) {
    z <- switch(tolower(x[i]),
                't2m'=expression(T[2 * m]),
                'tmax'=expression(T[x]),
                'tmin'=expression(T[n]),
                'tas'=expression(T[2 * m]))
    if (is.null(z)) y[i] <- x[i] else y[i] <- z
  }
  return(y)
}

## REB 2021-11-26: a plot routine for DS-objects that have used PCA (stations) as predictors
plot.ds.station.pca <- function(x,...,plot.type="multiple",what=NULL,new=TRUE,
                                lwd=1,type='l',pch=0,main=NULL,col=NULL,
                                colbar=list(pal=NULL,rev=FALSE,n=10,
                                            breaks=NULL,type="p",cex=2,show=TRUE,
                                            h=0.6, v=1,pos=0.05),
                                xlim=NULL,ylim=NULL,xlab="",ylab=NULL,verbose=FALSE) {
  if(is.null(what)) what <- c("map","ts",'xval')
  if (verbose) print(paste('plot.ds.station.pca',paste(what,collapse=',')))
  unit <- attr(x,'unit')
  if ( (is.na(unit) | is.null(unit)) ) unit <- " "
  for (i in 1:length(unit)) {
    if ((unit[i]=='degree Celsius') | (unit[i]=='deg C') | (unit[i]=='degC'))
      unit[i] <- 'degree*C'
  }
  
  if (is.null(ylab))
    ylab <- esd::ylab(x)
  
  if (verbose)  print(ylab)
  if (is.null(main)) main <- attr(x,'longname')[1]               
  if (is.null(col)) col <- rainbow(length(x[1,]))  
  
  cols <- rep("blue",100)
  model <- attr(x,'model')
  
  #mfrow <- c(2,2)
  #if (length(what)==2) mfrow <- c(2,1) else
  #  if (length(what)==1) mfrow <- c(1,1) 
  
  if (new) dev.new()
  if (plot.type=="single") new <- TRUE
  par(cex.axis=0.75,cex.lab=0.7,cex.main=0.8)
  
  if (sum(is.element(what,'map'))>0) {
    #par(fig=c(0,0.5,0.5,1))
    pattern <- attr(x,'pattern')
    if (!inherits(x,'radiosonde')) {
      if (verbose) print('Show map...')
      class(pattern) <- class(attr(x,'original_data'))
      attr(pattern,'longitude') <- lon(attr(x,'eof'))
      attr(pattern,'latitude') <- lat(attr(x,'eof'))
      map(pattern,FUN='mean',new=FALSE,colbar=list(show=FALSE),
          fig=c(0,0.5,0.5,1),verbose=verbose,...)
      points(lon(x),lat(x),lwd=3,cex=1.5)
    } else {
      pattern <- attr(pattern,'pattern') ## Here is a quick fix for a bug...
      if (verbose) print('Show vertical profiles (radiosonde)')
      if (verbose) str(pattern)
      par(fig=c(0,0.5,0.5,1))
      plot(pattern[,1],alt(attr(x,'eof')),type='l',col='grey',xlab='')
      np <- dim(pattern)[1]
      for (i in 2:np) lines(pattern[,i],alt(attr(x,'eof')),col='grey')
      lines(pattern[,1],alt(attr(x,'eof')))
    }
  }
  
  if ( (sum(is.element(what,'xval'))>0)  & (!is.null(attr(x,'evaluation'))) ) {
    par(new=TRUE,fig=c(0.5,1,0.5,1)) 
    plot(attr(x,'evaluation')[,1],attr(x,'evaluation')[,2],
         main='Cross-validation',xlab='original data',
         ylab='prediction',pch=19,col="grey")
    lines(range(c(attr(x,'evaluation')),na.rm=TRUE),
          range(c(attr(x,'evaluation')),na.rm=TRUE),lty=2)
    cal <- data.frame(y=coredata(attr(x,'evaluation')[,1]),
                      x=coredata(attr(x,'evaluation')[,2]))
    xvalfit <- lm(y ~ x, data = cal)
    abline(xvalfit,col=rgb(1,0,0,0.3),lwd=2)
    #par(bty="n",fig=c(0.6,0.95,0.48,0.52),mar=c(0,0,0,0),new=TRUE,
    #    xaxt='n',yaxt='n',cex.sub=0.7)
    #plot(c(0,1),c(0,1),type='n',xlab='',ylab='')
    ok <- is.finite(attr(x,'evaluation')[,1]) &
      is.finite(attr(x,'evaluation')[,2])
    text(par()$usr[1] + diff(range(par()$usr[1:2]))/24,
         par()$usr[4] - diff(range(par()$usr[3:4]))/12,
         paste('correlation=',
               round(cor(attr(x,'evaluation')[ok,1],attr(x,'evaluation')[ok,2]),2)),
         pos=4,cex=0.8,col='grey')
  }  else {
    xvalfit <- NULL
  }
  
  Y0 <- as.original.data(x)
  #print(index(Y0)); print(index(x))
  yX <- merge.zoo(Y0,x,all=FALSE)
  #print(summary(yX))
  y0 <- yX$Y0
  ## KMP 2021-02-26: Plot attribute 'fitted value' instead of coredata
  
  if (!is.null(attr(x,'n.apps'))) ns <- attr(x,'n.apps') else
    ns <- 0
  y.rng <- NA; x.rng <- NA
  if (ns > 0) {
    #print("Add other DS results")
    for (i in 1:ns) {
      eval(parse(text=paste("y <- attr(x,'appendix.",i,"')",sep="")))
      #print(summary(y))
      y.rng <- range(y.rng,y,na.rm=TRUE)
      x.rng <- range(x.rng,index(y),na.rm=TRUE)
    }
    ## KMP 2021-04-26: added following to solve problem when index is Date
    if(is.numeric(x.rng) & is.dates(index(x))) x.rng <- as.Date(x.rng)
  }
  
  if (is.null(ylim)) {
    #ylim <- range(coredata(x),coredata(y0),y.rng,na.rm=TRUE)
    ylim <- range(attr(x,"fitted_values"),coredata(y0),y.rng,na.rm=TRUE)
  }
  if (is.null(xlim)) {
    xlim <- range(index(x),index(y0),x.rng,na.rm=TRUE)
  }
  
  par(fig=c(0.025,1,0.025,0.475),new=TRUE)
  par(bty="n",fig=c(0,1,0.1,0.5),mar=c(1,4.5,1,1),new=TRUE, xaxt='s',yaxt='s')
  ds <- list(obs=y0)
  if (verbose) {print('zoo-plot'); print(main); print(xlab); print(ylab)}
  plot.zoo(y0,plot.type=plot.type,ylab=ylab,xlab=xlab,
           main=main,xlim=xlim,ylim=ylim,lwd=1,type='b',pch=19)
  par0 <- par()
  grid()
  if (verbose) print(c(class(index(x)),class(index(y0))))
  if ( (class(index(x))=='Date') & (class(index(y0))=='numeric') & inherits(x,'annual') ) 
    index(x) <- year(index(x))
  if ( (class(index(x))=='numeric') & (class(index(y0))=='Date') & inherits(x,'annual') ) 
    index(x) <- as.Date(paste(index(x),'01-01',sep='-'))
  #lines(x,col="red",type="l",lwd=lwd)
  lines(attr(x,"fitted_values"),col="red",type="l",lwd=lwd)
  
  cal0 <- data.frame(y=coredata(y0),t=year(y0))
  #cal1 <- data.frame(y=coredata(x),t=year(x))
  cal1 <- data.frame(y=attr(x,"fitted_values"),t=year(x))
  
  trend0 <- lm(y ~ t, data=cal0)
  trend1 <- lm(y ~ t, data=cal1)
  lines(zoo(predict(trend0),order.by=index(y0)),lty=2)
  lines(zoo(predict(trend1),order.by=index(x)),lty=2,col='red')
  
  st0 <- summary(trend0); st1 <- summary(trend1)
  obstrend <- paste('obs. trend: ', round(st0$coefficients[2],2),' (',
                    round(st0$coefficients[2]-2*st0$coefficients[4],2),',',
                    round(st0$coefficients[2]+2*st0$coefficients[4],2),')',
                    attr(x,'unit'),'/decade',sep='')
  dstrend <- paste('obs. trend: ', round(st1$coefficients[2],2),' (',
                   round(st1$coefficients[2]-2*st1$coefficients[4],2),',',
                   round(st1$coefficients[2]+2*st1$coefficients[4],2),')',
                   attr(x,'unit'),'/decade',sep='')
  
  if (is.null(attr(x,'source'))) attr(x,'source') <- 'ESD'
  if (is.na(attr(x,'source'))) attr(x,'source') <- 'ESD'
  legtext <- c("Observations",attr(x,'source')) 
  legcol <- c("black","red")
  
  if (ns > 0) {
    #print("Add other DS results")
    for (i in 1:ns) {
      eval(parse(text=paste("y <- attr(x,'appendix.",i,"')",sep="")))
      lines(zoo(coredata(y),order.by=index(y)),col=cols[i],lwd=lwd)
      legcol <- c(legcol,cols[i])
      legtext <- c(legtext,attr(y,'source'))
      eval(parse(text=paste("ds$result.",i,
                            " <- attr(x,'appendix.",i,"')",sep="")))
    }
    
  } else  legcol <- c("black","red")
  
  ## Replot observations and prediction for calibration period
  
  lines(y0,lwd=1,type='b',pch=19)
  #lines(x,col="red",type="l",lwd=lwd)
  lines(attr(x,"fitted_values"),col="red",type="l",lwd=lwd)
  #print(legcol)
  if (!is.null(attr(x,'appendix.1'))) legend <- c("Obs.","Cal.","Proj") else
    legend <- c("Obs.","Cal.")
  legend(x="topleft",legend=legend,bty="n",horiz=TRUE,
         col=c("black","red","blue"),lwd=c(1,1,1),pch=c(19,1,1))
  
  if (plot.type=="single") {
    par(fig=c(0,1,0,0.1),new=TRUE, mar=c(0,0,0,0),xaxt="s",yaxt="s",bty="n")
    plot(c(0,1),c(0,1),type="n",xlab="",ylab="")
    legend(0.01,0.95,c(paste(attr(x,'location'),": ",
                             #attr(x,'aspect'),
                             #attr(x,'longname')," - ",
                             round(attr(x,'longitude'),2),"E/",
                             round(attr(x,'latitude'),2),"N (",
                             attr(x,'altitude')," masl)",sep=""),
                       obstrend,dstrend),ncol=3,
           bty="n",cex=0.6,text.col="grey40",
           lwd=c(1,rep(2,ns),1,1),lty=c(1,rep(1,ns),2,2),
           col=c(col,'black',col[1]))
    
    par(fig=c(0,1,0.05,0.95),new=TRUE,mar=par0$mar,xaxt="n",yaxt="n",bty="n")
    #plot.zoo(x,plot.type=plot.type,type="n",ylab="",xlab="",xlim=xlim,ylim=ylim)
    plot.zoo(attr(x,"fitted_values"),plot.type=plot.type,type="n",
             ylab="",xlab="",xlim=xlim,ylim=ylim)
  }
  invisible(list(trend0=trend0,trend1=trend1,xvalfit=xvalfit))
}

#' @export plot.radiosonde
plot.radiosonde <- function(x,...,plot.type="single",new=TRUE,
                            lwd=3,type='l',pch=0,main=NULL,col=NULL,
                            xlim=NULL,ylim=NULL,xlab="",ylab=NULL,
                            errorbar=TRUE,legend.show=FALSE,
                            map.show=TRUE,map.type=NULL,map.insert=TRUE,
                            cex.axis=1.2,cex.lab=1.2,cex.main=1.2,
                            mar=c(4.5,4.5,0.75,0.5),fig=NULL, 
                            alpha=0.5,alpha.map=0.7,add=FALSE,
                            verbose=FALSE) {
  d <- dim(x)
  plot(zoo(x),plot.type='single',col=heat.colors(d[2]),xlab='',ylab=esd::unit(x),main=loc(x))
}

