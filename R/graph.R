#' InfoGraphics
#'
#' Various functions for visual display of data and statistics
#'
#' \code{graph} shows a fancy graph of output of \code{\link{DSensemble}}.
#'
#' @aliases graph.default graph.dsensemble graph.list graph.zoo
#' @seealso wheel cumugram visprob conf vis diagram scatter plot map
#' 
#' @param x an input object of class 'DSensemble'
#' @param img a 'raster' object, or an object that can be coerced to one by 'as.raster', to be used as background
#' @param pch see \code{\link{par}}
#' @param it see \code{\link{subset}}
#' @param col color
#' @param lwd line width
#' @param xlim range of x-axis
#' @param ylim range of y-axis
#' @param ensmean a boolean; if TRUE plot the ensemble mean, if FALSE show all members
#' @param new a boolean; if TRUE open new graphic device
#' @param col.obs color of markers representing observations
#' @param verbose a boolean; if TRUE print information about progress
#' @param \dots additional arguments
#'
#' @examples
#' data(dse.Oslo)
#' graph(dse.Oslo)
#'
#' @export
#' @export
graph <- function(x,...) UseMethod("graph")

#' @export
graph.default <- function(x,...,img=NULL,pch='fancy',it=NULL,col=rgb(0.5,0.5,0.5,0.5),lwd=5,
                          xlim=NULL,ylim=NULL,new=TRUE,col.obs='black',verbose=FALSE) {
    if(verbose) print('graph.default')
    ## Produce the graphics:
    
    if (new) dev.new()
    if (!is.null(img)) {
      par0 <- par()
      par(mar=rep(0,4))
      plot(c(0,1),c(0,1),type='n')
      rasterImage(img, -0.05, -0.05, 1.05, 1.05)
      par(new=TRUE,col.axis='white',col.lab='white',xaxt='n',yaxt='n',
          mar=par0$mar,bty='n',col.sub='white')
    }
    if (!is.null(it)) y <- subset(x,it=it) else y <- x
    plot.zoo(y,lwd=lwd,col=col,ylim=ylim,xlim=xlim,
             ylab=ylab(y),sub=loc(y))
    #balls(y)
    if (!is.null(pch)) if (pch=='fancy') balls(y,col=col.obs) else points(zoo(y),pch=pch,col=col.obs)
    par(xaxt='s',yaxt='s')
    axis(1,col='white')
    axis(2,col='white')
}

#' @export
graph.dsensemble <- function(x,img=NULL,pch='fancy',it=NULL,col=rgb(1,0.7,0.7,0.1),
                             lwd=5,xlim=NULL,ylim=NULL,add=FALSE,new=TRUE,ensmean=FALSE,
			     col.obs='black',verbose=FALSE) {
    if(verbose) print('graph.dsensemble')
    ## Produce the graphics:
    if ((!add) & (new)) dev.new()
    if (!is.null(img)) {
      par0 <- par()
      par(mar=rep(0,4))
      plot(c(0,1),c(0,1),type='n')
      rasterImage(img, -0.05, -0.05, 1.05, 1.05)
      par(new=TRUE,col.axis='white',col.lab='white',xaxt='n',yaxt='n',
          mar=par0$mar,bty='n',col.sub='white')
    }
    if (!is.null(it)) y <- subset(x,it=it) else y <- x
    index(y) <- year(y)
    index(attr(y,'station')) <- year(attr(y,'station'))
    if (is.null(xlim)) xlim <- range(index(y))
    if (is.null(ylim)) ylim <- range(coredata(y),na.rm=TRUE)
    
    if (!add) plot.zoo(attr(y,'station'),lwd=lwd,col=rgb(0.5,0.5,0.5,0.5),
                       ylim=ylim,xlim=xlim,ylab=ylab(attr(y,'station')),
                       sub=loc(x),plot.type='single',xlab='')
    for (i in 1:dim(x)[2]) lines(y[,i],lwd=7,col=col)
    if (ensmean) lines(index(x),apply(coredata(x),1,'mean',na.rm=TRUE),
                                 lwd=3,col='red')

    #balls(attr(y,'station'))
    if (!is.null(pch)) if (pch=='fancy') balls(attr(y,'station'),col=col.obs) else points(zoo(attr(y,'station')),pch=pch,col=col.obs)
    par(xaxt='s',yaxt='s')
    if (!is.null(img)) col.axis <- 'white' else col.axis <- 'black'
    axis(1,col=col.axis)
    axis(2,col=col.axis)
}

#' @export
graph.list <- function(x,img=NULL,pch='fancy',it=NULL,
                       col=c(rgb(1,1,0.5,0.05),rgb(1,0.5,0.5,0.05),rgb(0.5,1,0.5,0.05),
                             rgb(0.5,0.5,0.5,0.05) ),
                       lwd=5,xlim=NULL,ylim=NULL,add=FALSE,new=TRUE,ensmean=FALSE,
		       col.obs='black', verbose=FALSE) {
  if(verbose) print("graph.list")		       
  if ((!is.null(it)) & (inherits(x[[1]],'dsensemble')))
    y <- subset(x[[1]],it=it) else y <- x[[1]]
    index(y) <- year(y)
  graph(y,img=img,pch=pch,col=col[1],lwd=lwd,xlim=xlim,ylim=ylim,add=add,new=new,col.obs=col.obs)
  if (!is.null(attr(x,'obs')) & is.null(attr(y,'dsensemble'))) obs <- attr(x,'obs') else
                                                               obs <- attr(y,'station')
  index(obs) <- year(obs)
  
  for (j in c(2:length(x),1)) {
    if ((!is.null(it)) & (inherits(x[[j]],'dsensemble')))
         y <- subset(x[[j]],it=it) else y <- x[[j]]
         index(y) <- year(y)
    for (i in 1:dim(y)[2]) lines(y[,i],lwd=7,col=col[j])
  }
  lines(obs,lwd=3,col=rgb(0.5,0.5,0.5,0.25))

  if (ensmean) {
    emcol <- c('wheat','red','green','grey')    
    for (i in 1:length(x)) lines(index(x[[i]]),apply(coredata(x[[i]]),1,'mean',na.rm=TRUE),
                                 lwd=3,col=emcol[i])
    legend(index(y)[1],max(coredata(y)),names(x),col=emcol,lty=1,lwd=3,bty='n')
  }
  if (!is.null(pch)) if (pch=='fancy') balls(obs,col=col.obs) else points(obs,pch=pch,col=col.obs)
}

#' @export
graph.zoo <- function(x,img=NULL,it=NULL,col=rgb(1,0.7,0.7,0.1),pch=1,
                      lwd=5,xlim=NULL,ylim=NULL,xlab='',ylab='',add=FALSE,
                      new=TRUE,ensmean=FALSE,col.obs='black', verbose=FALSE) {
    if(verbose) print('graph.zoo')
    ## Produce the graphics:
    if ((!add) & (new)) dev.new()
    if (!is.null(img)) {
      par0 <- par()
      par(mar=rep(0,4))
      plot(c(0,1),c(0,1),type='n')
      rasterImage(img, -0.05, -0.05, 1.05, 1.05)
      par(new=TRUE,col.axis='white',col.lab='white',xaxt='n',yaxt='n',
          mar=par0$mar,bty='n',col.sub='white')
    }
    if (!is.null(it)) y <- subset(x,it=it) else y <- x
    index(y) <- year(y)
    if (is.null(xlim)) xlim <- range(index(y))
    if (is.null(ylim)) ylim <- range(coredata(y),na.rm=TRUE)
    if (!add) plot.zoo(y[,1],lwd=lwd,col=col,ylim=ylim,xlim=xlim,
                       ylab=ylab,xlab=xlab,plot.type='single')
    grid()
    if(is.null(dim(x))) {
      lines(y,lwd=7,col=col)
    } else {
      for (i in 1:dim(x)[2]) lines(y[,i],lwd=7,col=col)
    }

    if (!is.null(pch)) {
      if (pch=='fancy') {
        balls(attr(y,'station'),col=col.obs)
      } else {
        points(zoo(attr(y,'station')),pch=pch,col=col.obs)
      }
    }
    par(xaxt='s',yaxt='s')
    if (!is.null(img)) col.axis <- 'white' else col.axis <- 'black'
    axis(1,col=col.axis)
    axis(2,col=col.axis)
}
