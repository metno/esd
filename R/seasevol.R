#' Visualise the seasonal evolution of a daily time series
#'
#' Visualise the daily values of time series as a color scale on a plot
#' with the julian day on the y-axis and year on the x-axis. 
#'
#' @param x as \code{station} object with daily data
#' @param nv number of steps in color scale
#' @param verbose a boolean; if TRUE print information about progress
#'
#' @aliases seasevol seasevol.station
#'
#' @examples
#' data(ferder)
#' seasevol(ferder, new=FALSE)
#'
#' @export
seasevol <- function(x,nv=25,verbose=FALSE,...) UseMethod("seasevol")

#' @export seasevol.station
seasevol.station <- function(x, nv=25, new=TRUE, verbose=FALSE...) {

  stopifnot(inherits(x,'day'))
  yrs <- as.numeric(rownames(table(year(x))))
  #print(yrs)
  ny <- length(yrs)

  xn <- min(coredata(x),na.rm=TRUE)
  xx <- max(coredata(x),na.rm=TRUE)
  xi <- seq(floor(xn),ceiling(xx),length=nv)
  
  j <- 1:nv
  if (attr(x,'variable')=='t2m') {
    col <- rgb(j/nv,abs(sin(pi*j/nv)),(1-j/nv))
  } else if (attr(x,'variable')=='precip') {
    col <- rgb(1-j/nv,1-j/nv,1)
  } else {
    col <- rainbow(nv)
  }
  class(x) <- "zoo"

  Z <- matrix(rep(NA,ny*366),ny,366)
  
  if ( (attr(x,'unit') == "deg C") | (attr(x,'unit') == "degree Celsius") )
      unit <- expression(degree*C) else
      unit <- attr(x,'unit')
  
  for (i in 1:ny) {
    y <- window(x,start=as.Date(paste(yrs[i],'-01-01',sep='')),
                    end=as.Date(paste(yrs[i],'-12-31',sep='')))
    Z[i,1:length(y)] <- y
  }
  
  if(new) dev.new()
  par(bty="n",yaxt="n",fig=c(0.0,0.90,0.0,1.0))
  z <- coredata(x)
  image(yrs,1:366,Z,
       ylab="julian day",xlab="",
       col=col,
       main=paste(attr(x,'location'),attr(x,'variable'),attr(x,'unit')),
       sub=attr(x,'location'))
  grid()

  par(new=TRUE,fig=c(0.85,0.95,0.70,0.85),mar=c(0,3,0,0),
      cex.axis=0.7,yaxt="s",xaxt="n",las=1)
  colbar <- rbind(1:nv,1:nv)
  image(1:2,xi,colbar,col=col)

  par(new=TRUE,xaxt="n",yaxt="n",fig=c(0.0,0.90,0.0,1.0))
  plot(range(yrs),c(1,633),type="n")
  rownames(Z) <- yrs
  invisible(Z)
}
