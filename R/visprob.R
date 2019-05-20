#' InfoGraphics
#'
#' Various functions for visual display of data and statistics
#'
#' \code{visprob} displays the probability density function (PDF) of a precipitation time series (\code{x}) for each year.
#' If only one time series is provided (\code{y}=NULL), the color of the PDFs represent the year.
#' If a second time series is provided, the color scale shows the annual mean value of \code{y}.
#'
#' @aliases visprob.station visprob.station.precip
#' @seealso wheel cumugram climvar graph conf vis diagram scatter plot map
#' 
#' @param x an input object of class 'station'
#' @param y an input object of class 'station'
#' @param dy relative width of lines
#' @param threshold threshold defining a precipitation event 
#' @param breaks breaks in historgram 
#' @param pdf a boolean; if TRUE add pdfs estimated from wet-day mean
#' @param verbose a boolean; if TRUE print information about progress
#' @param \dots additional arguments
#'
#' @example
#' data(bjornholt)
#' visprob(bjornholt)
#'
#' @export
visprob <- function(x,...) UseMethod("visprob")

#' @export
visprob.default <- function(x,...) {

}

#' @export
visprob.station <- function(x,...,y=NULL,dy=0.01,verbose=FALSE) {
  if (is.precip(x)) visprob.station.precip(x,y=y,dy=dy,verbose=verbose,...) 
}

#' @export
visprob.station.precip <- function(x,...,y=NULL,threshold=1,dy=0.005,
                                   breaks=NULL,pdf=FALSE,verbose=FALSE) {
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

#' @export
histwet <- function(x,breaks=NULL,threshold=1) {
  if (is.null(breaks)) breaks=seq(0,1.1*max(x,na.rm=TRUE),by=5)
  h <- hist(x[x > threshold],breaks=breaks,plot=FALSE)
  return(h$density)
}
