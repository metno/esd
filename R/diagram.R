#' Visualise - different type of plotting...
#'
#' @param x input object 
#' @param \dots additional arguments
#'
#' @aliases diagram.dsensemble diagram.station
#'
#' @export
diagram <- function(x,...) UseMethod("diagram")

#' @exportS3Method
#' @export diagram.dsensemble
diagram.dsensemble <- function(x,...,it=0,verbose=FALSE) {
  if(verbose) print("diagram.dsensemble")
  stopifnot(inherits(x,'dsensemble'))
  if (!inherits(attr(x,'station'),'annual')) {
    z <- subset(x,it=it) 
  } else {
    z <- x
  }
  y <- attr(z,'station')
  pscl <- c(0.9,1.1)
  if (max(coredata(z),na.rm=TRUE) < 0) pscl <- rev(pscl)
  #print("...")
  plot(y,type="b",pch=19,
       xlim=range(year(z)),
       ylim=pscl*range(coredata(z),na.rm=TRUE))
  grid()
  #usr <- par()$usr; mar <- par()$mar; fig <- par()$fig
  t <- year(z); n <- dim(z)[2]
  col <- rgb(seq(1,0,length=n)^2,sin(seq(0,pi,length=n))^2,seq(0,1,length=n)^2,0.2)
  for (i in 1:n) lines(t,z[,i],col=col[i],lwd=2)
  points(y,pch=19,lty=1)
}

# not exported
diagram.ds <- function(x,...,verbose=FALSE) {
  if(verbose) print("diagram.ds - unfinished function")
}

# Show the temperatures against the day of the year. Use
# different colours for different year.
#' @exportS3Method
#' @export diagram.station
diagram.station <- function(x,...,it=NULL,new=TRUE,plot=TRUE,verbose=FALSE) {
  if(verbose) print("diagram.station")
  yrs <- as.numeric(rownames(table(year(x))))
  d <- dim(x)
  #print(yrs)
  ny <- length(yrs)
  j <- 1:ny
  Z <- matrix(rep(NA,ny*365),ny,365)
  col <- rgb(j/ny,abs(sin(pi*j/ny)),(1-j/ny),0.2)
  class(x) <- "zoo"
  
  if ( (attr(x,'unit') == "deg C") | (attr(x,'unit') == "degree Celsius") )
    unit <- expression(degree*C) else
      unit <- attr(x,'unit')
  main <- NULL
  eval(parse(text=paste("main <- expression(paste('Seasonal evaution: ',",
                        attr(x,'variable'),"))")))
  if (plot) {
    if (new) dev.new()
    par(bty="n")
  }
  z <- coredata(x)
  if (plot) {
    if (is.T(x)) ylab <- expression(T*(degree*C)) else 
      ylab <- paste(varid(x),' (',esd::unit(x),')',sep='')
    plot(c(0,365),1.25*range(z,na.rm=TRUE),
         type="n",xlab="",
         main=main,
         sub=attr(x,'location'),ylab=ylab)
    grid()
  }
  for (i in 1:ny) {
    y <- window(x,start=as.Date(paste(yrs[i],'-01-01',sep='')),
                end=as.Date(paste(yrs[i],'-12-31',sep='')))
    t <- julian(index(y)) - julian(as.Date(paste(yrs[i],'-01-01',sep='')))
    i1 <- is.element(0:364,t)
    i2 <- is.element(t,0:364)
    Z[i,i1] <- y[i2]
    if (plot) {
      if (is.null(d)) points(t,coredata(y),lwd=2,col=col[i],pch=19,cex=0.5) else
        points(rep(t,d[2]),coredata(y),lwd=2,col=col[i],pch=19,cex=0.5)
    }
  }
  if (!is.null(it)) {
    y <- window(x,start=as.Date(paste(it,'-01-01',sep='')),
                end=as.Date(paste(it,'-12-31',sep='')))
    t <- julian(index(y)) - julian(as.Date(paste(it,'-01-01',sep='')))
  }
  if (plot) {
    if (is.null(d)) points(t,coredata(y),col="black",cex=0.7) else
      points(rep(t,d[2]),coredata(y),col="black",cex=0.7)
    par0 <- par()
    par(new=TRUE,fig=c(0.70,0.85,0.70,0.85),mar=c(0,3,0,0),
        cex.axis=0.7,yaxt="s",xaxt="n",las=1)
    colbar <- rbind(1:ny,1:ny)
    image(1:2,yrs,colbar,col=col)
    par(fig=par0$fig,mar=par0$mar,cex.axis=par0$cex.axis,
        yaxt=par0$yaxt,xaxt=par0$xaxt,las=par0$las)
  }
  rownames(Z) <- yrs
  invisible(t(Z))
}
