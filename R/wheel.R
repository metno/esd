#' InfoGraphics
#'
#' Various functions for visual display of data and statistics
#'
#' \code{wheel} shows the seasonal cycle with different colors for different years
#'
#' @aliases wheel.station wheel.spell
#' @seealso graph visprob conf vis diagram cumugram scatter plot map
#' 
#' @param x an input object of class 'station' or 'spell'
#' @param y an input object of class 'station' or 'spell'
#' @param new if new create new graphic device
#' @param lwd relative line width
#' @param col color of line
#' @param type 'spiky' or 'flowy'
#' @param bg background color
#' @param verbose a boolean; if TRUE print information about progress
#' @param \dots additional arguments
#'
#' @examples
#' data(bjornholt)
#' wheel(bjornholt, new=FALSE)
#'
#' @export
wheel <- function(x,...) UseMethod("wheel")

# S3 method for 'station' object
#' @exportS3Method
#' @export wheel.station
wheel.station <- function(x,y=NULL,new=TRUE,lwd=2,col=NULL,type=NULL,
                          bg="grey90",verbose=FALSE,...) {

  if (verbose) print('wheel.station')
  ## Copied from visprob.station.precip:
  ## If y is provided, synchronise the two time series:
  if (!is.null(y)) {
    y <- subset(y,it=c(start(x),end(x)))
    x <- subset(x,it=c(start(y),end(y)))
    years <- annual(y)
  } else years <- year(annual(x))
  mx <- max(abs(coredata(x)),na.rm=TRUE)
#  r <- mean(coredata(x),na.rm=TRUE)
  r <- mx
  #jday <- julian(zoo(x))
  MD <- month(x)*100 + day(x)
  #REB2015-08-20 years <- as.integer(rownames(table(year(x))))
  #print(years)
  ny <- length(years)
  md <- as.integer(rownames(table(MD)))
  #print(md)
  m <- length(md)
  if (new) dev.new()
  par(bty="n",xaxt="n",yaxt="n") -> par0
  plot(c(-1.2,1.2)*(r+mx),c(-1.2,1.2)*(r+mx),type="n",
       xlab="",ylab="",
       main=paste("Seasonal 'wheel' for",attr(x,'location')[1]),
       sub=paste(attr(x,'variable'),collapse=', '))
  rect(-(r+mx),(r+mx),0,0,col=bg,border=NA)
  rect(0,0,(r+mx),-(r+mx),col=bg,border=NA)
  lines(c(-1,1)*(r+mx),c(0,0))
  lines(c(0,0),c(-1,1)*(r+mx))
  text(0,1.2*(r+mx),"January")
  text(0,-1.2*(r+mx),"July")
  text(1.2*(r+mx),0,"April",srt=-90)
  text(-1.2*(r+mx),0,"September",srt=90)
  w <- seq(0,2*pi,length=m)
  j <- 1:ny
  col <- rgb(j/ny,abs(sin(pi*j/ny)),(1-j/ny),0.2)
  ## REB 2015-08-20: if y is given, use it for setting the colours.
  ## e.g. according to an index like NINO3.4
  srtc <- order(years)
  col <- col[srtc]
  years <- year(annual(x))
  if (verbose) {print(srtc); print(col)}
  if (is.null(type) & !is.T(x)) type <- 'spiky' else type <- 'flowy'
  if(verbose) print(type)
  
  if (type=='spiky') {
    for (i in 1:m) {
     wj <- -w[i]
     s <- sin(wj)
     c <- cos(wj)
     ii <- is.element(MD,md[i])
    #print(sum(ii))
     yr <- year(x)[ii]
     y <- coredata(x)[ii]
     yn <- min(abs(y),na.rm=TRUE)
     srt <- order(abs(y),decreasing=TRUE)
    #print(srt); print(yr)
      yr <- yr[srt]; y <- y[srt]
      for (j in 1:sum(is.finite(y))) {
       jj <- (1:length(years))[is.element(years,yr[j])]
        lines(-c((r+yn)*s,(r+y[j])*s),
               c((r+yn)*c,(r+y[j])*c),
              lwd=lwd,col=col[jj])
     }  
   }
  } else {
    yr <- year(x); mo <- month(x); dy <- day(x)
    for (i in 2:length(index(x))) {
      jday1 <- as.numeric(as.Date(paste(yr[i-1],mo[i-1],dy[i-1],sep='-'))-as.Date(paste(yr[i-1],'01-01',sep='-')))
      jday2 <- as.numeric(as.Date(paste(yr[i],mo[i],dy[i],sep='-'))-as.Date(paste(yr[i],'01-01',sep='-')))
      theta1 <- 0.5*pi - 2*pi*jday1/366
      theta2 <- 0.5*pi - 2*pi*jday2/366
      jj <- (1:length(years))[is.element(years,yr[i])]
      lines(c((r+coredata(x)[i-1])*cos(theta1),(r+coredata(x)[i])*cos(theta2)),
            c((r+coredata(x)[i-1])*sin(theta1),(r+coredata(x)[i])*sin(theta2)),lwd=lwd,col=col[jj])
      if (verbose & (i %% 100 ==0)) print(c(theta1,theta2,
                           coredata(x)[i-1]*cos(theta1),coredata(x)[i]*cos(theta2),
                           coredata(x)[i-1]*sin(theta1),coredata(x)[i]*sin(theta2)))
    }
  }
  #clim <- coredata(climatology(x))
  #r <- approx(1:length(clim),clim,xout = 1:360)$y
  #lines(r[1:360]*cos(pi*(1:360)/180),r[1:360]*sin(pi*(1:360)/180),lwd=2,col=rgb(0.5,0.5,0.5,0.2))
  
  par(new=TRUE,fig=c(0.05,0.15,0.05,0.2),mar=c(0,3,0,0),
      cex.axis=0.7,yaxt="s",xaxt="n",las=1)
  colbar <- rbind(1:ny,1:ny)
  #print(years)
  image(1:2,years,colbar,col=col)
}

# S3 method for 'spell' object

#' @exportS3Method
#' @export wheel.spell
wheel.spell <- function(x,y=NULL,new=TRUE,lwd=2,col=NULL,verbose=FALSE,...) {

    if (verbose) print('wheel.spell')
    ## Copied from visprob.station.precip:
  ## If y is provided, synchronise the two time series:
  if (!is.null(y)) {
    y <- subset(y,it=c(start(x),end(x)))
    x <- subset(x,it=c(start(y),end(y)))
    years <- annual(y)
  } else years <- year(annual(x))
  mx <- max(abs(coredata(x)),na.rm=TRUE)
#  r <- mean(coredata(x),na.rm=TRUE)
  r <- mx
  #jday <- julian(zoo(x))
  MD <- month(x)*100 + day(x)
  #REB2015-08-20 years years <- as.integer(rownames(table(year(x))))
  #print(years)
  ny <- length(years)
  md <- as.integer(rownames(table(MD)))
  #print(md)
  m <- length(md)
  dev.new()
  par(bty="n",xaxt="n",yaxt="n") -> par0
  plot(c(-1.2,1.2)*(r+mx),c(-1.2,1.2)*(r+mx),type="n",
       xlab="",ylab="",
       main=paste("Seasonal 'wheel' for",attr(x,'location')[1]),
       sub=paste(attr(x,'variable'),collapse=', '))
  rect(-(r+mx),(r+mx),0,0,col="grey90",border=NA)
  rect(0,0,(r+mx),-(r+mx),col="grey90",border=NA)
  lines(c(-1,1)*(r+mx),c(0,0))
  lines(c(0,0),c(-1,1)*(r+mx))
  text(0,1.2*(r+mx),"January")
  text(0,-1.2*(r+mx),"July")
  text(1.2*(r+mx),0,"April",srt=-90)
  text(-1.2*(r+mx),0,"September",srt=90)
  w <- seq(0,2*pi,length=m)
  j <- 1:ny
  col <- rgb(j/ny,abs(sin(pi*j/ny)),(1-j/ny))
  ## REB 2015-08-20: if y is given, use it for setting the colours.
  ## e.g. according to an index like NINO3.4
  srtc <- order(years)
  col <- col[srtc]
  years <- year(annual(x))
  if (verbose) {print(srtc); print(col)}
  
  for (i in 1:m) {
    wj <- -w[i]
    s <- sin(wj)
    c <- cos(wj)
    ii <- is.element(MD,md[i])
    #print(sum(ii))
    yr <- year(x)[ii]
    y <- coredata(x)[ii,1]
    z <- coredata(x)[ii,2]
    yn <- min(abs(y),na.rm=TRUE)
    srt <- order(abs(y),decreasing=TRUE)
    #print(srt); print(yr)
    yr <- yr[srt]; y <- y[srt]
    for (j in 1:sum(is.finite(y))) {
      jj <- (1:length(years))[is.element(years,yr[j])]
      lines(-c((r+yn)*s,(r+y[j])*s),
             c((r+yn)*c,(r+y[j])*c),
             lwd=lwd,col=col[jj])
      lines(-c((r+yn)*s,(r-z[j])*s),
             c((r+yn)*c,(r-z[j])*c),
             lwd=lwd,col=col[jj])
    }  
  }
  lines(r*cos(pi*(1:360)/180),r*sin(pi*(1:360)/180),lwd=2)
  
  par(new=TRUE,fig=c(0.05,0.15,0.05,0.2),mar=c(0,3,0,0),
      cex.axis=0.7,yaxt="s",las=1)
  colbar <- rbind(1:ny,1:ny)
  #print(years)
  image(1:2,years,colbar,col=col)
}

