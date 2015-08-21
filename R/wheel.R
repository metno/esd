# A function that makes infographics - plots precipitation as a wheel
# with different colours for different years, and different angles for
# different julian days. The largest values are drawn first

wheel <- function(x,...) UseMethod("wheel")

wheel.station <- function(x,y=NULL,new=TRUE,lwd=2,col=NULL,
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
  lines(r*cos(pi*(1:360)/180),r*sin(pi*(1:360)/180),lwd=2)
  
  par(new=TRUE,fig=c(0.05,0.15,0.05,0.2),mar=c(0,3,0,0),
      cex.axis=0.7,yaxt="s",xaxt="n",las=1)
  colbar <- rbind(1:ny,1:ny)
  #print(years)
  image(1:2,years,colbar,col=col)
}

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

