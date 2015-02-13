## Author     K. Parding
## Last update   12.02.2015
## Tools for analyzing IMILAST stormtrack files

year.storm <- function(x) {
  start <- strptime(x[,colnames(x)=='start'],format="%Y%m%d%H")
  yr <- as.numeric(strftime(start,format="%Y"))
  invisible(yr)
}

month.storm <- function(x) {
  start <- strptime(x[,colnames(x)=='start'],format="%Y%m%d%H")
  mn <- as.numeric(strftime(start,format="%m"))
  invisible(mn)
}

season.storm <- function(x) {
  mn <- month.storm(x)
  mlist <- c(12,1:11)
  slist <- c(1,1,1,2,2,2,3,3,3,4,4,4)
  sn <- sapply(mn,function(x) slist[mlist==x])
  invisible(sn)
}
  
slp.storm <- function(x,FUN=min) {
  slp <- apply(x[,colnames(x)=='slp'],1,FUN)
  invisible(slp)
}

sort.storm <- function(x) {
  if (any('sorted' %in% attr(x,'aspect'))) {
    invisible(x)
  } else {
    start <- x[,colnames(x)=='start']
    date <- strptime(start,format="%Y%m%d%H")
    while (sum(duplicated(date))>0) {
      date[duplicated(date)] <- date[duplicated(date)]+60
    }
    y <- zoo(x,order.by=date)
    y <- attrcp(x,y)
    attr(y,'aspect') <- c('sorted',attr(x,'aspect'))
    class(y) <- c('zoo',class(x))
    invisible(y)
  }
}

anomaly.storm <- function(x) {
  if (any('anomaly' %in% attr(x,'aspect'))) {
    invisible(x)
  } else {
    ilat <- which(colnames(x)=='lat')
    ilon <- which(colnames(x)=='lon')
    lat.anomaly <- apply(x,1,function(x) x[ilat]-x[ilat[1]])
    lon.anomaly <- apply(x,1,function(x) x[ilon]-x[ilon[1]])
    x[,ilat] <- t(lat.anomaly)
    x[,ilon] <- t(lon.anomaly)
    y <- x
    y <- attrcp(x,y)
    attr(y,'aspect') <- c('anomaly',attr(x,'aspect'))
    attr(y,'history') <- history.stamp(x)
    invisible(y)
  }
}

count.storm <- function(x,by='year') {
  t <- strptime(x[,colnames(x)=="start"],format="%Y%m%d%H")
  if (by=='year') {
    fmt <- "%Y"
  } else if (by=='month') {
    fmt <- "%Y%m%d"
    t <- as.yearmon(t)
  }
  d <- strftime(t,format=fmt)
  n <- table(d)
  dn <- as.Date(strptime(dimnames(n)$d,format=fmt))
  nz <- zoo(n,order.by=dn)
  invisible(nz)
}

  
# Takes too long! 
polyfit.storm <- function(X) {
  Z <- apply(X,1,function(x) pfit(x[1:10],x[11:20])) 
  return(Z)
}

pfit <- function(x,y) {
  pfit <- lm(y ~ I(x) + I(x^2) + I(x^3) + I(x^4) + I(x^5))
  z <- predict(pfit)
  attr(z,'model') <- pfit$coef
  return(z)
}


# Should count storms in equal area bins (hexagonal?), not lat-lon grid
# Use binscatter.bin?
## bin.storm <- function(x,dlon=5,dlat=2,digits=6) {
##   lats <- x[,colnames(x)=='lat']
##   lons <- x[,colnames(x)=='lon']
##   D <- dim(lons)
##   if (length(D)==2) {
##     lons <- matrix(lons,1,D[1]*D[2])
##     lats <- matrix(lats,1,D[1]*D[2])
##   }
##   lons <- round(lons/dlon)*dlon
##   lats <- round(lats/dlat)*dlat
##   xy <- xy.coords(lons,lats,'lon','lat','')
##   tt <- xyTable(xy, digits=digits)
##   invisible(tt)
## }
