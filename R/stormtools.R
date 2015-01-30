## Author     K. Parding
## Last update   29.01.2015
## Tools for analysing IMILAST stormtrack files

year.stormmatrix <- function(x) {
  start <- strptime(x[,colnames(x)=='start'],format="%Y%m%d%H")
  yr <- as.numeric(strftime(start,format="%Y"))
  invisible(yr)
}

month.stormmatrix <- function(x) {
  start <- strptime(x[,colnames(x)=='start'],format="%Y%m%d%H")
  mn <- as.numeric(strftime(start,format="%m"))
  invisible(mn)
}

season.stormmatrix <- function(x) {
  mn <- month.stormmatrix(x)
  mlist <- c(12,1:11)
  slist <- c(1,1,1,2,2,2,3,3,3,4,4,4)
  sn <- sapply(mn,function(x) slist[mlist==x])
  invisible(sn)
}
  
slp.stormmatrix <- function(x,FUN=min) {
  slp <- apply(x[,colnames(x)=='slp'],1,FUN)
  invisible(slp)
}

stormcount <- function(x,...) UseMethod("stormcount.time")
  
stormcount.time <- function(x,by='year') {
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

stormcount.space <- function(x,dlon=5,dlat=2,digits=6) {
  lats <- x[,colnames(x)=='lat']
  lons <- x[,colnames(x)=='lon']
  D <- dim(lons)
  if (length(D)==2) {
    lons <- matrix(lons,1,D[1]*D[2])
    lats <- matrix(lats,1,D[1]*D[2])
  }
  lons <- round(lons/dlon)*dlon
  lats <- round(lats/dlat)*dlat
  xy <- xy.coords(lons,lats,'lon','lat','')
  tt <- xyTable(xy, digits=digits)
  invisible(tt)
}

  
# Takes forever! What's wrong? I am using apply...
polyfit.stormmatrix <- function(X) {
  Z <- apply(X,1,function(x) polyfit(x[1:10],x[11:20])) 
  return(Z)
}

polyfit <- function(x,y) {
  pfit <- lm(y ~ I(x) + I(x^2) + I(x^3) + I(x^4) + I(x^5))
  z <- predict(pfit)
  attr(z,'model') <- pfit$coef
  return(z)
}

pca.stormmatrix <- function(x) {

}

