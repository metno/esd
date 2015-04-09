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

length.storm <- function(x) {
  n <- x[,colnames(x)=='n']
  invisible(n)
}

lat.storm <- function(x,FUN=mean) {
  lat <- apply(x[,colnames(x)=='lat'],1,FUN)
  invisible(lat)
}

lon.storm <- function(x,FUN=mean) {
  lon <- apply(x[,colnames(x)=='lon'],1,FUN)
  invisible(lon)
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
    y <- x[order(date),]
    y <- attrcp(x,y)
    attr(y,'aspect') <- c('sorted',attr(x,'aspect'))
    class(y) <- class(x)
    invisible(y)
  }
}

anomaly.storm <- function(x,param=c('lon','lat','slp')) {
  if (any('lon' %in% param)) {
    i <- which(colnames(x)=='lon')
    dateline <- apply(x,1,function(x) (max(x[i])-min(x[i]))>180 )
    lon <- x[dateline,i]
    lon[lon<0] <- lon[lon<0]+360
    x[dateline,i] <- lon
  }
  m <- vector(mode="list", length=length(param))
  names(m) <- param
  for (p in param) {
    i <- which(colnames(x)==p)
    if (p=='slp') {
      p.anomaly <- x[,i]-mean(x[,i])
      m[param==p] <- list(mean(x[,i]))
    } else {
      p.anomaly <- apply(x,1,function(x) x[i]-x[i[1]])
      m[param==p] <- list(x[,i[1]])
    }
    x[,i] <- t(p.anomaly)
  }
  y <- x
  attr(y,'mean') <- m
  y <- attrcp(x,y)
  attr(y,'aspect') <- c('anomaly',attr(x,'aspect'))
  attr(y,'history') <- history.stamp(x)
  invisible(y)
}

count.storm <- function(x,it=NULL,is=NULL,by='year') {
  y <- subset.storm(x,it=it,is=is)
  if (by=='year') {
    i <- seq(min(which(min(month.storm(y))==month.storm(y))),
         max(which(max(month.storm(y))==month.storm(y))))
    y <- subset.storm(y,it=i)
  }
  t <- strptime(y[,colnames(y)=="start"],format="%Y%m%d%H")
  if (by=='year') {
    fmt <- "%Y"
    cls <- 'annual'
  } else if (by=='month') {
    fmt <- "%Y%m%d"
    t <- as.yearmon(t)
    cls <- 'month'
  } else if (by=='day') {
    fmt <- "%Y%m%d"
    cls <- 'day'
  }
  d <- strftime(t,format=fmt)
  n <- table(d)
  dn <- as.Date(strptime(dimnames(n)$d,format=fmt))
  nz <- zoo(n,order.by=dn)
  class(nz) <- c(class(nz),cls)
  attrcp(y,nz)
  attr(nz,'longname') <- 'storm count'
  attr(nz,'unit') <- 'storms/year'
  invisible(nz)
}

approx.lon <- function(lon,n=10) {
  if (!any(lon>0)|!any(lon<0)|(mean(lon[lon>0])-mean(lon[lon<0]))<120){
    x <- approx(lon,n=n)
  } else {
    lon[lon<0] <- lon[lon<0]+360
    x <- approx(lon,n=n)
    x$y[x$y>180] <- x$y[x$y>180]-360
  }
  return(x)
}

mean.lon <- function(lon) {
  if (!any(lon>0)|!any(lon<0)|(mean(lon[lon>0])-mean(lon[lon<0]))<120){
    x <- mean(lon)
  } else {
    lon[lon<0] <- lon[lon<0]+360
    x <- mean(lon)
    if (x>180) x <- x-360
  }
  return(x)
}

# Takes too long! 
polyfit.storm <- function(X) {
  ilon <- colnames(X,'lon')
  ilat <- colnames(X,'lat')
  Z <- apply(X,1,function(x) pfit(x[ilon],x[ilat]))
  return(Z)
}

pfit <- function(lon,lat) {
  OK <- !any(lon>0) | !any(lon<0) | (mean(lon[lon>0])-mean(lon[lon<0]))<120
  if (!OK) lon[lon<0] <- lon[lon<0]+360
  pfit <- lm(lat ~ I(lon) + I(lon^2) + I(lon^3) + I(lon^4) + I(lon^5))
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
