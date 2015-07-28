## Author     K. Parding
## Last update   24.04.2015
## Tools for analyzing trajectory objects

season.trajectory <- function(x) {
  stopifnot(!missing(x), inherits(x,"trajectory"))
  mn <- month(x)
  mlist <- c(12,1:11)
  slist <- c(1,1,1,2,2,2,3,3,3,4,4,4)
  sn <- sapply(mn,function(x) slist[mlist==x])
  invisible(sn)
}

param.trajectory <- function(x,param=NULL,FUN='mean') {
  stopifnot(!missing(x), inherits(x,"trajectory"))
  if (is.null(param)) {
    param <- colnames(x)[(!(colnames(x) %in%
          c('lon','lat','start','end','n')))]
    if (length(param)>1) param <- param[1]
    if (length(param)==0) param <- 'n'
  } else {
    if (!(param %in% colnames(x))) {
      print(paste('Input error! Unknown param',param))
      param <- 'n'
    }
  }
  if (!is.null(FUN) & !is.na(FUN)) {
    if (FUN=='first') FUN <- function(x) x[1] else
    if (FUN!='first' & param=='lon') FUN <- fnlon(FUN)
  }
  y <- x[,colnames(x)==param]
  if (sum(colnames(x)==param)>1) {
    y <- apply(x[,colnames(x)==param],1,FUN)
  }
  invisible(y)
}
  
sort.trajectory <- function(x) {
  stopifnot(!missing(x), inherits(x,"trajectory"))
  if (any('sorted' %in% attr(x,'aspect'))) {
    invisible(x)
  } else {
    start <- x[,colnames(x)=='start']
    date <- strptime(start,format="%Y%m%d%H")
    while (any(duplicated(date))) {
      date[duplicated(date)] <- date[duplicated(date)]+60
    }
    y <- x[order(date),]
    y <- attrcp(x,y)
    attr(y,'aspect') <- c('sorted',attr(x,'aspect'))
    class(y) <- class(x)
    invisible(y)
  }
}

anomaly.trajectory <- function(x,type='first',param=c('lon','lat'),
                               verbose=FALSE) {
  stopifnot(!missing(x), inherits(x,"trajectory"))

  if(is.null(param)) {
    param <- unique(colnames(x)[!(colnames(x) %in% c('start','end','n'))])
  } else if (any(!param %in% colnames(x))) {
    print(paste('Warning! Input error: param',
     paste(param[!param %in% colnames(x)],collapse=" "),"missing"))
  }
  if(verbose) print(param)

  y <- x
  y <- attrcp(x,y)  
  if(any('anomaly' %in% attr(x,'aspect')) &
     any(param %in% attr(x,'mean'))) x <- anomaly2trajectory(y)
  m <- vector(mode="list",length=length(param))
  names(m) <- param

  for (p in param) {
    i <- which(colnames(x)==p)
    xi <- x[,i]
    if (p=='lon') xi <- t(apply(x[,i],1,lontrack))
    if (type=='first') {
      m[param==p] <- list(xi[1])
      p.anomaly <- apply(xi,1,function(x) x-x[1])
    } else if (type=='mean') {
      m[param==p] <- list(mean(x[,i]))
      p.anomaly <- apply(x,1,function(x) xi-m[param==p][[1]])
    }
    if(verbose) print(p)
    if(verbose) print(dim(x[,i]))
    if(verbose) print(dim(p.anomaly))
    if(verbose) print(x[1,i])
    if(verbose) print(t(p.anomaly)[1,])
    y[,i] <- t(p.anomaly)
  }
  
  if(any('anomaly' %in% attr(y,'aspect')) &
  any(!names(attr(y,'mean')) %in% param)) {
    p <- names(attr(y,'mean'))
    p <- p[!p %in% c(param)]
    m0 <- vector(mode="list",length=length(p))
    names(m0) <- p
    for (q in p) {
      eval(parse(text=paste("m0[param=",q,"] <- attr(y,'mean')$",q,sep="")))
    }
    m <- c(m,m0)
  } 
  attr(y,'mean') <- m
  attr(y,'aspect') <- unique(c('anomaly',attr(x,'aspect')))
  attr(y,'history') <- history.stamp(x)
  invisible(y)
}

anomaly2trajectory <- function(x,verbose=FALSE) {
  stopifnot(!missing(x), inherits(x,"trajectory"))
  if (any("anomaly" %in% aspect(x))) {
    if(verbose) print(names(attr(x,'mean')))
    for (i in 1:length(attr(x,'mean'))) {
      param.i <- names(attr(x,'mean'))[i]
      mean.i <- unlist(attr(x,'mean')[i])
      if(verbose) print(paste(i,param.i,'length:',
                     paste(length(mean.i),collapse=" ")))
      if (length(mean.i)==1) {
        x[,colnames(x)==param.i] <- x[,colnames(x)==param.i] + mean.i
      } else if (length(mean.i)==sum(colnames(x)==param.i)) {
        x[,colnames(x)==param.i] <- x[,colnames(x)==param.i] +
         t(matrix( rep(array(mean.i),dim(x)[1]),
                  sum(colnames(x)==param.i), dim(x)[1]))
      } else if (length(mean.i)==dim(x)[1]) {
        x[,colnames(x)==param.i] <- x[,colnames(x)==param.i] +
         matrix( rep(array(mean.i),sum(colnames(x)==param.i)),
                dim(x)[1], sum(colnames(x)==param.i)) 
      }
    }
  }
  if (any('lon' %in% names(attr(x,'mean')))) {
    x[,colnames(x)=='lon'] <- t(apply(x[,colnames(x)=='lon'],1,lon2dateline))
  }
  attr(x,'aspect') <- attr(x,'aspect')[attr(x,'aspect')!='anomaly']
  invisible(x)
}

polyfit.trajectory <- function(x) {
  stopifnot(is.trajectory(x))
  if(!('anomaly' %in% attr(x,'aspect'))) x <- anomaly.trajectory(x)
  ilon <- which(colnames(x)=='lon')
  ilat <- which(colnames(x)=='lat')
  OK <- apply(x[,ilon],1,function(y) !any(y>0) | !any(y<0) |
              (mean(y[y>0])-mean(y[y<0]))<120)
  if(any(!OK)) {lon <- x[!OK,ilon]
             lon[lon<0] <- lon[lon<0]+360
             x[!OK,ilon] <- lon}
  latfit <- apply(x,1,function(y) polyfit(y[ilon],y[ilat]))
  coefs <- apply(x,1,function(y) attr(polyfit(y[ilon],y[ilat]),'model'))
  z <- x
  z[,ilat] <- t(latfit)
  attr(z,'aspect') <- unique(c("pfit",attr(x,'aspect')))
  attr(z,'coefs') <- coefs
  return(z)
}

polyfit <- function(x,y=NULL) {
  if(all(is.null(y))) y <- x; x <- seq(1,length(x))
  pfit <- lm(y ~ I(x) + I(x^2) + I(x^3) + I(x^4) + I(x^5))
  z <- predict(pfit)
  attr(z,'model') <- pfit$coef
  return(z)
}

count.trajectory <- function(x,it=NULL,is=NULL,by='year') {
  y <- subset(x,it=it,is=is)
  t <- as.Date(strptime(y[,colnames(y)=="start"],format="%Y%m%d%H"))
  if (by=='year') {
    fmt <- "%Y"
    if (inherits(y,'season')) {
      cls <- 'season'
      unit <- 'events/season'
    } else {
      cls <- 'annual'
      unit <- 'events/year'
    }
  } else if (by %in% c('month','4seasons')) {
    fmt <- "%Y%m%d"
    t <- as.yearmon(t)
    cls <- 'month'
    unit <- 'events/month'
  } else if (by=='day') {
    fmt <- "%Y%m%d"
    cls <- 'day'
    unit <- 'events/day'
  }
  d <- strftime(t,format=fmt)
  n <- table(d)
  if (by=='year') {dn <- dimnames(n)$d
  } else dn <- as.Date(strptime(dimnames(n)$d,format=fmt))
  nz <- zoo(n,order.by=dn)
  class(nz) <- c(cls,"zoo")
  if (by=='4seasons') nz <- as.4seasons(nz,FUN=sum)
  attrcp(y,nz)
  attr(nz,'longname') <- paste(attr(x,'longname'),'event count',sep=', ')
  attr(nz,'unit') <- unit
  invisible(nz)
}

approxlon <- function(lon,n=10) {
  x <- approx(lontrack(lon),n=n)
  x$y <- lon2dateline(x$y)
  return(x)
}

lontrack <- function(lon) {
  lon0 <- lon2dateline(lon)
  n <- length(lon)
  if (any(abs(diff(lon))>200)) {
    signs <- sign(lon0); signs[signs==0] <- 1
    i.change <- which(abs(diff(signs))>0 | abs(diff(lon))>200)
    lon <- lon2greenwich(lon)
    for (i in seq(1,length(i.change))) {
      i1 <- i.change[i]
      if (any(abs(diff(lon)[(i1-1):min(i1+1,n-1)])>200)) {
      	if (abs(diff(lon)[i1])<200) {
      	  if (i1>1 & abs(diff(lon)[max(1,i1-1)])>200) i1<i1-1 else i1<-i1+1
      	}
      	if (i<length(i.change)) i2<-i.change[i+1] else i2<-n
        add <- round(-diff(lon)[i1]/360)*360
        j <- sign(lon0)!=sign(lon0[i1]) & seq(1,n) %in% seq(i1,i2)
        lon[j] <- lon[j] + add
      }
    }
    if (mean(lon)>360) { lon <- lon-360
    } else if (mean(lon) < -180) lon <- lon+360
  }
  invisible(lon)	
}

lon2dateline <- function(lon) {
  while (any(lon > 180) | any(lon < -180)) {
    lon[lon > 180] = lon[lon > 180] - 360
    lon[lon < -180] = lon[lon < -180] + 360
  }
  invisible(lon)
}

lon2greenwich <- function(lon) {
  while (any(lon > 360) | any(lon < 0)) {
    lon[lon < 0] = lon[lon < 0] + 360
    lon[lon > 360] = lon[lon > 360] - 360
  }
  invisible(lon)
}

fnlon <- function(FUN=mean) {
  fn <- function(lon) {
    lon <- lon2dateline(lon)
    if (any(abs(diff(lon))>200)) {
      lon <- lon2greenwich(lon)
      x <- FUN(lon)
      return(lon2dateline(x))
    } else {
      return(FUN(lon))
    }
  }
  return(fn)
}
