# Rasmus Benestad
# Estimate the naomaly and climatology
# Store the monthly climatology as an attribute (12,np)

anomaly <-function(x,...) UseMethod("anomaly")

anomaly.default <- function(x,...) {
  if (inherits(x,'annual')) y <- anomaly.annual(x,...) else
  if (inherits(x,'month')) y <- anomaly.month(x,...) else
  if (inherits(x,'day')) y <- anomaly.day(x,...) else
  if (inherits(x,'season')) y <- anomaly.season(x,...) else
  y <- as.annual(x,...)
  return(y)
}

anomaly.field <- function(x,...) {
  stopifnot(inherits(x,"field"))
  x <- as.anomaly(x)
  return(x)
}

anomaly.comb <- function(x,...) {
  stopifnot(inherits(x,"field"),inherits(x,"comb"))
  y <- anomaly(x)
  n.apps <- attr(x,'n.apps')
  for (i in 1:n.apps) {
    eval(parse(text=paste("z <- attr(x,'appendix.",i,"')",sep="")))
    Z <- anomaly(z)
    eval(parse(text=paste("Z -> attr(x,'appendix.",i,"')",sep="")))
  }
  y <- attrcp(x,y)
  n.apps -> attr(x,'n.apps')
  attr(x,'history') <- history.stamp(x)
  return(y)
}


anomaly.station <- function(x,...) {
#  t <- index(X)[1:2]
#  datetype <- class(t)
#  if (!is.null(attr(x,'anomaly'))) {
#    orig <- coredata(X)
#    x <- zoo(attr(X,'anomaly'),order.by=index(X))
#    nattr <- softattr(X)
#    for (i in 1:length(nattr))
#      attr(x,nattr[i]) <- attr(X,nattr[i])
#    eval(parse(text=paste("attr(x,'",attr(X,'aspect'),"') <- orig")))
#    attr(x,'aspect') <- 'anomaly'
#    return(x)
#  }
#
#  if (datetype=="Date") {
#    dy <- diff(as.numeric(format(t,'%Y')))
#    dm <- diff(as.numeric(format(t,'%m')))
#    dd <- diff(as.numeric(format(t,'%d')))
#  } else if (datetype=="numeric") {
#    dy <- 1; dm <- 0; dd <- 0
#  }
#  if ((dy==1) & (dy==0) & (dd==0))
#    x <- anomaly.yearly(X) else
#  if ((dy==0) & (dm==1) & (dd==0))
#    x <- anomaly.monthly(X) else 
#  if ((dy==0) & (dm==0) & (dd==1))
#    x <- anomaly.daily(X)
#  attr(x,'history') <- history.stamp(X)
  x <- anomaly.default(x,...)
  return(x)
}

anomaly.annual <- function(x,ref=1961:1990,verbose=FALSE) {
  if (verbose) print('anomaly.annual')
  X <- x;  x <- coredata(X)
  t <- index(X)
  datetype <- class(t)
  if (verbose) print(datetype)
  if (datetype=="Date") years <- year(X) else
  if (datetype=="numeric") years <- t
  if (is.null(dim(x))) 
    clim <- mean(x[is.element(years,ref)],na.rm=TRUE) else
    clim <- apply(x,2,FUN=mean,na.rm=TRUE)
  if (verbose) print(clim)
  x <- t(t(x) - clim)
  x <- zoo(x,order.by=index(X))
  x <- attrcp(X,x)
  #nattr <- softattr(X)
  #for (i in 1:length(nattr))
  #  attr(x,nattr[i]) <- attr(X,nattr[i])
  attr(x,'climatology') <- clim
  attr(x,'aspect') <- 'anomaly'
  attr(x,'history') <- history.stamp(X)
  class(x) <- class(X)
  return(x)
}

anomaly.month <- function(x,ref=NULL,verbose=FALSE) {

  anomaly.month1 <- function(x,yr=NULL,ref=NULL) {
# This function computes the anomalies by removing the 12-month seasonal cycle
    l <- length(x); n <- ceiling(l/12)
    ## base-line period
    if (!is.null(yr) & !is.null(ref)) iref <- is.element(yr,ref) else
                                      iref <- rep(TRUE,n)
  # If the record is not full years, pad the extra months of the last year
    pad <- l %% 12
    if (pad>0) x <- c(rep(NA,pad),x)
  #Fast way to compute the climatology: clim
    dim(x) <- c(12,n)
    clim <- rowMeans(x[,iref],na.rm=TRUE)
    x <- c(x - clim)
    if (pad>0) x <- x[-(1:pad)]
    return(x)
  }
  
  X <- x
  if (verbose) print('anomaly.month')
  t <- index(x); yr <- year(x)
  
  if (is.null(dim(x))) Y <- anomaly.month1(coredata(x),yr,ref=ref) else
                       Y <- apply(coredata(x),2,FUN='anomaly.month1',yr=yr,ref=ref)
  y <- Y
  x <- zoo(y,order.by=t)
  x <- attrcp(X,x)
  #nattr <- softattr(X)
  #for (i in 1:length(nattr))
  #  attr(x,nattr[i]) <- attr(X,nattr[i])
  attr(x,'climatology') <- (X - Y)[1:12,]
  attr(x,'aspect') <- 'anomaly'
  class(x) <- class(X)
  return(x)
}


anomaly.season <- function(x,ref=NULL,verbose=FALSE) {

  anomaly.season1 <- function(x,yr=NULL,ref=NULL) {
# This function computes the anomalies by removing the 12-month seasonal cycle
    l <- length(x); n <- ceiling(l/4)
    ## base-line period
    if (!is.null(yr) & !is.null(ref)) iref <- is.element(yr,ref) else
                                    iref <- rep(TRUE,n)
    ## If the record is not full years, pad the extra months of the last year
    pad <- l %% 4
    if (pad>0) x <- c(rep(NA,pad),x)
    ##Fast way to compute the climatology: clim
    dim(x) <- c(4,n)
    clim <- rowMeans(x,na.rm=TRUE)
    x <- c(x - clim)
    if (pad>0) x <- x[-(1:pad)]
    return(x)
}
  X <- x
  if (verbose) print('anomaly.season')
  t <- index(x); yr <- year(x)
  if (is.null(dim(x))) y <- anomaly.season1(coredata(x),yr,ref=ref) else
                       y <- apply(coredata(x),2,FUN='anomaly.season1',yr=yr,ref=ref)
  x <- zoo(y,order.by=t)
  x <- attrcp(X,x)
  #nattr <- softattr(X)
  #for (i in 1:length(nattr))
  #  attr(x,nattr[i]) <- attr(X,nattr[i])
  attr(x,'climatology') <- clim
  attr(x,'aspect') <- 'anomaly'
  class(x) <- class(X)
  return(x)
}


anomaly.day <- function(x,ref=NULL,verbose=FALSE) {

  anomaly.day.1 <- function(x,t0,t,ref=NULL) {
    ## One station 
    c1 <- cos(pi*t0/365.25); s1 <- sin(pi*t0/365.25)
    c2 <- cos(2*pi*t0/365.25); s2 <- sin(2*pi*t0/365.25)
    c3 <- cos(3*pi*t0/365.25); s3 <- sin(3*pi*t0/365.25)
    c4 <- cos(4*pi*t0/365.25); s4 <- sin(4*pi*t0/365.25)
    C1 <- cos(pi*t/365.25);   S1 <- sin(pi*t/365.25)
    C2 <- cos(2*pi*t/365.25); S2 <- sin(2*pi*t/365.25)
    C3 <- cos(3*pi*t/365.25); S3 <- sin(3*pi*t/365.25)
    C4 <- cos(4*pi*t/365.25); S4 <- sin(4*pi*t/365.25)
    cal <- data.frame(y=coredata(x),c1=c1,c2=c2,c3=c3,c4=c4,
                      s1=s1,s2=s2,s3=s3,s4=s4)
    pre <- data.frame(c1=C1,c2=C2,c3=C3,c4=C4,
                      s1=S1,s2=S2,s3=S3,s4=S4)
    i1 <- is.element(year(x),year(x)[1])
    pre1 <- data.frame(c1=C1[i1],c2=C2[i1],c3=C3[i1],c4=C4[i1],
                       s1=S1[i1],s2=S2[i1],s3=S3[i1],s4=S4[i1])
    acfit <- lm(y ~ c1 + s1 + c2 + s2 + c3 + s3 + c4 + s4,data=cal)
    clim <- predict(acfit,newdata=pre)
    y <- zoo(coredata(x) - clim,order.by=index(x))
    return(y)
  }
  
  if (verbose) {print('anomaly.day');print(class(x))}
  yr <- year(x)
  if (is.null(ref)) ref <- seq(min(yr,na.rm=TRUE),max(yr,na.rm=TRUE),by=1)
  t0 <- julian(index(x)[is.element(yr,ref)]) -
    julian(as.Date(paste(yr[is.element(yr,ref)],"-01-01",sep="")))
  t <- julian(index(x)) - julian(as.Date(paste(yr,"-01-01",sep="")))
  if (is.null(dim(x))) 
    y <- anomaly.day.1(x=coredata(x),t0=t0,t=t,ref=ref) else 
    y <- apply(coredata(x),2,FUN='anomaly.day.1',t0=t0,t=t,ref=ref)
  y <- zoo(y,order.by=index(x))
  y <- attrcp(x,y)
  class(y) <- class(x)
  clim <- (x - y)[is.element(yr,yr[1]),] 
  attr(y,'climatology') <- clim
  attr(y,'aspect') <- 'anomaly'
  return(y)
}


climatology <- function(x,...) UseMethod("climatology")

climatology.default <- function(x) {
  x <- as.climatology(x)
  return(x)
}

climatology.field <- function(x) {
  x <- as.climatology(x)
  return(x)
}

climatology.station <- function(x) {
#  x <- X
#  orig <- coredata(X)
#  if (is.null(attr(x,'climatology'))) {
#    t <- index(X)[1:2]
#    dy <- diff(as.numeric(format(t,'%Y')))
#    dm <- diff(as.numeric(format(t,'%m')))
#    dd <- diff(as.numeric(format(t,'%d')))
#    if ((dy==1) & (dy==0) & (dd==0))
#      x <- anomaly.yearly(X) else
#    if ((dy==0) & (dm==1) & (dd==0)) 
#      y <- anomaly.monthly(X) 
#    if ((dy==0) & (dm==0) & (dd==1)) 
#      y <- anomaly.daily(X)
#    clim <- attr(y,'climatology')
#  } else clim <- attr(X,'climatology')
#
#  nc <- length(index(X))%/%length(clim)
#  pad <- length(index(X))%%length(clim)
#  clim <- rep(clim,nc)
#  if (pad>0) clim <- c(clim,clim[1:pad])
#    
#  x <- zoo(clim,order.by=index(X))
#  nattr <- softattr(X)
#  for (i in 1:length(nattr))
#      attr(x,nattr[i]) <- attr(X,nattr[i])
#  eval(parse(text=paste("attr(x,'",attr(X,'aspect'),"') <- orig")))
#  attr(x,'aspect') <- 'climatology'
  x <- as.climatology(x)
  return(x)
}





# Station data can be expressed as PCA where each of the EOFs represent one
# year. The PCs describe the seasonal variations

clim2pca <-function(x,...) UseMethod("clim2pca")

clim2pca.default <- function(x) {
}

clim2pca.month <- function(x) {
  X <- aggregate(x,year)
  ny <- length(x) %/% 12
  nm <- length(x) %% 12
  y <- coredata(x[1:(length(x)-nm)])
  dim(y) <- c(12,ny)
  ok <- is.finite(colMeans(y))
  pca <- svd(y[,ok])
  for (i in 1:12) {
    z <- zoo(pca$v[,i],order.by=index(X))
    if (i == 1) Z <- z else
                Z <- merge(Z,z)
  }
  season <- pca$u
  colnames(season) <- month.abb
  rownames(season) <- paste("pattern",1:12,sep=".")
  attr(Z,'season') <- season
  attr(Z,'d') <- pca$d
  return(Z)
}

clim2pca.day <- function(x) {
}


