## Author A. Mezghani
## Description Contains some rtools ...
## Created 14.11.2014

as.decimal <- function(x=NULL) {
    ## converts from degree min sec format to degrees ...
    ##x is in the form "49°17´38´´"
    if (!is.null(x)) {
        deg <-as.numeric(substr(x,1,2)) 
        min <- as.numeric(substr(x,4,5))
        sec <- as.numeric(substr(x,7,8))     
        x <- deg + min/60 + sec/3600
    }
    return(x)
}


## compute the percentage of missing data in x
missval <- function(x) sum(is.na(coredata(x)))/length(coredata(x))

## compute the quantile 95% of x
q95 <- function(x,na.rm=TRUE) quantile(x,probs=.95,na.rm=na.rm)

## compute the quantile 5% of x
q5 <- function(x,na.rm=TRUE) quantile(x,probs=.05,na.rm=na.rm)

## compute the quantile 5% of x
q995 <- function(x,na.rm=TRUE) quantile(x,probs=.995,na.rm=na.rm)

## compute the quantile 5% of x
q975 <- function(x,na.rm=TRUE) quantile(x,probs=.975,na.rm=na.rm)

## count the number of valid data points
nv <- function(x,...) sum(is.finite(x))

## Compute the coefficient of variation of x
cv <- function(x,na.rm=TRUE) {sd(x,na.rm=na.rm)/mean(x,na.rm=na.rm)}

stand <- function(x) (x - mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE)


## Compute the linear trend
trend.coef <- function(x,...) {
  t <- 1:length(x)
  model <- lm(x ~ t)
  y <- c(model$coefficients[2]*10)
  names(y) <- c("trend.coefficients")
  return(y)
}

## Compute the linear trend
trend.err <- function(x,...) {
  t <- 1:length(x)
  model <- lm(x ~ t)
  y <- c(summary(model)$coefficients[4]*10)
  names(y) <- c("trend.standard.error")
  return(y)
}


## Compute the p-value of the linear trend 
trend.pval <- function(x,...) {
    t <- 1:length(x)
    model <- lm(x ~ t)
    y <- anova(model)$Pr[1]
    names(y) <- c("trend.pvalue")
    return(y)
}

# Wrap-around for lag.zoo to work on station and field objects:
lag.station <- function(x,...) {
  y <- lag(zoo(x),...)
  y <- attrcp(x,y)
  class(y) <- class(x)
  invisible(y)
}

lag.field <- function(x,...) lag.station(x,...)
  
exit <- function() q(save="no")

filt <- function(x,...) UseMethod("filt")

filt.default <- function(x,n,type='ma',lowpass=TRUE) {
  
# A number of different filters using different window
# shapes.
#
# R.E. Benestad, July, 2002, met.no.
#
# ref: Press et al. (1989), Numerical Recipes in Pascal, pp. 466
#library(ts)

# Moving-average (box-car) filter
  ma.filt <- function(x,n) {
    y <- filter(x,rep(1,n)/n)
    y
  }

# Gaussian filter with cut-off at 0.025 and 0.975 of the area.
  gauss.filt <- function(x,n) {
    i <- seq(0,qnorm(0.975),length=n/2)
    win <- dnorm(c(sort(-i),i))
    win <- win/sum(win)
    y <- filter(x,win)
    y
  }

# Binomial filter
  binom.filt <- function(x,n) {
    win <- choose(n-1,0:(n-1))
    win <- win/max(win,na.rm=T)
    win[is.na(win)] <- 1
    win <- win/sum(win,na.rm=T)
    y <- filter(x,win)
    y
  }

# Parzen filter (Press,et al. (1989))
  parzen.filt  <-  function(x,n) {
    j <- 0:(n-1)
    win <- 1 - abs((j - 0.5*(n-1))/(0.5*(n+1)))
    win <- win/sum(win)
    y <- filter(x,win)
    y
  }

# Hanning filter (Press,et al. (1989))
  hanning.filt  <-  function(x,n) {
    j <- 0:(n-1)
    win <- 0.5*(1-cos(2*pi*j/(n-1)))
    win <- win/sum(win)
    y <- filter(x,win)
    y
  }

# Welch filter (Press,et al. (1989))
  welch.filt  <-  function(x,n) {
    j <- 0:(n-1)
    win <- 1 - ((j - 0.5*(n-1))/(0.5*(n+1)))^2
    win <- win/sum(win)
    y <- filter(x,win)
    y
  }

  y <- coredata(x)
  z <- eval(parse(text=paste(type,'filt(y,n)',sep='.')))
  hp <- as.numeric(y - coredata(z))
  if (!is.null(dim(x))) dim(hp) <- dim(x)
  if (lowpass) coredata(x) <- coredata(z) else
               coredata(x) <- hp
  attr(x,'history') <- history.stamp(x)
  return(x)
}
  
figlab <- function(x,xpos=0.001,ypos=0.001) {
  par(new=TRUE,pdx=NA,fig=c(0,1,0,1),xaxt='n',yaxt='n',bty='n',mar=rep(0,4))
  plot(c(0,1),c(0,1),type='n')
  text(xpos,ypos,x,type=2,cex=1.2,pos=4,col='grey30')
}

ensemblemean <- function(x,FUN='rowMeans') {
  if (inherits(x,'pca')) z <- as.station(x) else z <- x
  ## Estimate the ensemble mean
  zz <- unlist(lapply(coredata(z),FUN=FUN))
  zm <- matrix(zz,length(index(z[[1]])),length(z))
  zm <- zoo(zm,order.by=index(z[[1]]))
  zm <- as.station(zm,param=varid(z),unit=unit(z),
                   loc=unlist(lapply(z,loc)),lon=unlist(lapply(z,lon)),
                   lat=unlist(lapply(z,lat)),alt=unlist(lapply(z,alt)),
                   info='Ensemble mean ESD')
  invisible(zm)
}


propchange <- function(x,it0=c(1979,2013)) {
  z <- coredata(x)
  if (is.null(dim(z)))
      z <- 100*(z/mean(coredata(subset(x,it=it0)),na.rm=TRUE)) else
      z <- 100*t(t(z)/apply(coredata(subset(x,it=it0)),2,'mean',na.rm=TRUE))
  attributes(z) <- NULL
  z -> coredata(x)  
  x
}

## Triangulation of pressure measurements to estimate wind
TGW <- function(triangle,f=1.25e-4,rho=1.25,verbose=FALSE) {
  if (verbose) print("Get stations")
  stopifnot(is.station(triangle))
  if (verbose) print(loc(triangle))
    
  if (verbose) print(paste("Length of overlapping interval is ",length(index(triangle)),
                           " from",min(year(triangle)),"to",max(year(triangle))))
  lons <- lon(triangle); lats <- lat(triangle)
  scal <- rep(1,3)
  ## If the units are 'hPa', then need to convert to 'Pa'
  if (length(unit(triangle))==1) attr(triangle,'unit') <- rep(unit(triangle),3)
  scal[is.element(unit(triangle),'hPa')] <- 100
  if (verbose) print(rbind(lons,lats))

  p1 <- scal[1]*coredata(triangle[,1]); p2 <- scal[2]*coredata(triangle[,2]);
  p3 <- scal[3]*coredata(triangle[,3])
  
  ## quality check: only accept values between between 800hPa and 1200hPa
  p1[p1 < 80000] <- NA; p1[p1 > 120000] <- NA; 
  p2[p2 < 80000] <- NA; p2[p2 > 120000] <- NA; 
  p3[p3 < 80000] <- NA; p3[p3 > 120000] <- NA; 
  x2 <- distAB(lons[2],lats[1],lons[1],lats[1])
  y2 <- distAB(lons[1],lats[2],lons[1],lats[1])
  x3 <- distAB(lons[3],lats[1],lons[1],lats[1])
  y3 <- distAB(lons[1],lats[3],lons[1],lats[1])
  if (verbose) print(paste("x2,y2,x3,y3=",x2,y2,x3,y3))
  
  a <- ( (p3 - p1) - y3*(p2 - p1)/y2 ) / (x3 - x2 * y3/y2)
  b <- (p2 - p1 - a*x2)/y2
  c <- p1
  if (verbose) print(c(length(p1),length(p2),length(p3),length(a),length(b),length(c)))  
  ug <- -b/(f*rho); attr(ug,"unit") <- "m/s"
  vg <- a/(f*rho); attr(vg,"unit") <- "m/s"

  ## Also flag values exceeding 50m/s (180km/h) as bad.
  ug[abs(ug) > 50] <- NA
  vg[abs(vg) > 50] <- NA
  
  if (verbose) print('station object')
  wind <- zoo(cbind(ug,vg),order.by=index(triangle))
  wind <- as.station(wind,loc=rep(paste(loc(triangle),collapse='-'),2),
                     lon=rep(mean(lons),2),lat=rep(mean(lats),2),
                     param=c('u','v'),unit=rep('m/s',2),
                     longname=c('zonal geostrophic wind','meridional geostrophic wind'),
                     info="Derived from triangular geostropic method",
                     ref="Alexandersson et al. (1998), Glob. Atm. and Oce. Sys., vol 6, pp. 97-120")
  attr(wind,'history') <- history.stamp(triangle)
  invisible(wind)
}

geostrophicwind.station <- function(x,f=1.25e-4,rho=1.25,verbose=FALSE) {
  ## Estimates the geostrophic wind from mean sea-level pressure from stations
  n <- length(loc(x))
  ## Estimate the different combinations of 3 that is possible from the provided group of
  ## stations
  cn <- combn(1:n,3)
  d <- dim(cn)
  print(n,'stations gives ',paste(d[2],'of three')
  for (i in 1:d[2]) {
    wind <- TGW(subset(x,is=cn[,i]))
    if (i==1) Wind <- wind else Wind <- combine(Wind,wind)
  }   
  invisible(Wind)      
}

geostrophicwind.field <- function(x,f=1.25e-4,rho=1.25,verbose=FALSE) {
  ## Estimates the geostrophic wind from mean sea-level pressure field
  if (verbose) print('geostrophicwind')
  stopifnot(is.field(x))
  if (sum(is.element(varid(x),c('slp','psl'))==0))
    warning(paste('geostrophicwind: param=',varid(x)))
  if (sum(is.element(unit(slp),'millibars','hPa')>0) {
    x <- 100*x
    attr(x,'unit') <- 'Pa'
  dpdx <- dX(x,verbose=verbose)
  dpdy <- dY(x,verbose=verbose)
  v <- 1/(f*rho)*dpdx$dZ
  u <- -1/(f*rho)*dpdy$dZ
  ws <- sqrt(u^2+v^2)
  class(ws) <- class(v)
  ws <- attrcp(v,ws)
  attr(u,'variable') <- 'u'
  attr(u,'unit') <- 'm/s'
  attr(u,'longname') <- 'zonal geostrophic wind'
  attr(v,'variable') <- 'u'
  attr(v,'unit') <- 'm/s'
  attr(v,'longname') <- 'meridional geostrophic wind'
  attr(ws,'variable') <- 'windspeed'
  attr(ws,'unit') <- 'm/s'
  attr(u,'longname') <- 'geostrophic wind speed'
  attr(u,'history') <- history.stamp(x)
  attr(v,'history') <- history.stamp(x)
  attr(ws,'history') <- history.stamp(x)
      
  invisible(list(u=u,v=v,ws=ws))
}
