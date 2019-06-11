#' Functions to process and analyse storm trajectories
#' 
#' @aliases anomaly.trajectory polyfit.trajectory season.trajectory 
#' count.trajectory param.trajectory sort.trajectory polyfit.trajectory
#' anomaly2trajectory fnlon trajectory2density
#' 
#' @param x A trajectory object
#' @param \dots Other arguments
#' @param type type of anomaly: 'first' gives you the spatial anomaly with
#' regards to the first time step of the trajectory and 'mean' centers the
#' trajectories around the mean longitude and latitude
#' @param param parameters to calculate anomaly of
#' @param verbose if TRUE print information about progress
#' @return A trajecory object
#' @author Kajsa M. Parding
#' @examples
#' 
#' # Load trajectory example data
#' data('imilast.M03')
#' # Calculate anomaly
#' x <- anomaly.trajectory(imilast.M03)
#' # Show maps of original trajectories and spatial anomalies
#' map(imilast.M03)
#' map(x)
#' # Transform trajectory anomalies back to regular trajectories
#' y <- anomaly2trajectory(x)
#' # Print longitudes of first trajectory
#' imilast.M03[1,1:12]
#' x[1,1:12]
#' y[1,12]
#' 
#' # Fit polynomial to trajectories
#' y <- polyfit.trajectory(x)
#' # Show coefficients of first trajectory
#' print(attr(y,"coefs")[,1])
#' # Plot original trajectory and polynomial fit
#' ilon <- colnames(x)=="lon"
#' ilat <- colnames(x)=="lat"
#' plot(x[1,ilon], x[1,ilat], col="black", pch=1)
#' lines(y[1,ilon], y[1,ilat], col="blue", lty=2)
#'
#' @export
season.trajectory <- function(x, format="character",verbose=FALSE) {
  if(verbose) print("season.trajectory")
  stopifnot(!missing(x), inherits(x,"trajectory"))
  mn <- month(x)
  mlist <- c(12,1:11)
  slist <- c(rep(1,3),rep(2,3),rep(3,3),rep(4,4))
  if(format=="character") {
    slist <- sapply(slist,
      function(x) switch(x, "1"="djf", "2"="mam", "3"="jja", "4"="son"))
  }
  sn <- sapply(mn,function(x) slist[mlist==x])
  invisible(sn)
}

#' @export
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

#' @export
sort.trajectory <- function(x, decreasing=FALSE, ...) {
  stopifnot(!missing(x), inherits(x,"trajectory"))
  if (any('sorted' %in% attr(x,'aspect'))) {
    invisible(x)
  } else {
    start <- x[,colnames(x)=='start']
    date <- strptime(start,format="%Y%m%d%H")
    while (any(duplicated(date))) {
      date[duplicated(date)] <- date[duplicated(date)]+60
    }
    y <- x[order(date, decreasing=decreasing),]
    y <- attrcp(x,y)
    attr(y,'aspect') <- c('sorted',attr(x,'aspect'))
    class(y) <- class(x)
    invisible(y)
  }
}

#' @export
polyfit.trajectory <- function(x,verbose=FALSE) {
  if(verbose) print("polyfit.trajectory")
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
  coefs <- apply(x,1,function(y) attr(polyfit(y[ilon],y[ilat]),'coefs'))
  z <- x
  z[,ilat] <- t(latfit)
  attr(z,'aspect') <- unique(c("pfit",attr(x,'aspect')))
  attr(z,'coefs') <- coefs
  return(z)
}

polyfit <- function(x,y=NULL) {
  if (all(is.null(y))) {
    y <- x
    x <- seq(1,length(x))
  }
  pfit <- lm(y ~ I(x) + I(x^2) + I(x^3) + I(x^4) + I(x^5))
  z <- predict(pfit)
  attr(z,'coefs') <- pfit$coef
  return(z)
}

#' @export
count.trajectory <- function(x,it=NULL,is=NULL,by='year',verbose=FALSE) {
  if(verbose) print("count.trajectory")
  y <- subset(x,it=it,is=is)
  if(is.null(attr(x,"calendar"))) calendar <- "gregorian" else calendar <- attr(x,"calendar")
  if (requireNamespace("PCICt", quietly = TRUE)) {
    t <- PCICt::as.PCICt(as.character(y[,colnames(y)=="start"]),format="%Y%m%d%H",cal=calendar)
    nrt <- PCICt::as.PCICt(as.character(range(year(t))*1E4+range(month(t))*1E2+1),
                           format="%Y%m%d",cal=attr(x,"calendar"))
  } else {
    t <- as.POSIXct(as.character(y[,colnames(y)=="start"]),format="%Y%m%d%H")
    nrt <- as.POSIXct(as.character(range(year(t))*1E4+range(month(t))*1E2+1),format="%Y%m%d")
  }
  if (by=='year') {
    fmt <- "%Y"
    if (inherits(y,'season')) {
      cls <- 'season'
      unit <- 'events/season'
    } else {
      cls <- 'annual'
      unit <- 'events/year'
    }
    n0 <- zoo(,seq(from = year(nrt[1]), to = year(nrt[2])))
  } else if (by %in% c('month','4seasons')) {
    fmt <- "%Y%m%d"
    t <- as.Date(paste(format(t,format="%Y-%m"),"01",sep="-"))
    nrt <- as.Date(format(nrt,format="%Y-%m-%d"))
    cls <- 'month'
    unit <- 'events/month'
    n0 <- zoo(,seq(from = nrt[1], to = nrt[2], by = "month"))
  } else if (by=='day') {
    fmt <- "%Y%m%d"
    cls <- 'day'
    unit <- 'events/day'
    n0 <- zoo(,seq(from = nrt[1], to = nrt[2], by = "day"))
  }
  d <- format(t,format=fmt)
  n <- table(d)
  if (by=='year') {
    dn <- as.integer(dimnames(n)$d)
  } else if(by=="day" & requireNamespace("PCICt", quietly = TRUE)) {
    dn <- PCICt::as.PCICt(dimnames(n)$d,format=fmt,cal=calendar)
  } else {
    dn <- as.Date(strptime(dimnames(n)$d,format=fmt))
  }
  nz <- zoo(n,order.by=dn)
  n <- merge(nz, n0)
  n[is.na(n)] <- 0
  n <- attrcp(y,n)
  n <- as.station(n,param='trajectory count',unit=unit,
                  calendar=calendar,longname='trajectory count')
  if (inherits(y,"season")) {
    ok <- sapply(month(index(n0)),function(x) x %in% unique(month(t)))
    n <- subset(n,it=ok)
  }
  if (by=='4seasons') n <- as.4seasons(n,FUN=sum)
  n <- subset(n,it=format(c(min(t),max(t)),format="%Y-%m-%d"))
  if(!is.null(attr(y,"longitude"))) attr(n,"longitude") <- attr(y,"longitude")
  if(!is.null(attr(y,"latitude"))) attr(n,"latitude") <- attr(y,"latitude")
  invisible(n)
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


spherical2cartesian <- function(lon,lat,n=10,a=6.378e06,verbose=TRUE) {
  if (verbose) print("spherical2cartesian")
  x <- a * cos( lat*pi/180 ) * cos( lon*pi/180 )
  y <- a * cos( lat*pi/180 ) * sin( lon*pi/180 )
  z <- a * sin( lat*pi/180 )
  invisible(cbind(x,y,z))
}

cartesian2spherical <- function(x,y,z,a=6.378e06,verbose=TRUE) {
  if (verbose) print("cartesian2spherical")
  lon <- atan2( y, x )*180/pi
  lat <- asin( z/sqrt( x^2 + y^2 + z^2 ))*180/pi
  invisible(cbind(lon,lat))
}

#' @export
anomaly.trajectory <- function(x,...,type='first',param=c('lon','lat'),
                                verbose=FALSE) {
  if(verbose) print("anomaly.trajectory")
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

  #if (any('lon' %in% param)) {
  #  i.lon <- which(colnames(x)=='lon')
  #  i.lat <- which(colnames(x)=='lat')
  #  if (type=='first') {
  #    p.anomaly <- apply(x,1,function(x) {
  #      xyz <- spherical2cartesian(x[i.lon],x[i.lat])
  #                           })
  #    m[param==p] <- list(xi[1])
  #  } else if (type=='mean') {
  #    m[param==p] <- list(mean(x[,i]))
  #    p.anomaly <- apply(x,1,function(x) xi-m[param==p][[1]])
  #  }
  #}
  for (p in param) {
    i <- which(colnames(x)==p)
    xi <- x[,i]
    if (p=='lon') xi <- t(apply(x[,i],1,lontrack))
    if (type=='first') {
      m[param==p] <- list(apply(xi,1,function(x) x[1]))#list(xi[1])
      p.anomaly <- apply(xi,1,function(x) x-x[1])
    } else if (type=='mean') {
      m[param==p] <- list(mean(xi,na.rm=TRUE))
      p.anomaly <- apply(xi,1,function(x) x-m[param==p][[1]])
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

#' @export
trajectory2density <- function(x,it=NULL,is=NULL,dx=2,dy=2,radius=5E5,verbose=FALSE) {
  if(verbose) print("trajectory2density")
  y <- subset(x,it=it,is=is)
  if(!is.null(dim(y))) {
    A <- apply(y,1,function(x) trackdensity(x[colnames(y)=='lon'],
                                            x[colnames(y)=='lat'],
                                            dx=dx,dy=dy,radius=radius))
    lon <- unlist(lapply(A,function(x) factor2numeric(x$lon)))
    lat <- unlist(lapply(A,function(x) factor2numeric(x$lat)))
    hits <- as.data.frame(table(lon,lat))
    names(hits)[names(hits)=="Freq"] <- "density"
  } else {
    A <- trackdensity(y[names(y)=='lon'],y[names(y)=='lat'],dx=dx,dy=dy,radius=radius)
    lon <- factor2numeric(A$lon)
    lat <- factor2numeric(A$lat)
    hits <- as.data.frame(table(lon,lat))
    names(hits)[names(hits)=="Freq"] <- "density"
  }
  invisible(hits)
}

#' @export
anomaly2trajectory <- function(x,verbose=FALSE) {
  if(verbose) print("anomaly2trajectory")
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

#' @export
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

approxlon <- function(lon,n=10) {
  x <- approx(lontrack(lon),n=n)
  x$y <- lon2dateline(x$y)
  return(x)
}
