pieslice <- function(theta1,theta2,r=1,
                      col="grey",density=NULL,lwd=1,border=NULL) {
  s <- seq(theta1,theta2,by=1)
  if (is.null(border)) border <- col
  x <- c(0,cos(0.5*pi - pi*s/180),0); y <- c(0,sin(0.5*pi - pi*s/180),0)
  polygon(r*x,r*y,col=col,density=density,border=border)
}

## This code is based on an old version from clim.pact - can be rewritten to enhance efficiency.

windrose <- function(x,saw=10,max.scale=NULL,
            cols=c("grey90","yellow","green","red","blue","darkgreen",
                   "darkred","magenta","black"),param=c("u","v"),
                    simple=TRUE,verbose=FALSE) {
  if (verbose) print('windrose')
  ## Extract the zonal and merional components, stored as if they were different stations
  ## Make sure to extract matching records.
  u <- subset(x,is=is.element(varid(x),param[1]))
  v <- subset(x,is=is.element(varid(x),param[2]))
  ulonlat <- paste(lon(u),lat(u))
  vlonlat <- paste(lon(v),lat(v))
  i1 <- is.element(ulonlat,vlonlat)
  i2 <- is.element(vlonlat,ulonlat)
  u <- subset(u,is=i1); v <- subset(v,is=i2)
  ## If there are many stations, then run windrose recursively - once for each location
  if (!is.null(dim(u))) {
    if (verbose) print('Multiple locations')
    results <- list()
    for (is in 1:length(dim(u)[2])) {
      dev.new()
      uv <- combine.stations(subset(u,is=is),subset(v,is=is))
      ff <- windrose(uv,saw=saw,max.scale=max.scale,cols=cols,simple=simple,verbose=verbose)
      eval(parse(test=paste('results$ffdd.',is,' <- ff',sep='')))
    }
    return(results)
  }

  if (length(loc(u))!=1) {
    print(paste('Number of locations is ',length(loc(u))))
    print(loc(u))
  }
  
  ff <- sqrt(u^2 + v^2)
  dd <- 180/pi*atan2(u,v)
  ii <- is.finite(ff) & is.finite(dd)
  
  if (sum(ii)<100) {
    valid.stations <- stnr(param=param)
    print(paste("Too little (",sum(ii),") valid data was available"))
  }
  
  ff <- round(ff[ii]/5)*5; dd <- round(dd[ii]/saw)*saw
  if (sum(is.element(dd,360)>0) & sum(is.element(dd,0)>0)) {
    dd[is.element(dd,360)] <- 0
   }
  par(col.axis="white")
  plot(c(-1,1),c(-1,1),type="n",
       main=paste(loc(u),"wind rose; N=",sum(ii)),
       sub=paste("Lon=",round(lon(u),3),"E, lat=",round(lat(u),3),sep=""),xlab="S",ylab="W")

  for (ix in seq(-15,345,by=60)) {
    pieslice(ix,ix+30)
  }

  sectors <- as.numeric(rownames(table(dd)))
  speeds <-  as.numeric(rownames(table(ff)))
  N <- length(sectors); M <- length(speeds)
  categories <- matrix(rep(0,M*N),N,M)
  
  for (i in 1:N) {
    iv <- is.element(dd,sectors[i])
    tab <- table(ff[iv])
    rn <- as.numeric(rownames(tab))
    categories[i,is.element(speeds,rn)] <- as.numeric(tab)
  }

  if (is.null(max.scale)) maxr <- max(rowSums(categories)) else
                          maxr <- max.scale*sum(c(categories))/100
  nn <- min( round(100*maxr/(sum(c(categories)*5))*5),7 )
  
  ii <- 1
  for (ix in sectors) {
    for (il in seq(M,1,by=-1)) {
      r <- sum(categories[ii,1:il])/maxr
      if (sum(categories[ii,il])>0) pieslice(ix,ix+saw,col=cols[il],r=r)
    }
    ii <- ii+1
  }

  s <- seq(-2*pi,2*pi,length=360)
  lines(cos(s),sin(s))

  for (i in 1:nn) {
    iy <- i/nn
    lines(iy*cos(s),iy*sin(s),lty=3)
    portion <- maxr/sum(c(categories))*iy
    text(0,0-iy,round(portion*100),cex=0.5)
  }
  mtext("E",4); mtext("N",3)

  par(xpd=TRUE)
  legend(-1.2,-0.5,paste(rownames(table(ff)),"m/s"),
         col=cols,lty=1,lwd=3,cex=0.7,bg="grey95")

  if (verbose) print(table(ff,dd))
  par(col.axis="black")

  ivaldata <- is.finite(ff) & is.finite(dd)
  mtext(paste(start(x),end(x),sep="-"),1,adj=1,col="grey")
    
  if (!simple) {
    par0 <- par()
    par(fig=c(0.75,0.95,0.75,0.87),new=TRUE,mar=rep(0.5,4),
        xaxt="n",yaxt="n",cex.axis=0.5,cex.main=0.5)
    breaks <- seq(0,ceiling(max(ff,na.rm=TRUE))+5,by=5)
    ff[ff<0]<-0
    h <- hist(ff,breaks=breaks)
    par(fig=c(0.1,0.3,0.75,0.87),new=TRUE,mar=rep(0.5,4),
        xaxt="n",yaxt="n",cex.axis=0.5,cex.main=0.5)
    h <- hist(dd)
    par(par0)
  }
  invisible(list(ff=ff,dd=dd,categories=categories,
                 sectors=sectors))

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
  if (verbose) {
    print('Only accept pressures in the range 80000-120000Pa')
    print(paste('Set',sum(p1 < 80000 | p1 > 120000),
                sum(p2 < 80000 | p2 > 120000), sum(p3 < 80000 | p3 > 120000),
                ' outside this range to NA'))
  }
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
  if (verbose) {
    print(c(length(p1),length(p2),length(p3),length(a),length(b),length(c)))
    print(summary(a)); print(summary(b)); print(summary(c)) 
  }
  ug <- -b/(f*rho); attr(ug,"unit") <- "m/s"
  vg <- a/(f*rho); attr(vg,"unit") <- "m/s"

  ## Also flag values exceeding 50m/s (180km/h) as bad.
  if (verbose) print(paste('Set to NA:',sum(abs(ug) > 50),'for u;',sum(abs(vg) > 50),'for v'))
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



geostrophicwind<-function(x,...) UseMethod("geostrophicwind")

geostrophicwind.station <- function(x,f=1.25e-4,rho=1.25,verbose=FALSE,nmax=1000) {
  ## Estimates the geostrophic wind from mean sea-level pressure from stations
  n <- length(loc(x))
  ## Estimate the different combinations of 3 that is possible from the provided group of
  ## stations
  cn <- combn(1:n,3)
  d <- dim(cn)
  print(paste(n,'stations gives ',d[2],'combinations of three.'))
  if (!is.null(nmax)) {
    ## For very many combination, take a random sample by default.
    ## Turn this feature off by setting nmax to NULL.
    if (nmax < d[2]) {
      print(paste('Taking a random subsample of',nmax,'combinations.',
                  'If all are wnted, set argument nmax =NULL'))
      ii <- sample(1:d[2],nmax)
      cn <- cn[,ii]; d <- dim(cn)
    }
  }
  pb <- txtProgressBar(style=3)
  for (i in 1:d[2]) {
    wind <- TGW(subset(x,is=cn[,i]))
    setTxtProgressBar(pb,i/d[2]) 
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
  if (sum(is.element(unit(x),c('millibars','hPa')))>0) {
    x <- 100*x
    attr(x,'unit') <- 'Pa'
  }
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
