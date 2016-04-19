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

