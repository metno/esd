# K Parding, 29.05.2015

CCI <- function(Z,m=14,nsim=NULL,it=NULL,is=NULL,cyclones=TRUE,
                label=NULL,fname="cyclones.rda",mindistance=1E6,
                accuracy=NULL,pmax=NULL,rmin=1E4,rmax=2E6,dpmin=1,
                lplot=FALSE,verbose=FALSE) {
  stopifnot(inherits(Z,'field'))
  Z <- subset(Z,it=it,is=is)
  if (any(longitude(Z)>180)) Z <- g2dl(Z,greenwich=FALSE)

  ## Rearrange time index
  t <- as.numeric(strftime(index(Z),format="%Y%m%d%H%M"))
    
  ## Calculate first and second derivative
  if(verbose) print("Calculate first and second derivative")
  resx <- dX(Z,m=m,accuracy=accuracy,verbose=verbose)
  resy <- dY(Z,m=m,accuracy=accuracy,verbose=verbose)

  ## Spatial resolution
  dx <- abs(diff(resx$lon))[1]
  dy <- abs(diff(resx$lat))[1]
  
  ## Search for zero crossings of the first derivative:
  ## Reorganize and reshape
  if(verbose) print("Reshape arrays")
  nx <- length(resx$lon); ny <- length(resx$lat); nt <- dim(Z)[1]
  lonXY <- rep(0.5*(resx$lon[2:nx]+resx$lon[1:(nx-1)]),ny-1)
  dim(lonXY) <- c(nx-1,ny-1)
  latXY <- rep(0.5*(resx$lat[2:ny]+resx$lat[1:(ny-1)]),nx-1)
  dim(latXY) <- c(ny-1,nx-1); latXY <- t(latXY)
  Zx <- as.matrix(resx$Z.fit); dim(Zx) <- c(nt,nx,ny)
  Zy <- as.matrix(resy$Z.fit); dim(Zy) <- c(nt,nx,ny)
  dslpdy <- as.matrix(resy$dZ); dim(dslpdy) <- c(nt,nx,ny)
  dslpdx <- as.matrix(resx$dZ); dim(dslpdx) <- c(nt,nx,ny)
  dslpdx2 <- as.matrix(resx$dZ2); dim(dslpdx2) <- c(nt,nx,ny)
  dslpdy2 <- as.matrix(resy$dZ2); dim(dslpdy2) <- c(nt,nx,ny)
  px <- 0.25*(Zx[,1:(nx-1),2:ny] + Zx[,2:nx,2:ny] +
              Zx[,1:(nx-1),1:(ny-1)] + Zx[,2:nx,1:(ny-1)])
  py <- 0.25*(Zy[,1:(nx-1),2:ny] + Zy[,2:nx,2:ny] +
              Zy[,1:(nx-1),1:(ny-1)] + Zy[,2:nx,1:(ny-1)])

  ## Clear temporary objects from working memory
  rm("Zx","Zy","resx","resy"); gc(reset=TRUE)
   
  ## Search for zero crossing in y-direction
  if(verbose) print("Find zero crossing of first derivative in y-direction")
  P.lowy <- rep(0,nt*(nx-1)*(ny-1)); dim(P.lowy) <- c(nt,nx-1,ny-1) 
  dy11 <- 0.5*(dslpdy[,2:nx,2:ny]+dslpdy[,1:(nx-1),2:ny])
  dy12 <- 0.5*(dslpdy[,2:nx,1:(ny-1)]+dslpdy[,1:(nx-1),1:(ny-1)])
  dy21 <- 0.5*(dslpdy2[,2:nx,2:ny]+dslpdy2[,1:(nx-1),2:ny])
  dy22 <- 0.5*(dslpdy2[,2:nx,1:(ny-1)]+dslpdy2[,1:(nx-1),1:(ny-1)])
  if (cyclones) { i.low <- (dy11*dy12 < 0) & (dy21+dy22 > 0) &
              is.finite(dy11+dy12+dy21+dy22)
  } else { i.low <- (dy11*dy12 < 0) & (dy21+dy22 < 0) &
              is.finite(dy11+dy12+dy21+dy22)
  }
  P.lowy[i.low] <- 1
  DY <- 0.5*(dy11+dy12)
  DY2 <- 0.5*(dy21+dy22)

  ## Zero crossing in x-direction
  if(verbose) print("Find zero crossing of first derivative in x-direction")
  P.lowx <- rep(0,nt*(nx-1)*(ny-1)); dim(P.lowx) <- c(nt,nx-1,ny-1) 
  dx11 <- 0.5*(dslpdx[,2:nx,2:ny] + dslpdx[,2:nx,1:(ny-1)])
  dx12 <- 0.5*(dslpdx[,1:(nx-1),2:ny] + dslpdx[,1:(nx-1),1:(ny-1)])
  dx21 <- 0.5*(dslpdx2[,2:nx,2:ny] + dslpdx2[,2:nx,1:(ny-1)])
  dx22 <- 0.5*(dslpdx2[,1:(nx-1),2:ny] + dslpdx2[,1:(nx-1),1:(ny-1)])
  if (cyclones) { i.low <- (dx11*dx12 < 0) & (dx21+dx22 > 0) &
              is.finite(dx11 + dx12 + dx21 + dx22)
  } else { i.low <- (dx11*dx12 < 0) & (dx21+dx22 < 0) &
              is.finite(dx11 + dx12 + dx21 + dx22)
  }
  P.lowx[i.low] <- 1
  DX <- 0.5*(dx11+dx12)
  DX2 <- 0.5*(dx21+dx22)

  ## Clear temporary objects from working memory
  rm("dx11","dx12","dx21","dx22","dy11","dy12","dy21","dy22",
     "dslpdx","dslpdx2","dslpdy","dslpdy2","i.low"); gc(reset=TRUE)
  
  ## Find zero crossings in both directions
  ## Plowx & P.lowy are matrices with 0's and 1's.
  lows1 <- (P.lowy & P.lowx)
  qf <- matrix(rep(0,length(P.lowx)),dim(P.lowx))
  qf[lows1] <- 1
  
  ## Find cyclones that were missed in the first search
  ## using widened masks of zero crossings
  if(verbose) print("Find extra cyclones using widened masks")
  ## Mask the cylcones already selected
  P.lowx2 <- P.lowx; P.lowy2 <- P.lowy
  P.lowx2[lows1] <- 0; P.lowy2[lows1] <- 0
  ## Widen mask in x-direction
  P.lowx2[,2:(nx-1),] <- P.lowx2[,1:(nx-2),] | P.lowx2[,2:(nx-1),]
  P.lowx2[,1:(nx-2),] <- P.lowx2[,1:(nx-2),] | P.lowx2[,2:(nx-1),]
  ## Widen mask in y-direction
  P.lowy2[,,2:(ny-1)] <- P.lowy2[,,1:(ny-2)] | P.lowy2[,,2:(ny-1)]
  P.lowy2[,,1:(ny-2)] <- P.lowy2[,,1:(ny-2)] | P.lowy2[,,2:(ny-1)]
  ## Find new zero crossings
  lows2 <- (P.lowy2 & P.lowx2)
  qf[lows2] <- 2
 
  ## Lat, lon, and dates of cyclones
  lon<-rep(lonXY,nt); dim(lon)<-c(nx-1,ny-1,nt); lon<-aperm(lon,c(3,1,2)) 
  lat<-rep(latXY,nt); dim(lat)<-c(nx-1,ny-1,nt); lat<-aperm(lat,c(3,1,2))
  date<-rep(t,(nx-1)*(ny-1)); dim(date)<-c(nt,nx-1,ny-1)
  ## Cyclones found with ordinary method
  lon1 <- lon[lows1]; lat1<-lat[lows1]; date1 <- date[lows1]
  pcent1 <- 0.5*(px[lows1]+py[lows1])
  strength1 <- rank(pcent1)
  if (!cyclones) strength1 <- rank(-pcent1)
  ## Cyclones found after widening zero crossing masks
  lon2 <- lon[lows2]; lat2<-lat[lows2]; date2 <- date[lows2]
  pcent2 <- 0.5*(px[lows2]+py[lows2])
  strength2 <- rank(pcent2)
  if (!cyclones) strength2 <- rank(-pcent2)

  ## Keep only cyclones that are deeper than pmax
  del1 <- rep(TRUE,length(date1))
  del2 <- rep(TRUE,length(date2))
  if(!is.null(pmax)) {
    if(is.function(pmax)) {
      pmax <- pmax(Z,na.rm=T)
    } else if(is.character(pmax)) {
      pmax <- eval(parse(text=paste(pmax,"(Z,na.rm=T)",sep="")))
    }
    ok1 <- pcent1 < pmax
    if(!cyclones) ok1 <- pcent1 > pmax
    del1[!ok1] <- FALSE
    ok2 <- pcent2 > pmax
    if(!cyclones) ok2 <- pcent2 > pmax
    del2[!ok2] <- FALSE
  }
  lows1[lows1] <- del1
  lon1 <- lon1[del1]; lat1 <- lat1[del1]; date1 <- date1[del1]
  pcent1 <- pcent1[del1]; strength1 <- strength1[del1]
  lows2[lows2] <- del2
  lon2 <- lon2[del2]; lat2 <- lat2[del2]; date2 <- date2[del2]
  pcent2 <- pcent2[del2]; strength2 <- strength2[del2]

  ## Clear temporary objects from working memory
  rm("ok1","ok2","del1","del2"); gc(reset=TRUE)
  
  ## Remove secondary cyclones near a deeper one (same cyclonic system):
  if(verbose) print("Remove secondary cyclones")
  del1 <- rep(TRUE,length(date1))
  del2 <- rep(TRUE,length(date2))
  for (d in t) {
    if (inherits(index(Z),'Date')) d <- as.Date(d)
    ## Remove secondary cyclones identified with the standard method,
    ## located close (<1000km) to other stronger cyclones
    i1 <- which(date1==d)
    if (length(i1)>1) {
      distance <- apply(cbind(lon1[i1],lat1[i1]),1,
       function(x) suppressWarnings(distAB(x[1],x[2],lon1[i1],lat1[i1])))
      diag(distance) <- NA; distance[lower.tri(distance)] <- NA
      del.i1 <- which(distance<mindistance,arr.ind=TRUE)
      if(any(del.i1)) {
        col.del <- rep(1,length(del.i1)/2)
        if (is.null(dim(del.i1))) {
           s1 <- strength1[i1][,del.i1[1]]
           s2 <- strength1[i1][,del.i1[2]]
           col.del[s1<=s2] <- 2
           del.i1 <- unique(del.i1[col.del])
        } else {
           s1 <- strength1[i1][del.i1[,1]]
           s2 <- strength1[i1][del.i1[,2]]
           col.del[s1<=s2] <- 2
           del.i1 <- unique(del.i1[cbind(seq(1,dim(del.i1)[1]),col.del)])
        }
        del1[i1[del.i1]] <- FALSE
      }
      i1 <- i1[!1:length(i1) %in% del.i1]
    }
    ## Remove secondary cyclones identified with the widened masks,
    ## requiring longer distance between cyclones (2000 km)
    i2 <- which(date2==d)
    if(any(i2) & any(i1)) {
      distance <- apply(cbind(lon2[i2],lat2[i2]),1,
       function(x) suppressWarnings(distAB(x[1],x[2],lon1[i1],lat1[i1])))
      del.i2 <- unique(which(distance<2E6,arr.ind=TRUE)[,2])
      del2[i2[del.i2]] <- FALSE
      i2 <- i2[!1:length(i2) %in% del.i2]
    }
    if(length(i2)>1) {
      distance <- apply(cbind(lon2[i2],lat2[i2]),1,
       function(x) suppressWarnings(distAB(x[1],x[2],lon2[i2],lat2[i2])))
      diag(distance) <- NA; distance[lower.tri(distance)] <- NA
      del.i2 <- which(distance<2E6,arr.ind=TRUE)
      if(any(del.i2)) {
        col.del <- rep(1,length(del.i2)/2)
        if (is.null(dim(del.i2))) {
          s1 <- strength2[i2][del.i2[1]]
          s2 <- strength2[i2][del.i2[2]]
          col.del[s1<=s2] <- 2
          del.i2 <- unique(del.i2[col.del])
        } else {
          s1 <- strength2[i2][del.i2[,1]]
          s2 <- strength2[i2][del.i2[,2]]
          col.del[s1<=s2] <- 2
          del.i2 <- unique(del.i2[cbind(seq(1,dim(del.i2)[1]),col.del)])
        }
        del2[i2[del.i2]] <- FALSE
      }
    }
  }
  lows1[lows1] <- del1
  lows2[lows2] <- del2
 
  ## Add the two groups of cyclones together,
  ## keep track of which is which with the quality flag qf
  lows <- lows1 | lows2
  lon <- lon[lows]
  lat <- lat[lows]
  date <- date[lows]
  pcent <- 0.5*(px[lows]+py[lows])
  qf <- qf[lows]

  ## Clear temporary objects from working memory
  rm("lows1","lows2","lon1","lon2","lat1","lat2",
     "date1","date2","strength1","strength2",
     "pcent1","pcent2","del1","del2"); gc(reset=TRUE)

  ## Put cyclones in order of date
  i <- order(date)
  lon <- lon[i]
  lat <- lat[i]
  date <- date[i]
  pcent <- pcent[i]
  qf <- qf[i]
  strength <- rank(pcent)
  if (!cyclones) strength <- rank(-pcent)
 
  ## Keep only the nsim strongest cyclones each time step
  if(!is.null(nsim)) {
    s <- sapply(date,function(d) rank(strength[date==d]))
    lon <- lon[s<nsim]
    lat <- lat[s<nsim]
    date <- date[s<nsim]
    pcent <- pcent[s<nsim]
    qf <- qf[s<nsim]
  }

  ## Pressure gradient
  if(verbose) print("Pressure gradient")
  rho <- 1.2922
  dpsl <- sqrt(DX^2+DY^2)
  if (attr(Z,"unit")=="hPa") dpsl <- dpsl*100

  # Find points of inflexion (2nd derivative==0) to estimate
  # the storm radius and maximum speed and pressure gradient
  if(verbose) print("Find points of inflexion")
  NX <- dim(lonXY)[1]; NY <- dim(latXY)[2]
  lonXX <- rep(0.5*(lonXY[2:NX,2:NY]+lonXY[1:(NX-1),2:NY]),nt)
  dim(lonXX) <- c(NX-1,NY-1,nt); lonXX <- aperm(lonXX,c(3,1,2))
  latXX <- rep(0.5*(latXY[2:NX,2:NY]+latXY[2:NX,1:(NY-1)]),nt)
  dim(latXX) <- c(NX-1,NY-1,nt); latXX <- aperm(latXX,c(3,1,2))
  dateXX <- rep(t,(NX-1)*(NY-1)); dim(dateXX) <- c(nt,NX-1,NY-1)
  radius <- rep(NA,length(date))
  depth <- rep(NA,length(date))
  max.dslp <- rep(NA,length(date))
  max.speed <- rep(NA,length(date))
  max.vg <- rep(NA,length(date))
  for (i in seq(1,length(date))) {
    ilon <- abs(lonXX[1,,1]-lon[i])<dx
    ilat <- abs(latXX[1,1,]-lat[i])<dy
    inflx <- DX2[date[i]==t,2:NX,latXY[1,]==lat[i]]*
        DX2[date[i]==t,1:(NX-1),latXY[1,]==lat[i]]
    infly <- DY2[date[i]==t,lonXY[,1]==lon[i],2:NY]*
        DY2[date[i]==t,lonXY[,1]==lon[i],1:(NY-1)]
    lon.infl <- lonXY[inflx<0,1]
    lat.infl <- latXY[1,infly<0]
    dlon <- lon.infl-lon[i]
    dlat <- lat.infl-lat[i]
    ilon <- c()
    ilat <- c()
    if (any(dlat>0)) {
      ilon <- c(ilon,which(lonXY[,1]==lon[i]))
      ilat <- c(ilat,which(latXY[1,]==lat.infl[dlat==min(dlat[dlat>0])]))
    }
    if (any(dlat<0)) {
      ilon <- c(ilon,which(lonXY[,1]==lon[i]))
      ilat <- c(ilat,which(latXY[1,]==lat.infl[dlat==max(dlat[dlat<0])]))
    }
    if (any(dlon>0)) {
      ilon <- c(ilon,which(lonXY[,1]==lon.infl[dlon==min(dlon[dlon>0])]))
      ilat <- c(ilat,which(latXY[1,]==lat[i]))
    }
    if (any(dlon<0)) {
      ilon <- c(ilon,which(lonXY[,1]==lon.infl[dlon==max(dlon[dlon<0])]))
      ilat <- c(ilat,which(latXY[1,]==lat[i]))
    }
    dpi <- mapply(function(i1,i2) dpsl[t==date[i],i1,i2],ilon,ilat)
    ri <- distAB(lon[i],lat[i],lonXY[ilon,1],latXY[1,ilat])
    pxi <- mapply(function(i1,i2) px[t==date[i],i1,i2],ilon,ilat)
    pyi <- mapply(function(i1,i2) py[t==date[i],i1,i2],ilon,ilat)
    di <- abs(0.5*(pxi+pyi)-pcent[i])
    ## Geostrophic wind
    fi <- 2*7.29212*1E-5*sin(pi*latXY[1,ilat]/180)
    vg <- dpi/(fi*rho)
    ## Gradient wind
    v.grad <- -0.5*fi*pi*ri*(1 - sqrt(1 + 4*vg/(fi*ri)))
    radius[i] <- mean(ri)#ri[which.max(dpi)]
    max.dslp[i] <- mean(dpi)#dpi[which.max(dpi)]
    max.speed[i] <- mean(v.grad)#v.grad[which.max(dpi)]
    max.vg[i] <- mean(vg)#vg[which.max(dpi)]
    depth[i] <- mean(di)
  
    ## Plot examples of cyclones
    if (lplot & length(ilon)==4) {
      if(verbose) print("plot examples of cyclone identification")
      data(geoborders)
      pxi <- px[date[i]==t,,];  pyi <- py[date[i]==t,,]
      lon.i <- lonXY[ilon,1]; lat.i <- latXY[1,ilat] 
      xi <- lonXY[,1]; yi <- latXY[1,]; zi <- pxi
      if (!(all(diff(xi)>0))) {xi <- rev(xi); zi <- apply(zi,2,rev)}
      if (!(all(diff(yi)>0))) {yi <- rev(yi); zi <- t(apply(t(zi),2,rev))}
      dev.new()
      plot(lonXY[,1],pxi[,latXY[1,]==lat[i]],lty=1,type="l",main=date[i],
           xlab="lon",ylab="slp (hPa)")
      points(lon[i],pxi[lonXY[,1]==lon[i],latXY[1,]==lat[i]],col="blue",pch=19)
      points(lonXY[inflx<0,1],pxi[inflx<0,latXY[1,]==lat[i]],col="red",pch=1)
      points(lon.i[lat.i==lat[i]],pxi[ilon[lat.i==lat[i]],latXY[1,]==lat[i]],
             col="red",pch=20)
      dev.copy(jpeg,"cyclones.lon.jpg"); dev.off()
      dev.new()
      plot(latXY[1,],pyi[lonXY[,1]==lon[i],],lty=1,type="l",main=date[i],
           xlab="lat",ylab="slp (hPa)")
      points(lat[i],pyi[lonXY[,1]==lon[i],latXY[1,]==lat[i]],col="blue",pch=19)
      points(latXY[1,infly<0],pyi[lonXY[,1]==lon[i],infly<0],col="red",pch=1)
      points(lat.i[lon.i==lon[i]],pyi[lonXY[,1]==lon[i],ilat[lon.i==lon[i]]],
             col="red",pch=20)
      dev.copy(jpeg,"cyclones.lat.jpg"); dev.off()
      dev.new()
      image(xi,yi,zi,main=date[i],col=colscal(col="t2m",n=12,rev=FALSE),
            xlab="lon",ylab="lat")
      contour(xi,yi,zi,add=TRUE,col='Grey40',lty=1,zlim=c(950,1010),nlevels=6)
      contour(xi,yi,zi,add=TRUE,col='Grey40',lty=2,zlim=c(1020,1060),nlevels=5)
      lines(geoborders)
      a <- which(P.lowx[t==date[i],,]==1,arr.ind=TRUE)
      b <- which(P.lowy[t==date[i],,]==1,arr.ind=TRUE)
      lon.a <- mapply(function(i1,i2) lonXY[i1,i2],a[,1],a[,2])
      lat.a <- mapply(function(i1,i2) latXY[i1,i2],a[,1],a[,2])
      lon.b <- mapply(function(i1,i2) lonXY[i1,i2],b[,1],b[,2])
      lat.b <- mapply(function(i1,i2) latXY[i1,i2],b[,1],b[,2])
      points(lon.a,lat.a,col="black",pch=1,cex=0.5,lwd=0.5)
      points(lon.b,lat.b,col="black",pch="|",cex=0.5,lwd=0.5)
      j <- date==date[i]
      points(lon[j][qf[j]==1],lat[j][qf[j]==1],pch=21,lwd=2,
             bg="white",col="black",cex=2)
      points(lon[j][qf[j]==2],lat[j][qf[j]==2],pch=21,lwd=2,
             bg="white",col="black",cex=1)
      points(lon[i],lat[i],pch=4,lwd=2,col="black",cex=1)
      dev.copy(jpeg,"cyclones.map.jpg"); dev.off()
      lplot <- FALSE; rm("geoborders")
    }
  }
  ## Remove temporary variables and release the memory:
  rm('lonXX','latXX','dateXX','lonXY','latXY','inflx','infly'); gc(reset=TRUE)
  
  ## Keep only cyclones with radius within the range (rmin,rmax)
  ## and depth stronger than dpmin
  ok <- rep(TRUE,length(date))
  if(!is.null(rmin)) ok <- ok & radius>rmin
  if(!is.null(rmax)) ok <- ok & radius<rmax
  if(!is.null(dpmin)) ok <- ok & depth>dpmin
  lon <- lon[ok]
  lat <- lat[ok]
  date <- date[ok]
  pcent <- pcent[ok]
  qf <- qf[ok]
  max.dslp <- max.dslp[ok]
  max.speed <- max.speed[ok]
  radius <- radius[ok]
  depth <- depth[ok]
  strength <- rank(pcent)
  if (!cyclones) strength <- rank(-pcent)   
 
  ## Check units, Pa -> hPa
  max.dslp <- max.dslp*1E-2
  if (attr(Z,"unit")=='Pa') {
      pcent <- pcent*1E-2
      depth <- depth*1E-2
  }

  ## Arrange results
  date <- strptime(date,"%Y%m%d%H%M")
  dd <- as.numeric(strftime(date,"%Y%m%d"))
  hh <- as.numeric(strftime(date,"%H"))
  units <- c("date","hour CET","degrees","degrees","hPa","hPa/m",
             "m/s","km","hPa","quality")
  if (cyclones) {
    longname <- "low pressure systems identified with CCI method"
    variable <- "cyclones"
  } else {
    longname <- "high-pressure systems identified with CCI method"
    variable <- "anti-cyclones"
  }
  qflabel <- paste("1:",variable,"identified in cross-sections",
     "between EW and NS pressure gradient zero crossings;",
     "2: less accurate",variable,"identification from widened",
     "EW and NS pressure gradient zero crossings")
  X <- data.frame(date=dd,time=hh,lon=lon,lat=lat,pcent=pcent,
         max.dslp=max.dslp,max.speed=max.speed,
         radius=radius,depth=depth,quality=qf)
  X <- as.events(X,label=label,dx=dx,dy=dy,units=units,longname=longname,
         variable=variable,qflabel=qflabel,src=attr(Z,"source"),
         file=attr(Z,"file"),version="CCI in esd v1.0 (after August 25, 2015)")
  if(!is.null(fname)) save(file=fname,X)
  invisible(X)
}


 


#library(esd)
### MONTHLY SLP DATA:
##slp <- slp.ERAINT()
##cyclones <- CCI(slp,is=list(lon=c(-180,180),lat=c(0,90)),it='djf')
### 6-HOURLY SLP DATA:
## ncfile <- "/home/kajsamp/data/ecmwf/ERAinterim_slp_1979.nc"
## ncid <- open.ncdf(ncfile)
## lon <- get.var.ncdf(ncid,"longitude")
## lat <- get.var.ncdf(ncid,"latitude")
## time <- get.var.ncdf(ncid,"time")
## tunits <- att.get.ncdf(ncid, "time", "units")
## slp <- get.var.ncdf(ncid,"msl")
## close.ncdf(ncid)
## time <- as.POSIXlt(time*60*60,origin="1900-01-01")
## nt <- length(time); nlon <- length(lon); nlat <- length(lat)
## slp <- aperm(slp,c(3,1,2)); dim(slp) <- c(nt,nlon*nlat)
## slp <- zoo(slp,order.by=time)
## slp <- as.field(slp,lon=lon,lat=lat,
##               unit="Pa",longname="mean sea level pressure",
##               param="slp",quality=NULL,src="ERAinterim")
## #slp <- subset(slp,it="december")
## #Z <- slp
## cyclones <- CCI(slp,it="december",fname="cyclones.ERAint.1979.12.rda")




