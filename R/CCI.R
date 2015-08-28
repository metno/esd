# K Parding, 29.05.2015

CCI <- function(Z,m=14,it=NULL,is=NULL,cyclones=TRUE,
                accuracy=NULL,label=NULL,fname="cyclones.rda",
                plot=FALSE,verbose=FALSE) {
  stopifnot(inherits(Z,'field'))
  Z <- subset(Z,it=it,is=is)
  if (any(longitude(Z)>180)) Z <- g2dl(Z,greenwich=FALSE)
    
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

  ## Mask the cylcones already selected
  #P.lowy[lows] <- 0; P.lowx[lows] <- 0
  ## widen mask in x-direction
  #P.lowx[,2:(nx-1),] <- P.lowx[,1:(nx-2),] | P.lowx[,2:(nx-1),]
  #P.lowx[,1:(nx-2),] <- P.lowx[,1:(nx-2),] | P.lowx[,2:(nx-1),]
  ## widen mask in y-direction
  #P.lowy[,,2:(ny-1)] <- P.lowy[,,1:(ny-2)] | P.lowy[,,2:(ny-1)]
  #P.lowy[,,1:(ny-2)] <- P.lowy[,,1:(ny-2)] | P.lowy[,,2:(ny-1)]

  ## Find zero crossings in both directions
  ## Plowx & P.lowy are matrices with 0's and 1's.
  lows <- (P.lowy & P.lowx)
  pcent <- 0.5*(px[lows]+py[lows])
  strength <- order(pcent)
  if (!cyclones) strength <- reverse(strength)

  ## Handle time index
  t <- index(Z)
  if (inherits(t,"POSIXt")) t <- as.numeric(format(t,"%Y%m%d%H%M"))
    
  ## Remove secondary cyclones near a deeper one (same cyclonic system):
  if(verbose) print("Remove secondary cyclones")
  mindistance <- 5E4 # minimum distance between cyclones [m]
  lon<-rep(lonXY,nt); dim(lon)<-c(nx-1,ny-1,nt); lon<-aperm(lon,c(3,1,2)) 
  lat<-rep(latXY,nt); dim(lat)<-c(nx-1,ny-1,nt); lat<-aperm(lat,c(3,1,2))
  date<-rep(t,(nx-1)*(ny-1)); dim(date)<-c(nt,nx-1,ny-1)
  lon <- lon[lows]; date <- date[lows]; lat<-lat[lows]
  del <- rep(TRUE,length(date))
  for (d in unique(date)) {
    if (inherits(t,'Date')) d <- as.Date(d)
    i <- which(date==d)
    distance <- apply(cbind(lon[i],lat[i]),1,
     function(x) suppressWarnings(distAB(x[1],x[2],lon[i],lat[i])))
    diag(distance) <- NA; distance[lower.tri(distance)] <- NA
    del.i <- which(distance<1E6,arr.ind=TRUE)
    if(any(del.i)) {
      col.del <- rep(1,dim(del.i)[1])
      col.del[ strength[i][del.i[,1]] <= strength[i][del.i[,2]] ] <- 2
      del.i <- del.i[cbind(seq(1,dim(del.i)[1]),col.del)]
      del[i[del.i]] <- FALSE
    }
  }
  lows[lows] <- del
  lon <- lon[del]
  lat <- lat[del]
  date <- date[del]
  pcent <- pcent[del]
  
  # put cyclones in order of date
  i <- order(date)
  lon <- lon[i]
  lat <- lat[i]
  date <- date[i]
  pcent <- pcent[i]
  strength <-  order(pcent)
  
  # Pressure gradient
  if(verbose) print("Pressure gradient")
  rho <- 1.2922
  dpsl <- sqrt(DX^2+DY^2)
  if (attr(Z,"unit")=="hPa") dpsl <- dpsl*100

  # Find points of inflexion (2nd derivative==0) to estimate the storm radius
  # and maximum speed and pressure gradient
  if(verbose) print("Find points of inflexion")
  rmin <- 10; rmax <- 2000 # km
  NX <- dim(lonXY)[1]; NY <- dim(latXY)[2]
  lonXX <- rep(0.5*(lonXY[2:NX,2:NY]+lonXY[1:(NX-1),2:NY]),nt)
  dim(lonXX) <- c(NX-1,NY-1,nt); lonXX <- aperm(lonXX,c(3,1,2))
  latXX <- rep(0.5*(latXY[2:NX,2:NY]+latXY[2:NX,1:(NY-1)]),nt)
  dim(latXX) <- c(NX-1,NY-1,nt); latXX <- aperm(latXX,c(3,1,2))
  dateXX <- rep(t,(NX-1)*(NY-1)); dim(dateXX) <- c(nt,NX-1,NY-1)
  radius <- rep(NA,length(date))
  max.dslp <- rep(NA,length(date))
  max.speed <- rep(NA,length(date))
  max.vg <- rep(NA,length(date))
  if (plot) data(geoborders)
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
    ## geostrophic wind
    fi <- 2*7.29212*1E-5*sin(pi*latXY[1,ilat]/180)
    vg <- dpi/(fi*rho)
    ## gradient wind
    v.grad <- -0.5*fi*pi*ri*(1 - sqrt(1 + 4*vg/(fi*ri)))
    radius[i] <- mean(ri)#ri[which.max(dpi)]
    max.dslp[i] <- mean(dpi)#dpi[which.max(dpi)]
    max.speed[i] <- mean(v.grad)#v.grad[which.max(dpi)]
    max.vg[i] <- mean(vg)#vg[which.max(dpi)]

    if (plot & i<4) {
      pxi <- px[date[i]==t,,];  pyi <- py[date[i]==t,,]
      lon.i <- lonXY[ilon,1]; lat.i <- latXY[1,ilat] 
      xi <- lonXY[,1]; yi <- latXY[1,]; zi <- pxi
      if (!(all(diff(xi)>0))) {xi <- rev(xi); zi <- apply(zi,2,rev)}
      if (!(all(diff(yi)>0))) {yi <- rev(yi); zi <- t(apply(t(zi),2,rev))}
      dev.new()
      image(xi,yi,zi,main=t[i],col=colscal(col="t2m",n=12,rev=FALSE),
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
      points(lon.a,lat.a,col="black",pch="-",cex=1,lwd=0.5)
      points(lon.b,lat.b,col="black",pch="|",cex=0.5,lwd=0.5)
      points(lon[date==date[i]],lat[date==date[i]],pch=21,lwd=2,
             bg="white",col="black",cex=2)
      points(lon[i],lat[i],pch=4,col="black",cex=1)
      dev.new()
      plot(lonXY[,1],pxi[,latXY[1,]==lat[i]],lty=1,type="l",main=t[i],
           xlab="lon",ylab="slp (hPa)")
      points(lon[i],pxi[lonXY[,1]==lon[i],latXY[1,]==lat[i]],col="blue",pch=19)
      points(lonXY[inflx<0,1],pxi[inflx<0,latXY[1,]==lat[i]],col="red",pch=1)
      points(lon.i[lat.i==lat[i]],pxi[ilon[lat.i==lat[i]],latXY[1,]==lat[i]],
             col="red",pch=20)
      dev.new()
      plot(latXY[1,],pyi[lonXY[,1]==lon[i],],lty=1,type="l",main=t[i],
           xlab="lat",ylab="slp (hPa)")
      points(lat[i],pyi[lonXY[,1]==lon[i],latXY[1,]==lat[i]],col="blue",pch=19)
      points(latXY[1,infly<0],pyi[lonXY[,1]==lon[i],infly<0],col="red",pch=1)
      points(lat.i[lon.i==lon[i]],pyi[lonXY[,1]==lon[i],ilat[lon.i==lon[i]]],
             col="red",pch=20)
    }
  }

  ## Remove temporary variable and release the memory:
  rm('lonXX','latXX','dateXX','inflx','infly'); gc(reset=TRUE)
  
  # Check units
  if (attr(Z,"unit")=='Pa') pcent <- pcent*100 # Pa -> hPa
  max.dslp <- max.dslp*100 # Pa -> hPa

  # Arrange results
  if (inherits(index(Z),"POSIXt")) date <- strptime(date,"%Y%m%d%H%M")
  dd <- as.numeric(strftime(date,"%Y%m%d"))
  hh <- as.numeric(strftime(date,"%H"))
  results <- data.frame(date=dd,time=hh,lon=lon,lat=lat,
                  pcent=pcent,max.dslp=max.dslp,
                  max.speed=max.speed,radius=radius)
  attr(results,"label") <- label
  attr(results,"dx") <- dx
  attr(results,"dy") <- dy
  attr(results,"units") <- c("date","hour CET","degrees","degrees",
                             "hPa","hPa/m","m/s","km")  
  if (cyclones) {
    attr(results,"longname") <- "low pressure systems identified with CCI method"
    attr(results,"variable") <- "cyclones"
  } else {
    attr(results,"longname") <- "high-pressure systems identified with CCI method"
    attr(results,"variable") <- "anti-cyclones"
  }
  attr(results,"source") <- attr(Z,"source")
  attr(results,"file") <- attr(Z,"file")
  attr(results,"version") <- "CCI in esd v1.0 (after August 25, 2015)"
  class(results) <- c("events",class(results))
  if(!is.null(fname)) save(file=fname,results)
  invisible(results)
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
