# K Parding, 29.05.2015

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
## slp <- subset(slp,it="december")
## Z <- slp
##cyclones <- CCI(slp,it="december",fname="cyclones.ERAint.1979.12.rda")

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

  ## Find zero crossings in both directions
  lows <- (P.lowy & P.lowx)
  pcent <- 0.5*(px[lows]+py[lows])
  strength <- order(pcent)
  if (!cyclones) strength <- reverse(strength)

  ## Handle time index
  t <- index(Z)
  if (inherits(t,"POSIXt")) t <- as.numeric(format(t,"%Y%m%d%H%M"))
    
  ## Remove secondary cyclones near a deeper one (same cyclonic system):
  if(verbose) print("Remove secondary cyclones")
  mindistance <- 1E4 # minimum distance between cyclones [m]
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
  strength <-  order(pcent)
  
  # Geostrophic wind speed
  # CHECK THE UNITS
  if(verbose) print("Geostrophic wind")
  f <- 1.47e-04*sin(pi*latXY/180)
  f[abs(latXY)<10] <- NA # not valid close to equator
  rho <- 1.2922
  A <- rep(f*rho,nt); dim(A) <- c(nx-1,ny-1,nt); A <- aperm(A,c(3,1,2))
  dpsl <- sqrt(DX^2+DY^2)
  if (attr(Z,"unit")=="hPa") dpsl <- dpsl*100
  wind <- dpsl/A
 
  # Find points of inflexion (2nd derivative==0) to estimate the storm radius
  # and maximum 
  if(verbose) print("Find points of inflexion")
  rmin <- 10; rmax <- 2000 # km
  NX <- dim(lonXY)[1]; NY <- dim(latXY)[2]
  lonXX <- rep(0.5*(lonXY[2:NX,2:NY]+lonXY[1:(NX-1),2:NY]),nt)
  dim(lonXX) <- c(NX-1,NY-1,nt); lonXX <- aperm(lonXX,c(3,1,2))
  latXX <- rep(0.5*(latXY[2:NX,2:NY]+latXY[2:NX,1:(NY-1)]),nt)
  dim(latXX) <- c(NX-1,NY-1,nt); latXX <- aperm(latXX,c(3,1,2))
  dateXX <- rep(t,(NX-1)*(NY-1)); dim(dateXX) <- c(nt,NX-1,NY-1)
  inflx <- DX2[,2:NX,2:NY]*DX2[,1:(NX-1),2:NY]<0 &
           DX2[,2:NX,1:(NY-1)]*DX2[,1:(NX-1),1:(NY-1)]<0 &
           !is.na(DX2[,2:NX,2:NY]*DX2[,1:(NX-1),1:(NY-1)])
  infly <- DY2[,2:NX,2:NY]*DY2[,2:NX,1:(NY-1)]<0 &
           DY2[,1:(NX-1),2:NY]*DY2[,1:(NX-1),1:(NY-1)]<0 &
           !is.na(DY2[,2:NX,2:NY]*DY2[,1:(NX-1),1:(NY-1)])
  radius <- rep(NA,length(date))
  lon.radius <- rep(NA,length(date))
  lat.radius <- rep(NA,length(date))
  max.dpsl <- rep(NA,length(date))
  max.speed <- rep(NA,length(date))
  for (i in order(date)) {
    ilon <- abs(lonXX[1,,1]-lon[i])<dx
    ilat <- abs(latXX[1,1,]-lat[i])<dy
    lon.i <- c(lonXX[1,,ilat][inflx[t==date[i],,ilat]],
               lonXX[1,ilon,][infly[t==date[i],ilon,]])
    lat.i <- c(latXX[1,,ilat][inflx[t==date[i],,ilat]],
               latXX[1,ilon,][infly[t==date[i],ilon,]])
    distance <- distAB(lon[i],lat[i],lon.i,lat.i)*1E-3
    distance[distance<rmin | distance>rmax] <- NA
    i.min <- which.min(distance)
    if(any(i.min)) {
      radius[i] <- distance[i.min]
      lon.radius[i] <- lon.i[i.min]
      lat.radius[i] <- lat.i[i.min]
      i.x <- which(abs(lonXY[,1]-lon[i])<=abs(lon.radius[i]-lon[i]))
      i.y <- which(abs(latXY[1,]-lat[i])<=abs(lat.radius[i]-lat[i]))
      max.dpsl[i] <- max(dpsl[t==date[i],i.x,i.y],na.rm=T)
      max.speed[i] <- max(wind[t==date[i],i.x,i.y],na.rm=T)
    }
  }
  ## Remove temporary variable and release the memory:
  rm('lonXX','latXX','dateXX','inflx','infly'); gc(reset=TRUE)
  
  if(plot) {
    data(geoborders)
    for (i in seq(1,8)) {
      xi <- lonXY[,1]; yi <- latXY[1,]; zi <- px[i,,]
      if (!(all(diff(xi)>0))) {xi <- rev(xi); zi <- apply(zi,2,rev)}
      if (!(all(diff(yi)>0))) {yi <- rev(yi); zi <- t(apply(t(zi),2,rev))}
      dev.new(); image(xi,yi,zi,main=t[i])
      lines(geoborders)
      points(lon[date==t[i]],lat[date==t[i]],pch=1)
      matlines(rbind(lon[date==t[i]],lon.radius[date==t[i]]),
               rbind(lat[date==t[i]],lat.radius[date==t[i]]),
               lty=1,lwd=1,col='black')
    }
  } 

  # Check units
  if (attr(Z,"unit")=='Pa') pcent <- pcent*1E2
  if (inherits(index(Z),"POSIXt")) date <- strptime(date,"%Y%m%d%H%M")
  
  # Add attributes
  attr(lon,'units') <- 'degrees'
  attr(lat,'units') <- 'degrees'
  attr(pcent,'units') <- 'hPa'
  attr(max.dpsl,'units') <- 'Pa/m'
  attr(max.dpsl,'location') <-
    'at inflexion points at lon/lat lines through storm center'
  attr(max.speed,'units') <- 'm/s'
  attr(max.speed,'location') <-
    'at inflexion points at lon/lat lines through storm center'
  attr(radius,'units') <- 'km'
  
  results <- list(lon=lon,lat=lat,tim=date,pcent=pcent,
                  yy=year(date),mm=month(date),dd=day(date),
                  i=1:length(date),label=label,
                  max.dpsl=max.dpsl,max.speed=max.speed,
                  radius=radius,dx=dx,dy=dy,
                  version="cyclones v2.1-1 (after June 3, 2015)")
  save(file=fname,results) 
  invisible(results)
}
