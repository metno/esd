# K Parding, 29.05.2015

#library(esd)
#slp <- slp.ERAINT()
#Z <- subset(slp,is=list(lon=c(-180,180),lat=c(0,90)),it='djf')

CCI <- function(Z,m=14,nsim=10,it=NULL,is=NULL,
                cyclones=TRUE,accuracy=NULL,verbose=FALSE) {

  stopifnot(inherits(Z,'field'))
  Z <- subset(Z,it=it,is=is)
  
  ## Calculate first and second derivative
  if(verbose) print("Calculate first and second derivative")
  resx <- dX(Z,m=m,accuracy=accuracy,verbose=verbose)
  resy <- dY(Z,m=m,accuracy=accuracy,verbose=verbose)

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
  dy <- as.matrix(resy$dZ); dim(dy) <- c(nt,nx,ny)
  dx <- as.matrix(resx$dZ); dim(dx) <- c(nt,nx,ny)
  dx2 <- as.matrix(resx$dZ2); dim(dx2) <- c(nt,nx,ny)
  dy2 <- as.matrix(resy$dZ2); dim(dy2) <- c(nt,nx,ny)
  px <- 0.25*(Zx[,1:(nx-1),2:ny] + Zx[,2:nx,2:ny] +
              Zx[,1:(nx-1),1:(ny-1)] + Zx[,2:nx,1:(ny-1)])
  py <- 0.25*(Zy[,1:(nx-1),2:ny] + Zy[,2:nx,2:ny] +
              Zy[,1:(nx-1),1:(ny-1)] + Zy[,2:nx,1:(ny-1)])

  ## Search for zero crossing in y-direction
  if(verbose) print("Find zero crossing of first derivative in y-direction")
  P.lowy <- rep(0,nt*(nx-1)*(ny-1)); dim(P.lowy) <- c(nt,nx-1,ny-1) 
  dy11 <- 0.5*(dy[,2:nx,2:ny]+dy[,1:(nx-1),2:ny])
  dy12 <- 0.5*(dy[,2:nx,1:(ny-1)]+dy[,1:(nx-1),1:(ny-1)])
  dy21 <- 0.5*(dy2[,2:nx,2:ny]+dy2[,1:(nx-1),2:ny])
  dy22 <- 0.5*(dy2[,2:nx,1:(ny-1)]+dy2[,1:(nx-1),1:(ny-1)])
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
  dx11 <- 0.5*(dx[,2:nx,2:ny]+dx[,2:nx,1:(ny-1)])
  dx12 <- 0.5*(dx[,1:(nx-1),2:ny]+dx[,1:(nx-1),1:(ny-1)])
  dx21 <- 0.5*(dx2[,2:nx,2:ny]+dx2[,2:nx,1:(ny-1)])
  dx22 <- 0.5*(dx2[,1:(nx-1),2:ny]+dx2[,1:(nx-1),1:(ny-1)])
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
  
  ## Remove secondary cyclones near a deeper one (same cyclonic system):
  if(verbose) print("Remove secondary cyclones")
  mindistance <- 1E4 # minimum distance between cyclones [m]
  lon<-rep(lonXY,nt); dim(lon)<-c(nx-1,ny-1,nt); lon<-aperm(lon,c(3,1,2)) 
  lat<-rep(latXY,nt); dim(lat)<-c(nx-1,ny-1,nt); lat<-aperm(lat,c(3,1,2)) 
  date<-rep(index(Z),(nx-1)*(ny-1)); dim(date)<-c(nt,nx-1,ny-1)
  lon <- lon[lows]; date <- date[lows]; lat<-lat[lows]
  del <- rep(TRUE,length(date))
  for (d in unique(date)) {
    i <- which(date==as.Date(d))
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
  
  # Order lon, lat and pressure of lows
  if(verbose) print("Order cyclones by date")
  lon <- lon[strength]
  lat <- lat[strength]
  date <- date[strength]
  pcent <- pcent[strength]
  strength <- strength[strength]

  # Geostrophic wind speed
  # CHECK THE UNITS
   if(verbose) print("Geostrophic wind")
  f <- 1.47e-04*sin(pi*latXY/180)
  f[abs(latXY)<10] <- NA # not valid close to equator
  rho <- 1.2922
  A <- rep(f*rho,nt); dim(A) <- c(nx-1,ny-1,nt); A <- aperm(A,c(3,1,2))
  u <- -DX/A*100
  v <- DY/A*100
  wind <- sqrt(u^2 + v^2)
  # SLP in hPa -> Pa, but dx and dy in m.
 
  # Find points of inflexion (2nd derivative==0) to estimate the storm radius
  # CAN I DO THIS WITHOUT LOOPING OVER ALL STORMS?
  if(verbose) print("Find points of inflexion")
  NX <- dim(lonXY)[1]; NY <- dim(latXY)[2]
  inflx <- DX2[,2:NX,]*DX2[,1:(NX-1),]<0
  infly <- DY2[,,2:NY]*DY2[,,1:(NY-1)]<0
  date.old <- 0
  radius <- c()
  for (i in order(date)) {
    i.inflx <- inflx[index(Z)==date[i],,latXY[1,]==lat[i]]
    i.infly <- infly[index(Z)==date[i],lonXY[,1]==lon[i],]
    i.x <- which(i.inflx)[which.min(abs(lonXY[i.inflx,1]-lon[i]))]
    i.y <- which(i.infly)[which.min(abs(latXY[1,i.infly]-lat[i]))]
    r.x <- distAB(lon[i],lat[i],lonXY[i.x,1],lat[i])
    r.y <- distAB(lon[i],lat[i],lon[i],latXY[1,i.y])
    r.i <- sqrt(r.x^2 + r.y^2)
    radius <- c(radius,r.i)
    if(plot) {
      if(i==1 | date[i]!=date.old) {
        dev.new()
        image(resx$lon,resy$lat,Zy[index(Z)==date[i],,],main=date[i])
      }
      points(lon[i],lat[i])
      lines(c(lon[i],lonXY[i.x,1]),c(lat[i],lat[i]))
      lines(c(lon[i],lon[i]),c(lat[i],latXY[1,i.y]))
    }
    date.old <- date[i]
  }
 

}
