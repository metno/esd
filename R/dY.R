# Calculates the x-derivatives for gridded data in longitude-
# latitude coordinates. After Gill (1982) p. 94
#
# dA/dx = 1/(r) * d/d THETA  
#
# where PHI is the latitude in radians and THETA the longitude.
# R.E. Benestad & Kajsa Parding, 2015-05-26

#library(esd)
#slp <- slp.ERAINT()
#slp.dx <- dX(slp,verbose=TRUE)

dY <- function(Z,m=10,mask.bad=TRUE,plot=FALSE,r=6.371e06,
               chk.conf=1,accuracy=NULL,verbose=FALSE) {
# REB: 08.02.2007: To improve accuracy/spatial resolution of the fit

  ## Convert the field object into 3D objects with lon-lat dimensions
  ## seperated.

  if (verbose) print('dY')
  z <- as.pattern(Z)
  lon <- lon(Z)
  lat <- lat(Z)
  
  if (is.null(accuracy)) accuracy <-  max(diff(lon))
  LON <- seq(min(lon),max(lon),by=accuracy)
  NX <- length(LON)

  ny <- length(lat)
  nx <- length(lon)
  nt <- length(index(Z))
  
  if (is.null(m)) m <- nx
  m <- min(nx,m)
  theta <- pi*lon/180
  phi <- pi*lat/180
  mask <- !is.finite(z)

  ## Generate frequencies for the siusoids: rows of matrix Wi[1..m,n..nt]
  ## contain the harmonics 1..m using the outer product '%o%'
  Wi <- 2*pi/ny*(1:m)
  Wii <- Wi %o% seq(1,ny,by=1) 
  rownames(Wii) <- paste('harmonic',1:m)

  ## Generate the terms to include in the regression:
  terms <- paste('c',(1:m), ' + s',(1:m),sep='',collapse=' + ')

  ## generate data.frame containing the original data and the harmonics
  cal.dat <- eval(parse(text=paste('data.frame(',
                          paste('c',(1:m),'= cos(Wii[',(1:m),',]),',
                                's',(1:m),'= sin(Wii[',(1:m),',])',
                                sep='',collapse=', '),')')))

  ## Set up matrices containing the coefficients:
  a <- rep(0,nx*m*nt)
  dim(a) <- c(m,nx,nt)
  b <- a; 
  z0 <- matrix(rep(NA,nt*nx),nx,nt)

  ## Loop over time steps and apply the harmonic fit to each latitude:
  t1 <- Sys.time()
  pb <- txtProgressBar(style=3)
  for ( it in 1:nt ) {
    setTxtProgressBar(pb,it/nt) 
    ## Create a matrix containing m harmonic fits for ny latitudes:
    beta <- apply(z[,,it],1,regfit,cal.dat=cal.dat,terms=terms)
    ## The constant
    z0[,it] <- beta[1,]
    a[,,it] <- beta[seq(2,2*m+1,by=2),]
    b[,,it] <- beta[seq(3,2*m+1,by=2),]
  }
  t2 <- Sys.time()
  if (verbose) print(paste('Taking dY of the field took',
                round(as.numeric(t2-t1,units="secs")),'s'))

  ## Reorganise the data and take the inner product:
  a2d <- a;  dim(a2d) <- c(m,nt*nx)
  b2d <- b;  dim(b2d) <- c(m,nt*nx)
  zz <- t(cos(Wii))%*%a2d  + t(sin(Wii))%*%b2d
  dim(zz) <- c(ny,nx,nt); zz <- aperm(zz,c(2,1,3))
  zz0 <- rep(z0,ny); dim(zz0) <- c(nx,nt,ny); zz0 <- aperm(zz0,c(1,3,2))
  
  ## Find the regression fit:
  if (verbose) print('Find the best-fit')
  z.fit <- zz0 + zz
  dim(z.fit) <- c(nx*ny,nt)
  Z.fit <- zoo(t(z.fit),order.by=index(Z))
  Z.fit <- as.field(Z.fit,lon=lon(Z),lat=lat(Z),param=varid(Z),unit=unit(Z),
                    longname=paste('fitted',attr(Z,'longname')),
                    greenwich = attr(Z,'greenwich'),aspect='fitted')
  ## Remove temporary variable and release the memory:
  rm('zz0','zz','z.fit'); gc(reset=TRUE)

  dy <- diff(lat)[1]*r

  ## Derive the first derivative:
  if (verbose) print('Find the first derivative')
  a2d <- a/dy;  dim(a2d) <- c(m,nx*nt)
  b2d <- b/dy;  dim(b2d) <- c(m,nx*nt)
  dz <- t(-Wi*sin(Wii))%*%a2d +  t(Wi*cos(Wii))%*%b2d
  dim(dz) <- c(ny,nx,nt); dz <- aperm(dz,c(2,1,3)); dim(dz) <- c(nx*ny,nt)
  dZ.fit <- zoo(t(dz),order.by=index(Z))
  dZ.fit <- as.field(dZ.fit,lon=lon(Z),lat=lat(Z),
                    param=paste('d*',varid(Z)),
                    unit=paste(unit(Z),'/dx'),
                    longname=paste('fitted',attr(Z,'longname')),
                    greenwich = attr(Z,'greenwich'),aspect='fitted')  
  
  ## Derive the second derivative:
  if (verbose) print('Find the second derivative')
  a2d2 <- a/dy^2;  dim(a2d2) <- c(m,nx*nt)
  b2d2 <- b/dy^2;  dim(b2d2) <- c(m,nx*nt)
  dz2 <- t(-Wi^2*cos(Wii))%*%a2d2 + t(-Wi^2*sin(Wii))%*%b2d2
  dim(dz2) <- c(ny,nx,nt); dz2 <- aperm(dz2,c(2,1,3)); dim(dz2) <- c(nx*ny,nt)
  dZ2.fit <- zoo(t(dz2),order.by=index(Z))
  dZ2.fit <- as.field(dZ2.fit,lon=lon(Z),lat=lat(Z),param=varid(Z),
                    unit=paste(unit(Z),'^2/dy^2'),
                    longname=paste('fitted',attr(Z,'longname')),
                    greenwich = attr(Z,'greenwich'),aspect='fitted')  

  if (mask.bad) {
    mask.fit <- aperm(mask,c(3,1,2))
    dim(mask.fit) <- dim(Z.fit)
    Z.fit[mask.fit] <- NA
    dZ.fit[mask.fit] <- NA
    dZ2.fit[mask.fit] <- NA 
  }
  results <- list(Z=Z,a=a,b=b,z0=z0,dZ=dZ.fit,dZ2=dZ2.fit,Z.fit=Z.fit,
                  lon=lon,lat=lat,dy=diff(lat)[1],span=range(lat))
  class(results) <- "map"
  attr(results,"long_name") <- "y-derivative"
  attr(results,"spatial units") <- "km"
  attr(results,"descr") <- "dY()"
  invisible(results)
}

