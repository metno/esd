#' @export
regfit <- function(z,cal.dat,terms) {
  cal.dat$Z <- z
  model <- eval(parse(text=paste('lm(Z ~ ',terms,',data=cal.dat)')))
  ln <- length(unlist(strsplit(terms,split='\\+')))
  modelcoefs <- matrix(rep(0,(ln+1)*4),ln+1,4)
  cf <- summary(model)$coefficients
  colnames(modelcoefs) <- colnames(cf)
  modelcoefs[seq(nrow(cf)),] <- cf
  #modelcoefs <- summary(model)$coefficients
  return(modelcoefs)
}

#' Derivatives 
#' 
#' \code{dX}, \code{dY}, and \code{dT} are functions to estimate derivatives for
#' gridded field objects based on a fit to truncated Fourier series.
#' The three functions give the x-, y- and time derivatives respectively.
#' See Benestad & Chen (2006) 'The use of a
#' Calculus-based Cyclone Identification method for generating storm
#' statistics' (Tellus A 58A, 473-486, doi:10.1111/j.1600-0870.2006.00191) for
#' more details.
#'
#' \code{regfit} is a help function for generating a model for fitting the profile.
#'
#' @aliases dX dY dT regfit
#' @param Z A field object 
#' @param m number of harmonics for fitting the Fourier series 
#' @param mask.bad mask missing data 
#' @param plot if TRUE show plot 
#' @param r radius of the Earth (m) 
#' @param accuracy resolution of output
#' @param progress show the progress  
#' @param verbose show diagnostics of the progress 
#' 
#' @return a list with several comonents:
#' 
#' \item{Z}{original data} \item{a}{Fourier coefficients for cosine}
#' \item{b}{Fourier coeffieicnes for sine} \item{z0}{defunct?} \item{dZ}{The
#' component contains the first derivative.} \item{dZ2}{The component contains
#' the second derivative (quicker to do both in one go).} \item{lon}{longitude}
#' \item{lat}{latitude} \item{dx}{spatial resolution} \item{span}{spatial
#' extent}
#' 
#' @examples
#' data(slp.ERA5)
#' slp.dx <- dX(slp.ERA5,verbose=TRUE)
#' map(slp.dx$Z) # map of SLP 
#' map(slp.dx$dZ) # map of first derivative in longitude direction
#' map(slp.dx$dZ2) # map of second derivative in longitude direction
#' \dontrun{
#' u10 <- retrieve('~/Downloads/Jan2018_ERAINT_uvp.nc',param='u10')
#' v10 <- retrieve('~/Downloads/Jan2018_ERAINT_uvp.nc',param='v10')
#' ## Estimate the vorticity
#' zeta <- dX(v10)$dZ - dY(u10)$dZ
#' zeta <- attrcp(u10,zeta)
#' class(zeta) <- class(u10)
#' attr(zeta,'variable') <- 'vorticity'
#' attr(zeta,'unit') <- '1/s'
#' map(subset(zeta,it=1),projection='np')
#' }
#' 
#' @export dX
dX <- function(Z,m=10,mask.bad=TRUE,plot=FALSE,r=6.378e06,
               accuracy=NULL,progress=TRUE,verbose=FALSE) {

  ## Convert the field object into 3D objects with lon-lat dimensions
  ## seperated.

  if (verbose) print('dX')
  z <- as.pattern(Z)
  lon <- lon(Z)
  lat <- lat(Z)
  
  if (is.null(accuracy)) accuracy <-  max(diff(lon))
  LON <- seq(min(lon),max(lon),by=accuracy)
  NX <- length(LON)

  ny <- length(lat)
  nx <- length(lon)
  if (verbose) str(Z)
  if (is.matrix(z)) {
    dim(z) <- c(dim(z),1)
    nt <- 1
    t <- 1 
  } else if (!is.null(index(Z))) {
    nt <- length(index(Z))
    t <- index(Z)
  } else {
    nt <- 1
    t <- 1
  }
  
  if (is.null(m)) m <- nx
  m <- min(nx,m)
  if(verbose) print(paste("Number of harmonics:",m))
  theta <- pi*lon/180
  phi <- pi*lat/180
  mask <- !is.finite(z)

  ## Generate frequencies for the siusoids: rows of matrix Wi[1..m,n..nt]
  ## contain the harmonics 1..m using the outer product '%o%'
  Wi <- 2*pi/nx*(1:m)
  Wii <- Wi %o% seq(1,nx,by=1)
  rownames(Wii) <- paste('harmonic',1:m)

  ## Generate the terms to include in the regression:
  terms <- paste('c',(1:m), ' + s',(1:m),sep='',collapse=' + ')

  ## generate data.frame containing the original data and the harmonics
  cal.dat <- eval(parse(text=paste('data.frame(',
                          paste('c',(1:m),'= cos(Wii[',(1:m),',]),',
                                's',(1:m),'= sin(Wii[',(1:m),',])',
                                sep='',collapse=', '),')')))

  ## Set up matrices containing the coefficients:
  a <- rep(0,ny*m*nt)
  dim(a) <- c(m,ny,nt)
  b <- a; 
  z0 <- matrix(rep(NA,nt*ny),ny,nt)

  ## Loop over time steps and apply the harmonic fit to each latitude:
  t1 <- Sys.time()
  if(progress) pb <- txtProgressBar(style=3)
  if (verbose) print(paste('nt=',nt))
  for ( it in 1:nt ) {
    if(progress) setTxtProgressBar(pb,it/nt) 
    ## Create a matrix containing m harmonic fits for ny latitudes:
    #beta <- apply(z[,,it],2,regfit,cal.dat=cal.dat,terms=terms) ## KMP 2016-02-01
    beta <- apply(z[,,it],2,function(x) {
        y <- regfit(x,cal.dat=cal.dat,terms=terms)
        invisible(y[,"Estimate"]) })
    ## The constant
    z0[,it] <- beta[1,]
    a[1:floor(dim(beta)[1]/2),,it] <- beta[seq(2,dim(beta)[1],by=2),] ## KMP 2016-02-01
    b[1:floor(dim(beta)[1]/2),,it] <- beta[seq(3,dim(beta)[1],by=2),] ## KMP 2016-02-01
    #a[,,it] <- beta[seq(2,2*m+1,by=2),] ## KMP 2016-02-01
    #b[,,it] <- beta[seq(3,2*m+1,by=2),] ## KMP 2016-02-01
  }
  t2 <- Sys.time()
  if (verbose) print(paste('Taking dX of the field took',
                round(as.numeric(t2-t1,units="secs")),'s'))

  ## Reorganise the data and take the inner product:
  a2d <- a;  dim(a2d) <- c(m,nt*ny)
  b2d <- b;  dim(b2d) <- c(m,nt*ny)
  zz <- t(cos(Wii))%*%a2d  + t(sin(Wii))%*% b2d
  dim(zz) <- c(nx,ny,nt)
  zz0 <- rep(z0,nx); dim(zz0) <- c(ny,nt,nx)
  zz0 <- aperm(zz0,c(3,1,2))
  
  ## Find the regression fit:
  if (verbose) print('Find the best-fit')
  z.fit <- zz0 + zz
  dim(z.fit) <- c(nx*ny,nt)
  Z.fit <- zoo(t(z.fit),order.by=t)
  Z.fit <- as.field(Z.fit,lon=lon(Z),lat=lat(Z),param=varid(Z),unit=unit(Z),
                    longname=paste('fitted',attr(Z,'longname')),
                    greenwich = attr(Z,'greenwich'),
                    calendar = attr(Z,'calendar'),aspect='fitted')
  ## Remove temporary variable and release the memory:
  rm('zz0','zz','z.fit'); gc(reset=TRUE)

  dx <- rep(mapply(distAB,lon[1],lat,lon[2],lat),m*nt)
  dim(dx) <- c(ny,m,nt); dx <- aperm(dx,c(2,1,3))
  dx[,cos(phi)<1E-3,] <- NA
  
  ## Derive the first derivative:
  if (verbose) print('Find the first derivative')
  a2d <- a/dx;  dim(a2d) <- c(m,ny*nt)
  b2d <- b/dx;  dim(b2d) <- c(m,ny*nt)
  dz <- t(-Wi*sin(Wii))%*%a2d + t(Wi*cos(Wii))%*%b2d
  dim(dz) <- c(nx*ny,nt)
  dZ.fit <- zoo(t(dz),order.by=t)
  dZ.fit <- as.field(dZ.fit,lon=lon(Z),lat=lat(Z),
                    param=paste('d*',varid(Z)),
                    unit=paste0(unit(Z),'/dx'),
                    longname=paste('fitted',attr(Z,'longname')),
                    greenwich = attr(Z,'greenwich'),
                    calendar = attr(Z,'calendar'),aspect='fitted')  
  
  ## Derive the second derivative:
  if (verbose) print('Find the second derivative')
  a2d2 <- a/dx^2;  dim(a2d2) <- c(m,nt*ny)
  b2d2 <- b/dx^2;  dim(b2d2) <- c(m,nt*ny)
  dz2 <- t(-Wi^2*cos(Wii))%*%a2d2 + t(-Wi^2*sin(Wii))%*%b2d2
  dim(dz2) <- c(nx*ny,nt)
  dZ2.fit <- zoo(t(dz2),order.by=t)
  dZ2.fit <- as.field(dZ2.fit,lon=lon(Z),lat=lat(Z),param=varid(Z),
                    unit=paste0(unit(Z),'2/dx2'),
                    longname=paste('fitted',attr(Z,'longname')),
                    greenwich = attr(Z,'greenwich'),
                    calendar = attr(Z,'calendar'),aspect='fitted')  

  if (mask.bad) {
    mask.fit <- aperm(mask,c(3,1,2))
    dim(mask.fit) <- dim(Z.fit)
    Z.fit[mask.fit] <- NA
    dZ.fit[mask.fit] <- NA
    dZ2.fit[mask.fit] <- NA 
  }
  results <- list(Z=Z,a=a,b=b,z0=z0,dZ=dZ.fit,dZ2=dZ2.fit,Z.fit=Z.fit,
                  lon=lon,lat=lat,dx=diff(lon)[1],span=range(lon))
  class(results) <- "map"
  attr(results,"long_name") <- "x-derivative"
  attr(results,"spatial units") <- "hPa/m"
  attr(results,"descr") <- "dX()"
  invisible(results)
}

