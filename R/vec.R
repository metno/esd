#' Draw vectors on maps for esd objects
#' 
#' Plot arrows for e.g. a wind field or ocean currents. 
#' @aliases vec
#' @seealso map
#' 
#' @param x the zonal comonents of the vectors
#' @param y the meridional comonents of the vectors
#' @param new TRUE: create a new graphic device.
#' @param projection Projections: c("lonlat","sphere","np","sp") - the latter
#' gives stereographic views from the North and south poles.
#' @param xlim see \code{\link{plot}} - only used for 'lonlat' and 'sphere'
#' projections.
#' @param ylim see \code{\link{plot}} - only used for 'lonlat' and 'sphere'
#' projections.
#' @param n The number of colour breaks in the color bar
#' @param breaks graphics setting - see \code{\link{image}}
#' @param type graphics setting - colour shading or contour
#' @param gridlines Only for the lon-lat projection
#' @param lonR Only for the spherical projection 
#' @param latR Only for the spherical projection
#' @param axiR Only for the spherical projection 
#' @param a used  to scale the length of the arrows
#' @param r used to make a 3D effect of plotting the arrows up in the air.
#' @param it see \code{\link{subset}}
#' @param is see \code{\link{subset}}
#' @param nx is the number of vectors/points along the x-axis
#' @param ny is the number of vectors/points along the y-axis
#' 
#' @examples
#'\dontrun{
#' sst <- sst.NCEP()
#' download.file('ftp://ftp.cdc.noaa.gov/Datasets/ncep.reanalysis.derived/surface/uwnd.mon.mean.nc','ncep_u.nc')
#' download.file('ftp://ftp.cdc.noaa.gov/Datasets/ncep.reanalysis.derived/surface/vwnd.mon.mean.nc','ncep_v.nc')
#' u <- retrieve('ncep_u.nc')
#' v <- retrieve('ncep_v.nc')
#' ## Make a map of SSt and wind vectors:
#' map(sst)
#' vec(u,v,a=1.5,length=0.03,new=FALSE)
#' 
#' ## Make a fancy map
#' map(sst,projection='sphere')
#' vec(u,v,a=0.5,r=1.1,length=0.03,projection='sphere',new=FALSE)
#'}
#' 
#' 
#' @export
vec <- function(x,y,it=NULL,a=NULL,r=1,ix=NULL,iy=NULL,new=TRUE,nx=150,ny=80,
                projection='lonlat',lonR=NULL,latR=NULL,axiR=0,verbose=FALSE,length=NULL,...) {
  if (verbose) print(paste('Add vector to a a map or sphere:',projection))
  if (!is.null(it)) {
    x <- subset(x,it=it); y <- subset(y,it=it)
  } else {
    if (verbose) print(class(x))
    x <- map(x,plot=FALSE,verbose=verbose)
    y <- map(y,plot=FALSE) 
    attr(x,'dimensions') <- dim(x)
  }
  d <- attr(x,'dimensions')
  if (is.null(d)) stop('vec needs esd field objects with proper attributes')
  if (verbose) {print(d); print(dim(x))}
  if (is.null(ix)) ix <- pretty(lon(x),n=nx)
  if (is.null(iy)) iy <- pretty(lat(x),n=ny)
  if (is.null(lonR)) lonR <- mean(lon(x),na.rm=TRUE)  # Logitudinal rotation
  if (is.null(latR)) latR <- mean(lat(y),na.rm=TRUE)  # Latitudinal rotation
  #print(c(d[2],d[1]))
  if (verbose) {print('---pretty coordinates: ---');print(ix); print(iy)}
  X <- coredata(x); Y <- coredata(y)
  dim(X) <- c(d[1],d[2])
  dim(Y) <- c(d[1],d[2])
  #X <- t(X); Y <- t(Y)
  x0 <- rep(ix,length(iy))
  y0 <- sort(rep(iy,length(ix)))
  ij <- is.element(ix,lon(x))
  ji <- is.element(iy,lat(x))
  if (verbose) {print(ix); print(lon(x)); print(sum(ij))
    print(iy); print(lat(x)); print(sum(ji))}
  dim(x0) <- c(length(ij),length(ji)); dim(y0) <- dim(x0)
  x0 <- x0[ij,ji]
  y0 <- y0[ij,ji]
  ii <- is.element(lon(x),ix)
  jj <- is.element(lat(x),iy)
  if (is.null(a)) a <- 0.25*max(diff(range(lon(x))),diff(range(lat(x))))/(max(d))
  if (is.null(length)) length= 0.5*a
  ## Set the magnitude of the vectors: use values from a & y
  x1 <- a*X[ii,jj]; y1 <- a*Y[ii,jj]
  #print(dim(x1)); print(c(length(x0),sum(ii),sum(jj)))
  x1 <- x0 + x1; y1 <- y0 + y1
  if (projection=='sphere') {
    if (verbose) print('o--- Projection on a sphere ---o')
    # Rotate data grid:
    # The coordinate of the arrow start:
    
    theta <- pi*x0/180; phi <- pi*y0/180
    x <- r*c(sin(theta)*cos(phi))
    y <- r*c(cos(theta)*cos(phi))
    z <- r*c(sin(phi))
    if (verbose) {print(dim(rotM(x=0,y=0,z=lonR))); print(dim(rbind(x,y,z))); print(c(lonR,latR))}
    A <- rotM(x=0,y=0,z=lonR) %*% rbind(x,y,z)
    A <- rotM(x=latR,y=0,z=0) %*% A
    if (verbose) print(dim(A))
    x0 <- A[1,]; y0 <- A[3,]
    invisible <- A[2,] < 0
    x0[invisible] <- NA; y0[invisible] <- NA
    #The coordinate of the arrow end:
    theta <- pi*x1/180; phi <- pi*y1/180
    x <- r*c(sin(theta)*cos(phi))
    y <- r*c(cos(theta)*cos(phi))
    z <- r*c(sin(phi))
    A <- rotM(x=0,y=0,z=lonR) %*% rbind(x,y,z)
    A <- rotM(x=latR,y=0,z=0) %*% A
    x1 <- A[1,]; y1 <- A[3,]
    invisible <- A[2,] < 0
    x1[invisible] <- NA; y1[invisible] <- NA
  }    
  
  if (verbose) {print('x:'); print(x0); print(x1); print('y:'); print(y0); print(y1)}
  if (new) {
    dev.new()
    plot(range(x0,x1,na.rm=TRUE),range(y0,y1,na.rm=TRUE),xlab='',ylab='')
    data(geoborders, envir = environment())
    lines(geoborders$x,geoborders$y)
  }
  arrows(x0, y0, x1, y1,length=length,...)
}
