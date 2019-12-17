# Documentation in map.R
#' @export
vec <- function(x,y,it=NULL,a=NULL,r=1,ix=NULL,iy=NULL,new=TRUE,nx=150,ny=80,
                projection='lonlat',lonR=NULL,latR=NULL,axiR=0,verbose=FALSE,length=NULL,...) {
  if (verbose) print(paste('Add vector to a a map or sphere:',projection))
  if (!is.null(it)) {x <- subset(x,it=it); y <- subset(y,it=it)} else
  {x <- map(x,plot=FALSE); y <- map(y,plot=FALSE); attr(x,'dimensions') <- dim(x)}
  d <- attr(x,'dimensions')
  if (is.null(d)) stop('vec needs esd field objects with proper attributes')
  if (verbose) {print(d); print(dim(x))}
  if (is.null(ix)) ix <- pretty(lon(x),n=nx)
  if (is.null(iy)) iy <- pretty(lat(x),n=ny)
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
  x1 <- a*X[ii,jj]; y1 <- a*Y[ii,jj]
  #print(dim(x1)); print(c(length(x0),sum(ii),sum(jj)))
  x1 <- x0 + x1; y1 <- y0 + y1
  if (projection=='sphere') {
    # Rotate data grid:
    # The coordinate of the arrow start:
    
    theta <- pi*x0/180; phi <- pi*y0/180
    x <- r*c(sin(theta)*cos(phi))
    y <- r*c(cos(theta)*cos(phi))
    z <- r*c(sin(phi))
    a <- rotM(x=0,y=0,z=lonR) %*% rbind(x,y,z)
    x0 <- a[1,]; y0 <- a[3,]
    invisible <- a[2,] < 0
    x0[invisible] <- NA; y0[invisible] <- NA
    #The coordinate of the arrow end:
    theta <- pi*x1/180; phi <- pi*y1/180
    x <- r*c(sin(theta)*cos(phi))
    y <- r*c(cos(theta)*cos(phi))
    z <- r*c(sin(phi))
    a <- rotM(x=0,y=0,z=lonR) %*% rbind(x,y,z)
    x1 <- a[1,]; y1 <- a[3,]
    invisible <- a[2,] < 0
    x1[invisible] <- NA; y1[invisible] <- NA
  }    
  
  if (verbose) {print('x:'); print(x0); print(x1); print('y:'); print(y0); print(y1)}
  if (new) {
    dev.new()
    plot(range(x0,x1),range(y0,y1),xlab='',ylab='')
    data(geoborders, envir = environment())
    lines(geoborders$x,geoborders$y)
  }
  arrows(x0, y0, x1, y1,length=length,...)
}
