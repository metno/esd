# internal function - no need to export
regrid.irregweights <- function(xo,yo,xn,yn,verbose=FALSE) {
# Compute the weights for irregular grids (xo,yo) - the station class

  mindist <- function(x,X,verbose=FALSE) {
    d <- distAB(x[1],x[2],X[1,],X[2,])
    if (min(d) > 0) {
      gridind <- order(d)[1:4]
      weight <- 1/d[gridind]
      weight <- weight/sum(weight)
      nearest4 <- c(weight,gridind)
    } else {
      n <- length(d)
      nearest4 <- c(1,rep(0,3),(1:n)[d==0][1],rep(1,3))
    }
    names(nearest4) <- c(rep('weight',4),rep('index',4))
    return(nearest4)
  }
  
  if (verbose) print('regrid.irregweights')
  t1 <- Sys.time()
  nx <- length(xn); ny <- length(yn)
  xy <- rbind(xo,yo)
  XY <- rbind(rep(xn,ny),sort(rep(yn,nx)))
  dx <- apply(XY,2,'mindist',xy)
  beta <- t(dx[1:4,])
  attr(beta,'index') <- t(dx[5:8,])
  
  t2 <- Sys.time()
  if (verbose) {print("Computation time:"); print(t2 - t1)}
  invisible(beta)
}


#' @exportS3Method
#' @export
regrid.station <- function(x,is=NULL,...,approach="station",verbose=FALSE) {

  stopifnot(inherits(x,'station'))
  if (verbose) print('regrid.station')  
  if (approach=="pca2station") {
    y <- regrid.pca2station(x,is=is)
    return(y)
  }
  

  
  # The coordinates which to grid: lon.new & lat.new
  if ( (is.data.frame(is)) | (is.list(is)) )
    {lon.new <- is[[1]]; lat.new <- is[[2]]} else
  if ( (inherits(is,'station')) | (inherits(is,'field')) |
      (inherits(is,'eof')) | (inherits(is,'pca')) ) {
    lon.new <- attr(is,'longitude'); lat.new <- attr(is,'latitude')
  }
  
  
  if (verbose) {print("New coordinates"); print(lon.new); print(lat.new)}
  if (verbose) {print("Longitude range"); print(range(lon.new)); print(range(lon(x)))}
  if (verbose) {print("Latitude range"); print(range(lat.new)); print(range(lat(x)))}
  
  mlon <- MATCH(attr(x,'longitude'), lon.new)
  mlat <- MATCH(attr(x,'latitude'),  lat.new)
 
  if ( sum(is.na(c(mlon,mlat))) ==0 ) {
    if (verbose) print(summary(mlon,mlat))
     if ( max( c(diff(mlon),diff(mlat)) ) == 1) return(x)
  }
  
  ns <- length(attr(x,'longitude'))
  if (verbose) print(paste("regrid.station: from",
                           ns,"locations to", length(lon.new),"x",length(lat.new),'regular grid'))  
  
# Apply the grid-based subset - only use the nearest longitudes/latitudes
#  for the regridding. 

# Clever solution (not default): fix edge problems...  
 
  lon.old <- attr(x,'longitude')
  lat.old <- attr(x,'latitude')
 
  #print("...");print(lonj2w)
  #print(lon.old); print(lat.old)
  if (verbose) print(paste("subsample the original grid and use",length(lon.old),
                      "longitudes and",length(lat.old),"latitudes"))

  beta <- regrid.irregweights(lon.old,lat.old,lon.new,lat.new,verbose=verbose)
                                          
  d <- dim(x)
  X <- x; 
  D <- c(length(lon.new),length(lat.new))
  y <- matrix(rep(NA,D[1]*D[2]*d[1]),D[1]*D[2],d[1])
  if (verbose) {
    print("Weight matrix");  print(dim(beta))
    print('input data dimensions'); print(c(D,D[1]*D[2]))
    print('output data dimensions'); print(dim(y))
    print('index:'); print(dim(attr(beta,'index')))
    print(range(c(attr(beta,'index'))))
  }
  
  pb <- txtProgressBar(style=3)
  #str(X); str(beta)
  for (i in 1:d[1]) {
    setTxtProgressBar(pb,i/d[1])  
    M <- as.matrix(X[i,attr(beta,"index")])
    dim(M) <- c(D[1] * D[2],4)
    y[, i] <- rowSums(as.matrix(beta) * M)      
  }
  
  print("set attributes:")
  y <- zoo(t(y),order.by=index(x))
  class(y) <- c('field',class(x)[2:3]) 

  y <- attrcp(x,y)
  attr(y,'unit') <- unit(y)[1]
  attr(y,'variable') <- varid(y)[1]
  attr(y,'longitude') <- lon.new
  attr(y,'latitude') <- lat.new
  attr(y,'altitude') <- lat.new
  attr(y,'location') <- NULL
  attr(y,'dimensions') <- c(D,d[1])
  if (verbose) print(paste("New dimensions:",attr(y,'dimensions')[1],"x",
             attr(y,'dimensions')[2],"x",attr(y,'dimensions')[3] ))

  attr(y,'history') <- history.stamp(x)
 
  invisible(y)
}

#' @exportS3Method
#' @export
regrid.pca <- function(x,is=NULL,...,verbose=FALSE) {
  stopifnot(inherits(x,'pca'))
  if (is.list(is)) {
    lon.new <- is[[1]]
    lat.new <- is[[2]]
  } else if ( (inherits(is,'station')) | (inherits(is,'field')) |
              (inherits(is,'eof')) ) {
    lon.new <- attr(is,'longitude')
    lat.new <- attr(is,'latitude')
  }
  lon.old <- attr(x,'longitude')
  lat.old <- attr(x,'latitude')

  greenwich <- attr(x,'greenwich')
  if ( greenwich & (min(lon.new < 0)) ) x <- g2dl(x,greenwich=FALSE)
  
  if (sum(is.na(c( MATCH(lon.old, lon.new),
                   MATCH(lat.old, lat.new) ))) == 0) return(x)
  
  if (verbose) print(paste("regrid.field: from",length(lon.old),"x",
              length(lat.old),"to",
              length(lon.new),"x",length(lat.new)))  

  beta <- regrid.irregweights(lon.old,lat.old,lon.new,lat.new)
                                          
  #print(dim(beta))
  X <- attr(x,'pattern')
  d <- dim(X)
  dim(X) <- c(d[1]*d[2],d[3])
  D <- c(length(lon.new),length(lat.new))
  Y <- matrix(rep(NA,D[1]*D[2]*d[3]),D[1]*D[2],d[3])
  for (i in 1:d[3]) {
    cat(".")
#    z <- apply(cbind(beta,attr(beta,'index')),1,sparseMproduct,coredata(X[i,]))
#   print(c(i,d[1],length(z),length(y[,i]),NA,dim(X),dim(y)))
#    Y[,i] <- z
 #   old lines
 #z <- apply(b, 1, sparseMproduct, X[, i])
 #print(c(i, d[1], length(z), length(Y[, i]), NA, dim(X), dim(Y)))
 #Y[, i] <- z
 # New lines
    M <- as.matrix(X[attr(beta,"index"),i])
    dim(M) <- c(D[1] * D[2],4)
    Y[, i] <- rowSums(as.matrix(beta) * M)      
    #print(c(i, d[1], length(z), length(Y[, i]), NA, dim(X), dim(Y)))
    #Y[, i] <- z
  }
  #print("set attributes:")

  y <- x
  if (attr(y,'greenwich') != greenwich) y <- g2dl(y,greenwich)
  #mostattributes(y) <- attributes(x)
  #nattr <- softattr(x,ignore=c('longitude','latitude','dimensions'))
  #for (i in 1:length(nattr))
  #  attr(y,nattr[i]) <- attr(x,nattr[i])
  #dim(Y) <- c(D[1],D[2],d[3])
  y <- attrcp(x,y)
  Y ->  attr(y,'pattern')
  attr(y,'longitude') <- lon.new
  attr(y,'latitude') <- lat.new
  attr(y,'dimensions') <- c(D,d[3])
  if (verbose) print(paste("New dimensions:",attr(y,'dimensions')))
  #attr(y,'history') <- c(attr(y,'history'),'regrid.field')
  #attr(y,'date-stamp') <- date()
  #attr(y,'call') <- match.call()
  attr(y,'history') <- history.stamp(x)
  invisible(y)
}


regrid.pca2station <- function(x,is=NULL) {
  stopifnot(inherits(x,'station'),inherits(is,'list'))
  print("regrid.pca2station - estimate EOFs")

 #greenwich <- attr(x,'greenwich')
  if ( greenwich & (min(lon.new < 0)) ) x <- g2dl(x,greenwich=FALSE)
  pca0 <- PCA(x)
  pca1 <- regrid(pca0,is=is)
  print("Reconstruct field from EOFs")
  y <- pca2station(pca1)
  if (attr(y,'greenwich') != greenwich) y <- g2dl(y,greenwich)
  #attr(y,'history') <- c(attr(X,'history'),'regrid.eof2field')
  #attr(y,'date-stamp') <- date()
  #attr(y,'call') <- match.call()
  attr(y,'history') <- history.stamp(x)
  return(y)
}
