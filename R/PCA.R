# Documentation in EOF.R
#' @export
PCA <- function(X,...) UseMethod("PCA")

#' @exportS3Method
#' @export
PCA.default <- function(X,...) {
  stop("Don't know how to handle objects other than station")
}

#' @exportS3Method
#' @export
PCA.matrix <- function(X,...,verbose=FALSE) {
  if(verbose) print("PCA.matrix")
  if(length(dim(X))==3) {
    product <- function(x) x[1]*x[2]
    Z <- unlist(X)
    dim(Z) <- c(product(dim(X)[1:2]),dim(X)[3])
  }
  if(length(dim(X))!=2) {
    if(verbose) print("Warning! Wrong input dimensions. Input should be a matrix of 2 or 3 dimensions.")
    Z.pca <- NULL
  } else {
    Z.pca <- svd(Z)
  }
  return(Z.pca)
}


#' @exportS3Method
#' @export
PCA.station <- function(X,...,n=20,na.action='fill',verbose=FALSE,it=NULL,is=NULL,anomaly=TRUE) {
  if (verbose) print('PCA.station')
  if (!is.null(it) | !is.null(is))
    X <- subset(X,it=it,is=is)
  
  if (na.action=='fill') {
    if (verbose) print('Fill missing data gaps')
    # Use interpolation to till missing data gaps. OK for small glitches.
    d <- dim(X)
    for (i in 1:d[2]) {
      x <- coredata(X[,i]); ok <- is.finite(x)
      if (sum(ok)>0) X[,i] <- approx(y=x[ok],x=(1:d[1])[ok],xout=1:d[1])$y else
                     X[,i] <- NA
    }
    if (sum(!is.finite(X))>0) stop('PCA.station: detected invalid data (NA, etc)')
  }
  
  # Re-order dimensions: space x time
  x <- t(coredata(X))
  neofs <- n
  D <- dim(x)
  neofs <- min(neofs,D[1])
  ok.time <- is.finite(colMeans(x))
  z <- x[,ok.time]
  ok.site <- is.finite(rowMeans(z))
  z <- z[ok.site,]
  if (verbose) print(paste('Extract',sum(ok.site),'stations',
                           'and',sum(ok.time),'time steps',
                           'dimension=',dim(z)[1],'x',dim(z)[2]))
  X.clim <- rowMeans(x,na.rm=TRUE)
  pca <- svd(z - rowMeans(z))
  #str(pca); print(neofs)
  
  autocor <- 0
  #print(paste("Find max autocorr in the grid boxes."))
  #print(dim(y))
  for (i in 1:dim(z)[2]) {
    vec <- as.vector(coredata(z[,i]))
    ar1 <- acf(vec[],plot=FALSE)
    autocor <- max(c(autocor,ar1$acf[2,1,1]),na.rm=TRUE)
  }  
  
  # Recover the original matrix size, and insert the results
  # where there was valid data
  U <- matrix(rep(NA,D[1]*neofs),D[1],neofs)
  U[ok.site,] <- pca$u[,1:neofs]
  V <- matrix(rep(NA,D[2]*neofs),D[2],neofs)
  #print(D); print(dim(V)); print(dim(pca$v))
  #print(dim(V[ok.time,])); print(dim(pca$v[1:neofs,]))
  V[ok.time,] <- pca$v[,1:neofs]
  y <- zoo(V,order.by=index(X))
  names(y) <- paste("X.",1:neofs,sep="")
  
  invert <- apply(U,2,mean) < 0
  U[,invert] <- -U[,invert]
  y[,invert] <- -y[,invert]
  
  y <- attrcp(X,y)
  #nattr <- softattr(X)
  #for (i in 1:length(nattr))
  #  attr(y,nattr[i]) <- attr(X,nattr[i])
  attr(y,'pattern') <- U
  attr(y,'dimensions') <- D
  attr(y,'mean') <- X.clim
  attr(y,'max.autocor') <- autocor
  attr(y,'eigenvalues') <- pca$d[1:neofs]
  attr(y,'sum.eigenv') <- sum(pca$d)
  attr(y,'tot.var') <- sum(pca$d^2)
  attr(y,'aspect') <- 'anomaly'
  attr(y,'history') <- history.stamp(X)
  class(y) <- c("pca",class(X))
  invisible(y)
}
