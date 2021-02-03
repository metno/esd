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
  if (!is.null(it) | !is.null(is))
    X <- subset(X,it=it,is=is)
  
  if (na.action=='fill') {
    if (verbose) print('Fill missing data gaps')
    # Use interpolation to till missing data gaps. OK for small glitches.
    d <- dim(X)
    for (i in 1:d[2]) {
      x <- coredata(X[,i]); ok <- is.finite(x)
      X[,i] <- approx(y=x[ok],x=(1:d[1])[ok],xout=1:d[1])$y
    }
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

# Transfer PCA back to station data
#' @export
pca2station <- function(X,is=NULL,ip=NULL,anomaly=FALSE,what='pca',verbose=FALSE) {
  stopifnot(!missing(X), inherits(X,"pca"))
  if (inherits(X,'ds')) class(X) <- class(X)[-1]
  if (verbose) print('pca2station')
  
  pca <- X
  cls <- class(pca)
  
  ## REB 2016-11-03
  ## If there is only one single station, avoid collapse of dimension
  if (is.null(dim(attr(pca,'pattern')))) {
    dim(attr(pca,'pattern')) <- c(1,length(attr(pca,'pattern')))
  }
  
  U <- attr(pca,'pattern')
  W <- attr(pca,'eigenvalues')
  d <- dim(U)
  if (what=='pca') V <- coredata(pca) else if (what=='xval') {
    if (verbose) print('Retrieve the cross-validation data')
    V <- coredata(attr(X,'evaluation')[,seq(2,dim(attr(X,'evaluation'))[2],by=2)])
  } else if (what=='test') {
    if (verbose) print('Retrieve the cross-validation observations')
    V <- coredata(attr(X,'evaluation')[,seq(1,dim(attr(X,'evaluation'))[2]-1,by=2)])
  }
  V[!is.finite(V)] <- 0
  if (!is.null(ip)) {
    dU <- dim(U)
    dV <- dim(V)
    U <- U[,ip]
    W <- W[ip]
    V <- V[,ip]
    dim(U) <- c(dU[1], length(ip))
    dim(V) <- c(dV[1], length(ip))
  }

  if (verbose) {str(U); str(W); str(V)}
  if(length(W)>1) diag.W <- diag(W) else diag.W <- W
  x <- U %*% diag.W %*% t(V)
  if (verbose) str(x)
  
  if (verbose) print(paste('anomaly=',anomaly))
  if (!anomaly) {
    x <- x + c(attr(pca,'mean'))
  }
  x <- zoo(t(x),order.by=index(pca))
  
  if (anomaly) {
    attr(x,'aspect') <- 'anomaly' 
  } else {
    attr(x,'aspect') <- 'original'
  }
  #  attr(x,'call') <- match.call()
  #  class(x) <- class(pca)[-1]
  
  if (verbose) print(attr(pca,'location'))
  names(x) <- attr(pca,'location') # AM 30.07.2013 added
  
  x <- attrcp(attr(pca,'station'),x)
  
  # REB 2014-10-27: if the object is DS-results, then look for
  # cross-validation
  if (!is.null(attr(pca,'evaluation'))) {
    if (verbose) print('include evaluation')
    cval <- attr(pca,'evaluation')
    d.cval <- dim(cval)
    V.x <- coredata(cval)
    # The evaluation data are stored as the original calibration
    # predictand followed by the prediction, i.e. station 1,1,2,2,3,3
    if (is.null(ip)) {
      ii1 <- seq(1,d.cval[2]-1,by=2)
      ii2 <- seq(2,d.cval[2],by=2)
    } else {
      ii1 <- ip*2-1#seq(1,d.cval[2]-1,by=2) ip = [1,2,3]
      ii2 <- ip*2#seq(2,d.cval[2],by=2)
    }
    # Recover the station data from the original data x and the
    # cross-validation prediction z
    # seperately using the same spatial PCA pattern and eigenvalues:
    x.cvalx <- U %*% diag.W %*% t(V.x[,ii1])
    x.cvalz <- U %*% diag.W %*% t(V.x[,ii2])
    # Combine the two together and then sort so that the prediction
    # of the first station follows the observation
    # from the first station:
    ii <- order(seq(1,d.cval[1],by=1),seq(1,d.cval[1],by=1)+0.5)
    x.cval <- rbind(x.cvalx,x.cvalz)[,ii]
    mpca <- c(attr(pca,'mean'))
    jj <- order(c(1:length(mpca),1:length(mpca)+0.5))
    if (!anomaly) x.cval <- x.cval + rep(mpca,2)[jj]
    if (anomaly) attr(x.cval,'aspect') <- 'anomaly' else
      attr(x.cval,'aspect') <- 'original'
    attr(x,'evaluation') <- zoo(t(x.cval),order.by=index(cval)) 
  }
  
  ## REB 2016-05-09: if the object is DS-results, then look for
  ## common EOFs
  if (!is.null(attr(pca,'n.apps'))) {
    if (verbose) print(paste('include',attr(pca,'n.apps')))
    for (i.app in 1:attr(pca,'n.apps')) {
      cval <- attr(pca,paste('appendix.',i.app,sep=''))
      d.cval <- dim(cval)
      V.x <- coredata(cval)
      if(!is.null(ip)) {
        V.x <- V.x[,ip,]
      }
      x.cval <-U %*% diag.W %*% t(V.x)
      mpca <- c(attr(pca,'mean'))
      if (!anomaly) x.cval <- x.cval + mpca
      if (anomaly) attr(x.cval,'aspect') <- 'anomaly' else
        attr(x.cval,'aspect') <- 'original'
      attr(x,paste('appendix.',i.app,sep='')) <- zoo(t(x.cval),order.by=index(cval))
    }
  }
  
  #nattr <- softattr(pca)
  #for (i in 1:length(nattr))
  #  attr(x,nattr[i]) <- attr(pca,nattr[i])
  # Add meta data as attributes:
  attr(x,'longitude') <- attr(pca,'longitude')
  attr(x,'latitude') <- attr(pca,'latitude')
  attr(x,'history') <- history.stamp(pca)
  class(x) <- cls[-1]
  if(!is.null(is)) x <- subset(x, is=is)
  if (verbose) print(class(x))
  if (verbose) print("exit pca2station")
  invisible(x)
}
