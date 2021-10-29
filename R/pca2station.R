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
