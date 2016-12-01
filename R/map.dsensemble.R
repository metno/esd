## Map downscaled ensembles (dse-objects)
## Map the result according to time (it), space (is) or member (im)
## Select a set of PCs and then use these in matrix product to reproduce
## physical elements.

expandpca <- function(x,it=NULL,FUN=NULL,FUNX='mean',verbose=FALSE,anomaly=FALSE,test=FALSE) {
  ## Get the spatial weights
  if (verbose) print('expandpca')
  if (test) print('--TEST ON ONE GCM simulation--')
  if (inherits(x,'pca')) UWD <- x$pca else UWD <- x$eof
  if (!is.null(FUN)) {
    if (FUN != 'mean') anomaly <- TRUE; 
    if (verbose) print(c(FUN,anomaly))
    }
  if (verbose) print(names(attributes(UWD)))
  D <- attr(UWD,'eigenvalues')
  ## Create a matrix with only the GCM time series
  if (verbose) print('PCA/EOF-based ensemble')
  X <- x
  X$info <- NULL; X$pca <- NULL; X$eof <- NULL
  V <- lapply(X,FUN='subset.pc',it=it)
  if (!test) {
    n <- length(names(V))
    d <- dim(V[[1]])
    V <- unlist(V)
    dim(V) <- c(d[1]*d[2],n)
    V <- apply(V,1,FUN=FUNX)
  } else {
    V <- V[[1]] # Pick one member for testing ## Testing
    n <- 1
    d <- dim(V)
  }
  if (verbose) print(c(n,d))
  dim(V) <- d
  
  ## REB 2016-12-01: Can also aggregate in time to speed things up and create a vector  
  if (!is.null(FUN)) {  
      if (FUN=='trend') FUN <- 'trend.coef'
      if (verbose) print(paste('FUN=',FUN,!is.null(dim(V))))
      if (is.null(dim(V))) dim(V) <- c(1,length(V))
      V <- apply(V,2,FUN=FUN)
      
      if (!is.null(dim(V))) d <- dim(V) else {
                         d <- c(1,length(V)) # If there is only one single time point
                         dim(V) <- d
      }
      if (verbose) print(V)
  }
    
  ## Aggregate statistics over ensemble members
  if (verbose) print('Aggregate ensemble statistics')
  ## Apply FUNX to each of the PCs across all members
  #
  
  U <- attr(UWD,'pattern')
  if (!is.null(dim(U))) dU <- dim(U) else {
                        dU <- c(1,length(U)) # If there is only one single station
                        dim(U) <- dU
                      }
  if (verbose) {print(d); print(dU)}
  if (inherits(x,'eof')) {
    if (verbose) {print('eof'); print(dU)}
    dim(U) <- c(dU[1]*dU[2],dU[3])
  }
  dim(V) <- d
  if (verbose) {
    print('Matrix multiplication')
    str(U); str(D); str(V)
  }
  Y <- V %*% diag(D) %*% t(U)
  ## Add mean and insert into zoo frame
  if (!anomaly) {
    Y <- t(t(Y) + c(attr(UWD,'mean')))
    if (verbose) print('add mean field')
  }
  Y <- zoo(Y,order.by=index(subset(X[[1]],it=it)))
  Y <- attrcp(UWD,Y)
  class(Y) <- class(UWD)[-1]
  if (inherits(x,'eof')) attr(Y,'dimensions') <- c(attr(x$eof,'dimensions')[1:2],length(index(V)))
  attr(Y,'mean') <- NULL
  if (verbose) {print('exit expandpca'); print(dim(Y))}
  return(Y)
}

## Function for extracting the subset from PCs stored as zoo
subset.pc <- function(x,ip=NULL,it=NULL,verbose=FALSE) {
  if (verbose) print('subset.pc')
  d <- dim(x)
  if (!is.null(it)) {
    if (verbose) print('subset it')
    if (is.numeric(it) | is.integer(it)) 
      it <- c(as.Date(paste(it,'01-01',sep='-')),
              as.Date(paste(it,'12-31',sep='-')))
    x <- window(x,start=min(it),end=max(it))
  }
  if (!is.null(ip)) {
    if (verbose) print('subset pattern')
    x <- x[,ip]
    d <- dim(x)
  }
  dim(x) <- c(length(index(x)),d[2])
  if (verbose) print(dim(x))
  return(x)
}



map.dsensemble <- function(x,it=c(2000,2099),is=NULL,im=NULL,ip=NULL,
                           colbar=list(pal=NULL,rev=FALSE,n=10,breaks=NULL,pos=0.05,
                                   show=TRUE,type="p",cex=2,h=0.6,v=1),
                           FUN='mean',FUNX='mean',verbose=FALSE,anomaly=FALSE,test=FALSE,plot=TRUE,...) {
  ## PCA/EOF objects

  if (verbose) print('map.dsensemble')
  
  if (inherits(x,c('pca','eof'))) {
    ## Extract a subset of the data
    if (verbose) print(names(x)[2])
      x <- subset(x,is=is,im=im,ip=ip,verbose=verbose)
    ## REB 2016-12-01: Do all the analysis on the PC weights to speed up. Linearity.  
#    Y <- expandpca(x,it=it,FUNX=FUNX,verbose=verbose,anomaly=anomaly,test=test)
    Y <- expandpca(x,it=it,FUN=FUN,FUNX=FUNX,verbose=verbose,anomaly=anomaly,test=test)
    if (verbose) {str(x[[2]]); str(Y)}
#    if (plot) map(Y,FUN=FUN,colbar=colbar,verbose=verbose,...)
    if (plot) map(Y,FUN="mean",colbar=colbar,verbose=verbose,...)
    invisible(Y)
  } else return(NULL)
}



## Tools to subset or reduce the size of a dsensemble, e.g. removing the
## high-order modes of PCA/EOF that represent noise.
subset.dsensemble.multi <- function(x,ip=NULL,it=NULL,is=NULL,im=NULL,
                              verbose=FALSE,...) {
 
  if (verbose) print('subset.dsensemble.multi')
  cls <- class(x)
  
  Y <- list()
  Y$info <- x$info
  if (inherits(x,'pca')) {
    if (verbose) print('subset pca')
    Y$pca <- subset(x$pca,it=it,is=is,ip=ip,verbose=verbose)
  }
  if (inherits(x,'eof')) {
    if (verbose) print('subset eof')
    Y$eof <- subset(x$eof,it=it,is=is,ip=ip,verbose=verbose)
  }
  X <- x

  X$info <- NULL; X$pca <- NULL; X$eof <- NULL
  n <- length(names(X))
  if (verbose) print('subset gcm-zoo')
  y <- lapply(X,FUN='subset.pc',ip=ip,it=it)
  if (verbose) print(dim(y[[1]]))

  if (!is.null(im)) {
    ## Subset ensemble members
    if(verbose) print(paste('subset im',length(y)))
    if (is.logical(im)) im <- (1:n)[im]
    for (i in rev(setdiff(1:n,im))) y[[i]] <- NULL
    if(verbose) print(paste('subset im',length(y)))
  }
  Y <- c(Y,y)
  class(Y) <- cls
  if (verbose) print('exit subset.dsensemble.multi')
  return(Y)
}


