## Map downscaled ensembles (dse-objects)
## Map the result according to time (it), space (is) or member (im)
## Select a set of PCs and then use these in matrix product to reproduce
## physical elements.

expandpca <- function(x,it=NULL,FUNX='mean',verbose=FALSE,anomaly=FALSE,test=FALSE) {
  ## Get the spatial weights
  if (verbose) print('expandpca')
  if (test) print('--TEST ON ONE GCM simulation--')
  if (inherits(x,'pca')) UWD <- x$pca else UWD <- x$eof
  if (verbose) print(names(attributes(UWD)))
  D <- attr(UWD,'eigenvalues')
  ## Create a matrix with only the GCM time series
  if (verbose) print('PCA/EOF-based ensemble')
  X <- x
  X$info <- NULL; X$pca <- NULL; X$eof <- NULL
  V <- lapply(X,FUN='subsetzoo',it=it)
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
  
  ## Aggregate statistics over ensemble members
  if (verbose) print('Aggregate ensemble statistics')
  ## Apply FUNX to each of the PCs across all members
  #
  U <- attr(UWD,'pattern'); dU <- dim(U)
  if (length(dU==3)) {
    dim(U) <- c(dU[1]*dU[2],dU[3])
  }
  dim(V) <- d
  if (verbose) {
    print('Matrix multiplication')
    str(U); str(D); str(V)
  }
  Y <- V %*% diag(D) %*% t(U)
  ## Add mean and insert into zoo frame
  if (!anomaly) Y <- t(t(Y) + c(attr(UWD,'mean')))
  Y <- zoo(Y,order.by=index(V))
  Y <- attrcp(UWD,Y)
  class(Y) <- class(UWD)[-1]
  attr(Y,'mean') <- NULL
  if (verbose) print('expandpca done')
  return(Y)
}

## Function for extracting the subset from PCs stored as zoo
subsetzoo <- function(x,pattern=NULL,it=NULL,verbose=FALSE) {
  if (verbose) print('subsetzoo')
  if (!is.null(it)) {
    if (verbose) print('subset it')
    if (is.numeric(it) | is.integer(it)) 
      it <- as.Date(paste(it,'01-01',sep='-'))
    x <- window(x,start=min(it),end=max(it))
  }
  if (!is.null(pattern)) {
    if (verbose) print('subset pattern')
    x <- x[,pattern]
  }
  return(x)
}



map.dsensemble <- function(x,it=c(2000,2099),is=NULL,im=NULL,pattern=NULL,colbar=NULL,
                           FUN='mean',FUNX='mean',verbose=FALSE,anomaly=FALSE,test=FALSE) {
  ## PCA/EOF objects

  if (verbose) print('map.dsensemble')
  
  if (inherits(x,c('pca','eof'))) {
    ## Extract a subset of the data
    if (verbose) print(names(x)[2])
    x <- subset(x,is=is,im=im,pattern=pattern,verbose=verbose)
    Y <- expandpca(x,it=it,FUNX=FUNX,verbose=verbose,anomaly=anomaly,test=test)
    
    map(Y,FUN=FUN,colbar=colbar,verbose=verbose)
    invisible(Y)
  } else return(NULL)
}



## Tools to subset or reduce the size of a dsensemble, e.g. removing the
## high-order modes of PCA/EOF that represent noise.
subset.dsensemble.multi <- function(x,pattern=NULL,it=NULL,is=NULL,im=NULL,
                              verbose=FALSE,...) {
 
  if (verbose) print('subset.dsensemble.multi')
  cls <- class(x)
  
  Y <- list()
  Y$info <- x$info
  if (inherits(x,'pca')) {
    if (verbose) print('subset pca')
    Y$pca <- subset(x$pca,it=it,is=is,pattern=pattern,verbose=verbose)
  }
  if (inherits(x,'eof')) {
    if (verbose) print('subset eof')
    Y$eof <- subset(x$eof,it=it,is=is,pattern=pattern,verbose=verbose)
  }
  X <- x

  X$info <- NULL; X$pca <- NULL; X$eof <- NULL
  n <- length(names(X))
  if (verbose) print('subset gcm-zoo')
  y <- lapply(X,FUN='subsetzoo',pattern=pattern,it=it)
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
  return(Y)
}


