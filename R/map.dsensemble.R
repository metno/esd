## Map downscaled ensembles (dse-objects)
## Map the result according to time (it), space (is) or member (im)
## Select a set of PCs and then use these in matrix product to reproduce
## physical elements.

#' Expand PCA to obtain station data
#'
#' @param x an object of type 'pca'
#' @param it time index (see \code{\link{subset}})
#' @param FUN function applied to aggregate in time
#' @param FUNX function applied aggregate ensemble members
#' @param verbose if TRUE print progress
#' @param anomaly if FALSE add the mean value stored as attribute in x
#' @param test if TRUE perform test on one GCM simulation
#'
#' @export
expandpca <- function(x,it=NULL,FUN=NULL,FUNX='mean',verbose=FALSE,anomaly=FALSE,test=FALSE,eof=TRUE) {
  ## Get the spatial weights
  if (verbose) print('expandpca')
  if ((eof) & (!is.null(x$eof))) x$pca <- x$eof
  if (test) print('--TEST ON ONE GCM simulation--')
  if (inherits(x,'pca')) UWD <- x$pca else UWD <- x$eof
  if (!is.null(FUN)) {
    if (FUN != 'mean') anomaly <- TRUE; 
    if (verbose) print(c(FUN,anomaly))
    }
  if (verbose) print(names(attributes(UWD)))
  ## Eigenvalues
  D <- attr(UWD,'eigenvalues')
  ## Create a matrix with only the GCM time series
  if (verbose) print('PCA/EOF-based ensemble')
  X <- x
  X$info <- NULL; X$pca <- NULL; X$eof <- NULL
  if (verbose) for (ii in 1:length(X)) print(dim(X[[ii]]))
  ## Dimension of each downscaled GCM results
  if (verbose) print(paste('subset.pc, it=',it))
  ## Check if the ensemble members have the same size - if not, only keep the ones with most common sizes
  if (verbose) print('Check ensemble member size')
  n <- length(names(X))
  if (verbose) print(paste('Original length of X is',n))
  memsiz <- rep("?",n)
  for (i in 1:n) memsiz[i] <- paste(dim(X[[i]]),collapse='x')
  memsiztab <- table(memsiz)
  if (verbose) print(memsiztab)
  memkeep <- rownames( memsiztab)[as.numeric(memsiztab)==max(as.numeric(memsiztab))]
  if (verbose) print(memkeep)
  im <- sort((1:n)[-grep(memkeep,memsiz)],decreasing = TRUE)
  if (verbose) print(im)
  for (ix in im) X[[ix]] <- NULL
  n <- length(names(X))
  if (verbose) print(paste('New length of X is',n))
  V <- lapply(X,FUN='subset.pc',it=it)
  memsiz <- rep("?",n)
  for (i in 1:n) memsiz[i] <- paste(dim(V[[i]]),collapse='x')
  memsiztab <- table(memsiz)
  if (verbose) print(memsiztab)
  d <- dim(V[[1]])
  
  if (verbose) {print(names(V)); print(c(d,n,length(unlist(V)))); print(paste('FUNX=',FUNX))}
  ## Apply function FUNX
  if (!test) {
    V <- unlist(V)
    dim(V) <- c(d[1]*d[2],n)
    if (verbose) print(FUNX)
    V <- apply(V,1,FUN=FUNX)
  } else {
    V <- V[[1]] # Pick one member for testing ## Testing
    n <- 1
    d <- dim(V)
  }
  if (verbose) {print('FUNX done'); print(c(n,d))}
  dim(V) <- d
  ## REB 2016-12-01: Can also aggregate in time to speed things up and create a vector  
  if (!is.null(FUN)) {  
    if (verbose) print(paste('FUN=',FUN))
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
  if (!is.null(dim(U))) {
    dU <- dim(U) 
  } else {
    dU <- c(1,length(U)) # If there is only one single station
    dim(U) <- dU
  }
                      
  if (verbose) {print(d); print(dU)}
  if (inherits(UWD,'eof')) {
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
  # Not right if FUN is defined and time mean has been applied:
  if(nrow(V)==length(index(subset(X[[1]],it=it)))) {
    Y <- zoo(Y,order.by=index(subset(X[[1]],it=it)))
  } else {
    Y <- zoo(Y,order.by=seq(nrow(V)))
  }
  Y <- attrcp(UWD,Y)
  attr(Y,'time') <- range(index(subset(X[[1]],it=it)))
  class(Y) <- class(UWD)[-1]
  if (inherits(UWD,'eof')) {
    if (verbose) print('Use dimensions and lon/lat from EOFs')
    attr(Y,'dimensions') <- c(attr(x$eof,'dimensions')[1:2],length(index(V)))
    attr(Y,'longitude') <- lon(UWD)
    attr(Y,'latidude') <- lat(UWD)
  }
  attr(Y,'mean') <- NULL
  if (verbose) {print('exit expandpca'); print(dim(Y))}
  return(Y)
}

#' @exportS3Method
#' @export
map.dsensemble <- function(x,it=c(2000,2099),is=NULL,im=NULL,ip=NULL,
                           colbar=list(pal=NULL,rev=FALSE,n=10,breaks=NULL,pos=0.05,
                                   show=TRUE,type="p",cex=2,h=0.6,v=1),
                           FUN='mean',FUNX='mean',verbose=FALSE,anomaly=FALSE,test=FALSE,plot=TRUE,...) {
  ## PCA/EOF objects

  if (verbose) print('map.dsensemble')

  if (inherits(x,c('pca','eof'))) {
    ## Extract a subset of the data
    if (verbose) print(names(x)[2])
    x <- subset(x,is=is,im=im,ip=ip,verbose=TRUE)#verbose)
    ## REB 2016-12-01: Do all the analysis on the PC weights to speed up. Linearity.  
#    Y <- expandpca(x,it=it,FUNX=FUNX,verbose=verbose,anomaly=anomaly,test=test)
    Y <- expandpca(x,it=it,FUN=FUN,FUNX=FUNX,verbose=verbose,anomaly=anomaly,test=test)
    if (verbose) {str(x[[2]]); str(Y)}
#    if (plot) map(Y,FUN=FUN,colbar=colbar,verbose=verbose,...)
    if (plot) map(Y,FUN="mean",colbar=colbar,verbose=verbose,...)
    invisible(Y)
  } else return(NULL)
}

## Function for extracting the subset from PCs stored as zoo
# not exported
subset.pc <- function(x,ip=NULL,it=NULL,verbose=FALSE) {
  if (verbose) print('subset.pc')
  d <- dim(x)
  if (!is.null(it)) {
    if (verbose) print('subset it')
    if ((is.numeric(it) | is.integer(it)) & is.dates(index(x))) {
        it <- c(as.Date(paste(it,'01-01',sep='-')),
                as.Date(paste(it,'12-31',sep='-')))
    }
    x <- window(x,start=min(it),end=max(it))
  }
  if (!is.null(ip)) {
    if (verbose) print('subset pattern')
    x <- x[,ip]
    d <- dim(x)
  }
  dim(x) <- c(length(index(x)),d[2])
  if(verbose) print(dim(x))
  return(x)
}





