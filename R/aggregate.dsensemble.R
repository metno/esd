## Aggregate ensemble performs similar function to expandpca but applies the functions to the 
## results after the transform from PCA to field and hence is a bit slower.

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

#' aggregate.dsensemble
#' 
#' The function \code{aggregate.dsensemble} is used to extract and aggregate information from PCA-based dsensemble-objects, 
#' to expand PCA into station data or EOFs into field data.
#' It applies operations through the argument 'FUNX' to estimate quantities such as the ensemble mean and 
#' takes care of extra house keeping, such as attributes with meta data. \code{aggregate.dsensemble}
#' is similar to \code{expandpca} but applies 'FUNX' for each time step after it has been expanded by matrix multiplication 
#' of the singular vectors (X = U %*% diag(W) %*% t(V)) and is therefore slower. 
#' \code{map.dsensemble} is a wrapper for \code{map} that uses one of the two former routines to distill selected
#' information.
#' @seealso aggregate expandpca map.dsensemble
#' @aliases aggregate.dsensemble
#'
#' @param x an object of type 'dsensemble'
#' @param it time index (see \code{\link{subset}})
#' @param FUN function applied to aggregate in time
#' @param FUNX function applied aggregate ensemble members
#' @param verbose if TRUE print progress
#' @param anomaly if FALSE add the mean value stored as attribute in x
#' @param test if TRUE perform test on one GCM simulation
#' @param eof if TRUE and if the dsensemble object contains a spatially aggregated pattern, x$eof, 
#' expand with x$eof into an ensemble of field objects to which FUNX is applied. 
#' If FALSE, expand using the original pattern x$pca into an ensemble of station objects to be aggregated. 
#'
#' @exportS3Method
#' @export aggregate.dsensemble
aggregate.dsensemble <- function(x,...,it=NULL,FUN=NULL,FUNX='mean',verbose=FALSE,anomaly=FALSE,test=FALSE,eof=TRUE) {
  ## Get the spatial weights
  if (verbose) print('aggregate.ensemble')
  if ((eof) & (!is.null(x$eof))) x$pca <- x$eof
  ##x$eof <- gridmap(x$pca)
  if (test) print('--TEST ON ONE GCM simulation--')
  if (inherits(x,'pca')) UWD <- x$pca else UWD <- x$eof
  if (!is.null(FUNX)) {
    if (FUNX != 'mean') anomaly <- TRUE; 
    if (verbose) print(c(FUNX,anomaly))
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
  ## Only select the selected time interval - saves time
  if (is.null(it)) it <- range(index(X[[1]]))
  V <- lapply(X,FUN='subset.pc',it=it)
  if (verbose) print(paste('Interval',paste(range(index(V[[1]])),collapse=' - ')))
  memsiz <- rep("?",n)
  for (i in 1:n) memsiz[i] <- paste(dim(V[[i]]),collapse='x')
  memsiztab <- table(memsiz)
  if (verbose) print(memsiztab)
  d <- dim(V[[1]])
  gcmnames <- attr(X, "model_id")
  ## Quality control
  if (verbose) print(paste('Before quality control: original number of members=',n))
  for (i in seq(n,1,by=-1)) {
    #print(range(V[[i]],na.rm=TRUE)); print(dim(V[[i]]))
    if (max(abs(V[[i]]),na.rm=TRUE) > 10)  {
      print(paste(i,'Remove suspect results: ',gcmnames[i]))
      V[[i]] <- NULL
      gcmnames <- gcmnames[-i]
    }
  }
  n <- length(V)
  if (verbose) print(paste('After quality control: new number of members=',n))
  if (verbose) {print(names(V)); print(c(d,n,length(unlist(V)))); print(paste('FUNX=',FUNX))}
  lengths <- rep(NA,length(V))
  for (i in 1:length(V)) lengths[i] <- length(index(V[[i]]))
  if (length(table(lengths))>1) browser()
  
  ## Aggregate statistics over ensemble members
  if (verbose) print('Aggregate ensemble statistics')
  U <- attr(UWD,'pattern')
  if (!is.null(dim(U))) {
    dU <- dim(U) 
  } else {
    dU <- c(1,1,length(U)) # If there is only one single station
    dim(U) <- dU
  }
  
  if (verbose) {print(d); print(dU)}
  if (inherits(UWD,'eof')) {
    if (verbose) {print('eof'); print(dU)}
    dim(U) <- c(dU[1]*dU[2],dU[3])
  }
  if (verbose) {
    print('Matrix multiplication')
    str(U); str(D); str(V)
  }
  ## Loop through each time step - aggregate ensemble statistics for each time step
  Y <- matrix(rep(NA,dU[1]*dU[2]*d[1]),d[1],dU[1]*dU[2])
  for (it in 1:d[1]) { 
    ## loop through each ensemble member
    z <- matrix(rep(NA,dU[1]*dU[2]*n),dU[1]*dU[2],n)
    for (im in 1:n) { 
      v <- V[[im]]
      z[,im] <- v[it,] %*% diag(D) %*% t(U)
    }
    Y[it,] <- apply(z,1,FUNX)
  }
  ## Add mean and insert into zoo frame
  if (!anomaly) {
    if (verbose) print('add mean field')
    Y <- t(t(Y) + c(attr(UWD,'mean')))
  }
  # Not right if FUNX is defined and time mean has been applied:
  if(nrow(V[[1]])==length(index(V[[1]]))) {
    Y <- zoo(Y,order.by=index(V[[1]]))
  } else {
    Y <- zoo(Y,order.by=seq(nrow(V[[1]])))
  }
  Y <- attrcp(UWD,Y)
  attr(Y,'time') <- range(index(subset(X[[1]],it=it)))
  class(Y) <- class(UWD)[-1]
  if (inherits(UWD,'eof')) {
    if (verbose) print('Use dimensions and lon/lat from EOFs')
    attr(Y,'dimensions') <- c(attr(x$eof,'dimensions')[1:2],length(index(V)))
    attr(Y,'longitude') <- lon(UWD)
    attr(Y,'latidude') <- lat(UWD)
    class(Y)[1] <- 'field'
  }
  attr(Y, 'model_id') <- gcmnames
  attr(Y,'mean') <- NULL
  if(!is.null(FUN)) Y <- map(Y, FUN=FUN, plot=FALSE)
  if (verbose) {print('exit aggregate.dsensemble'); print(dim(Y))}
  return(Y)
}
