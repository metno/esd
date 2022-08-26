#' expandpca
#' 
#' The function \code{expandpca} is used to extract information from PCA-based dsensemble-objects,
#' but takes care of extra house keeping, such as attributes with meta data. It applies operations onto the 
#' PCs through the argument 'FUNX' to estimate quantities such as the ensemble mean. The argument 'FUN' is used to 
#' aggregate the results from e.g. the ensemble mean to provide maps of means or trends. 
#' Hence these functions are used to expand PCA to obtain station data or EOFs into field data.
#' \code{aggregate.dsensemble} is used for aggregating ensembles in a similar fashion, but is slower and applies 
#' the operation onto the data after it has been expanded by matrix multiplication of the singular vectors 
#' (X = U %*% diag(W) %*% t(V)). 
#' \code{map.dsensemble} is a wrapper for \code{map} that uses \code{expandpca} or \code{aggregate.dsensemble} to distill selected
#' information.
#' @seealso aggregate map.dsensemble aggregate.dsensemble
#' @aliases expandpca
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
    #if (!FUN %in% c('min','mean','max') anomaly <- TRUE
    if (FUN %in% c('trend','trend.coef','sd')) anomaly <- TRUE
    if (verbose) print(c(FUN,anomaly))
  }
  if (verbose) print(names(attributes(UWD)))
  ## Eigenvalues
  D <- attr(UWD,'eigenvalues')
  ## Create a matrix with only the GCM time series
  if (verbose) print('PCA/EOF-based ensemble')
  X <- x
  X$info <- NULL; X$pca <- NULL; X$eof <- NULL
  n <- length(X)
  dims <- rep('?',n)
  for (ii in 1:n) dims[ii] <- paste(dim(X[[ii]]),collapse='x')
  if (verbose) print(table(dims))
  ## Dimension of each downscaled GCM results
  if (verbose) print(paste('subset, it=',it))
  if (!is.null(it)) X <- subset(X,it=it)
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
  if (is.null(it)) it <- range(index(X[[1]]))
  memsiz <- rep("?",n)
  for (i in 1:n) memsiz[i] <- paste(dim(X[[i]]),collapse='x')
  memsiztab <- table(memsiz)
  if (verbose) print(memsiztab)
  ## Quality control
  if (verbose) print(paste('Before quality control: original number of members=',n))
  for (i in seq(n,1,by=-1)) {
    if (verbose) {print(range(X[[i]],na.rm=TRUE)); print(dim(X[[i]]))}
    if (max(abs(X[[i]]),na.rm=TRUE) > 10)  {
      print(paste(i,'Remove suspect results')); X[[i]] <- NULL
    }
  }
  n <- length(X)
  if (verbose) print(paste('After quality control: new number of members=',n))
  
  d <- dim(X[[1]])
  if (verbose) {print(names(X)); print(c(d,n,length(unlist(X)))); print(paste('FUNX=',FUNX))}
  ## Apply FUNX to each of the PCs across all members
  lengths <- rep(NA,length(X))
  for (i in 1:length(X)) lengths[i] <- length(index(X[[i]]))
  #if (length(table(lengths))>1) 
  if (!test) {
    V <- unlist(X)
    V[abs(V) > 10] <- NA
    dim(V) <- c(d[1]*d[2],n)
    if (verbose) print(FUNX)
    V <- apply(V,1,FUN=FUNX)
    dim(V) <- d
    # if (verbose) print('alternative way...')
    # VV <- V[[1]]
    # for (i in 2:length(V)) {
    #   print(range(V[[i]],na.rm=TRUE)); print(dim(V[[i]]))
    #   VV <- VV + V[[i]]
    # }
    # V <- VV/length(V)
  } else {
    V <- X[[1]] # Pick one member for testing ## Testing
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
  #if(nrow(V[[1]])==length(index(subset(X[[1]],it=it)))) {
  # if(dim(V)[1]==length(index(subset(X[[1]],it=it)))) {
  #   Y <- zoo(Y,order.by=index(subset(X[[1]],it=it)))
  # } else {
  #   Y <- zoo(Y,order.by=seq(nrow(V[[1]])))
  # }
  Y <- zoo(Y,order.by=index(X[[1]]))
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
  attr(Y,'mean') <- NULL
  if (verbose) {print(range(it)); print(dim(Y)); print('exit expandpca')}
  return(Y)
}
