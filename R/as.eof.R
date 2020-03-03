#' Coerce input to an \code{eof} object
#' 
#' Transform an input object into the esd class \code{eof} (see \code{\link{EOF}}).
#' 
#' \code{as.eof} is an S3 method and will redirect to a fitting function depending on
#' the class of the input object. 
#'
#' \code{as.eof.dsensemble.pca} converts PCA-based DSensemble objects to EOF-based results (gridded)
#'
#' @aliases as.eof as.eof.field as.eof.comb as.eof.list as.eof.zoo as.eof.ds
#' as.eof.eof as.eof.dsensemble as.eof.appendix as.eof.dsensemble.pca 
#' 
#' @param x the input object
#' @param iapp index of appendix
#' @param ... other arguments
#' 
#' @return an \code{eof} object
#' 
#' @export as.eof
as.eof <- function(x,...) UseMethod("as.eof")

#' @export as.eof.zoo
as.eof.zoo <- function(x,...,verbose=FALSE) {
  if(verbose) print("as.eof.zoo")
  class(x) <- c('eof','zoo')
  return(x)
}

#' @export as.eof.ds
as.eof.ds <- function(x,...,iapp=NULL,verbose=FALSE) {
  if(verbose) print("as.eof.ds")
  y <- as.eof(attr(x,'eof'),iapp=iapp) 
  return(y)
}

#' @export as.eof.eof
as.eof.eof <-function(x,...,iapp=NULL) {
  if (inherits(x,'comb')) {
    x <- as.eof.comb(x,iapp=iapp) 
  } 
  return(x)
}

#' @export as.eof.comb
as.eof.comb <- function(x,...,iapp=NULL,verbose=FALSE) {
  if(verbose) print("as.eof.comb")
  stopifnot(inherits(x,'comb'))

  # if x is a 'field'
  if (!inherits(x,'eof')) x <- EOF(x)

  # assume x from now on is an 'eof'
  if (!is.null(iapp)) {
    y <- as.eof.appendix(x,...,iapp=iapp)
    return(y)
  }
  class(x) <- class(x)[-grep('comb',class(x))]
  napps <- attr(x,'n.apps')
  for (i in seq(napps)) {
    eval(parse(text=paste("attr(x,'appendix.",i,"') <- NULL",sep="")))
  }
  attr(x,'n.apps') <- NULL
  attr(x,'history') <- history.stamp(x)
  return(x)
}

#' @export as.eof.field
as.eof.field <- function(x,...,iapp=NULL,verbose=FALSE) {
  if(verbose) print("as.eof.field")
  y <- EOF(x,...)
  if (!is.null(iapp)) y <- as.eof.appendix(y,iapp=iapp)
  return(y)
}

#' @export as.eof.appendix
as.eof.appendix <- function(x,...,iapp=1,verbose=FALSE) {
  if (verbose) print("as.eof.appendix")
  stopifnot(inherits(x,'comb'))
  clim <- eval(parse(text=paste("attr(attr(x,'appendix.",iapp,"'),'climatology')",sep="")))
  aveg <- eval(parse(text=paste("attr(attr(x,'appendix.",iapp,"'),'mean')",sep="")))
  y <- eval(parse(text=paste("attr(x,'appendix.",iapp,"')",sep="")))
  x <- as.eof.comb(x)
  y <- attrcp(x,y)
  if (!is.null(clim)) attr(y,'climatology') <- clim 
  if (!is.null(aveg)) attr(y,'mean') <- aveg
  attr(y,'history') <- history.stamp(x)
  class(y) <- class(x)
  return(y)
}

#' @export as.eof.list
as.eof.list <- function(x,...,verbose=FALSE) {
  if (verbose) print('as.eof.list')
  stopifnot(inherits(x,'list'),inherits(x[[1]],'eof'))
  
  wPC <- function(z,iapp=NULL) {
    eigv <- attr(z,'eigenvalues')
    w <- eigv/sum(eigv)
    if (is.null(iapp)) Z <- z %*% diag(w) else
                       Z <- attr(z,paste('appendix.',iapp,sep='')) %*% diag(w)
    Z <- zoo(Z,order.by=index(z))
    return(Z)
  }

  if (verbose) try(print(summary(x)))
  if (inherits(x[[1]],'character')) x[[1]] <- NULL
  if (inherits(x[[1]],'eof')) {eof <- x[[1]]; x[[1]] <- NULL}
  X.list <- lapply(x,wPC)
  X <- do.call("merge", X.list)
  if (verbose) print(summary(X))
  t <- index(X)
  udv <- svd(coredata(X))
  eof <- zoo(udv$u[,1:20],order.by=t)
  attr(eof,'eigenvalues') <- udv$d
  pattern <- rep(1,dim(udv$v)[1])
  names(pattern) <- names(X)
  attr(eof,'pattern') <- pattern
  if (inherits(x[[1]],'comb')) {
    if (verbose) print('Combined field: appendix.1')
    for (i in 1:attr( attr(x[[1]],'n.apps'))) {
      z.list <- lapply(x,wPC,iapp=i)
      udv1 <- svd(coredata(do.call("merge", z.list)))
      attr(eof,paste('appendix.',i,sep='')) <- zoo(udv1$u[,1:20],
               order.by=index(attr(x,paste('appendix.',i,sep=''))))
      names(attr(eof,paste('appendix.',i,sep=''))) <- paste("X.",1:20,sep="")
    }
  }
  attr(eof,'original.list.of.eofs') <- x
  attr(eof,'udv') <- udv
  id <- c()
  for (i in 1:length(x)) id <- c(id,rep(i,length(attr(x[[i]],'eigenvalues'))))
  attr(eof,'id') <- id
  names(eof) <- paste("X.",1:20,sep="")
  class(eof) <- class(x[[1]])
  return(eof)
}

#' @export as.eof.dsensemble
as.eof.dsensemble <- function(x,...,FUN='mean',verbose=FALSE) {
  ## R.E. Benestad, 2017-05-19
  ## Convert the dsensemble object to an EOF of the multi-model mean
  if (verbose) print('as.eof.dsensemble')
  stopifnot(inherits(x,'dsensemble'),inherits(x[[2]],'eof')|inherits(x[[2]],'pca'))
  eof0 <- x[[2]]
  x[[2]] <- NULL
  x[[1]] -> info
  x[[1]] <- NULL
  d <- c(dim(x[[1]]),length(x))
  y <- unlist(x)
  dim(y) <- c(d[1]*d[2],d[3])
  Y <- apply(y,1,FUN)
  dim(Y) <- c(d[1],d[2])
  eof <- zoo(Y,order.by=index(x[[1]]))
  eof <- attrcp(eof0,eof)
  class(eof) <- class(eof0)
  attr(eof,'info') <- info
  attr(eof,'history') <- history.stamp()
  return(eof)
}

#' @export as.eof.dsensemble.pca
as.eof.dsensemble.pca <- function(x,...,is=NULL,it=NULL,ip=NULL,verbose=FALSE) {
  if (verbose) print('as.eof.dsensemble.pca')
  stopifnot(inherits(x,"dsensemble") & inherits(x,"pca"))
  if (inherits(x,"eof")) {
    invisible(x)
  } else {
    eof <- pca2eof(x$pca)
    eof <- subset(eof,ip=ip)
    if (!is.null(is)) eof <- subset(eof,is=is,it=it,verbose=verbose)
    x$eof <- eof 
    class(x) <- c("dsensemble", "eof", "list")
    invisible(x)
  }
}
