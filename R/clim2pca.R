#' PCA analysis of the seasonal cycle 
#'
#' Express station data as \code{\link{PCA}} where each of the EOFs
#' (attribute 'pattern' of output object)
#' represent one year and the PCs
#' (main part of ouput object) describe the seasonal variations.
#'
#' @aliases clim2pca.default clim2pca.month clim2pca.day
#'
#' @export
clim2pca <- function(x,verbose=FALSE,...) UseMethod("clim2pca")

#' @export
clim2pca.default <- function(x,verbose=FALSE,...) {
  if(verbose) print("clim2pca.default - unfinished function returning input object")
  return(x)
}

#' @export
clim2pca.month <- function(x,verbose=FALSE,...) {
  if(verbose) print("clim2pca.month")
  X <- aggregate(x,year)
  ny <- length(x) %/% 12
  nm <- length(x) %% 12
  y <- coredata(x[1:(length(x)-nm)])
  dim(y) <- c(12,ny)
  ok <- is.finite(colMeans(y))
  pca <- svd(y[,ok])
  for (i in 1:12) {
    z <- zoo(pca$v[,i],order.by=index(X))
    if (i == 1) Z <- z else
                Z <- merge(Z,z)
  }
  season <- pca$u
  colnames(season) <- month.abb
  rownames(season) <- paste("pattern",1:12,sep=".")
  attr(Z,'season') <- season
  attr(Z,'d') <- pca$d
  return(Z)
}

#' @export
clim2pca.day <- function(x,verbose=FALSE,...) {
  if(verbose) print("clim2pca.day - unfinished function returning input object")
  return(x)
}