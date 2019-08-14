#' Coerce input to a \code{pca} object
#' 
#' Transform an input object into the esd class \code{pca} which is the output of principle component analysis (\code{PCA}). 
#' \code{PCA} decomposes a group of time series (a \code{station} object) into a set of spatial patterns 
#' (stored as the attribute \code{'pattern'} of the \code{pca} object), 
#' corresponding time series (the core of the \code{pca} object often referred to as principle components),
#' and eigenvalues that represent the relative strength of each spatial pattern.

#' \code{as.pca} is an S3 method and will redirect to a fitting function depending on the output. 
#' The way in which the transformation is performed depends on the type of input data.
#' 
#' @seealso PCA as.pca DS
#' 
#' @param x the input object
#' @param ... other arguments
#' 
#' @return a \code{pca} object
#' 
#' @export as.pca
as.pca <- function(x,verbose=FALSE,...) UseMethod("as.pca")

#' Coerce input to a \code{pca} object
#' 
#' Coerce a \code{ds} \code{pca} object into a \code{pca} object by replacing the class.
#'
#' @seealso 
#' 
#' @param x the input object
#' @param verbose a boolean; if TRUE print information about progress
#' @param ... other arguments
#' 
#' @return a \code{pc} object
#' 
#' @export as.pca.ds
as.pca.ds <- function(x,verbose=FALSE,...) {
  if(verbose) print("as.pca.ds")
  stopifnot(inherits(x,'pca'))
  class(x) <- class(x)[-1]
  return(x)
}

#' Coerce input to a \code{pca} object
#' 
#' Transform a \code{station} object into a \code{pca} object using the function \code{PCA}.
#' \code{PCA} decomposes a group of time series (a \code{station} object) into a set of spatial patterns 
#' (stored as the attribute \code{'pattern'} of the \code{pca} object), 
#' corresponding time series (the core of the \code{pca} object often referred to as principle components),
#' and eigenvalues that represent the relative strength of each spatial pattern.
#'
#' @seealso PCA
#' 
#' @seealso as.pca plot.pca map.pca as.station.pca DS.pca
#' 
#' @param x the input object
#' @param ... other arguments
#' 
#' @return a \code{pc} object
#' 
#' @export as.pca.station
as.pca.station <- function(x,verbose=FALSE,...) {
  if(verbose) print("as.pca.station")
  y <- PCA(x,verbose=verbose,...)
  return(y)
}
