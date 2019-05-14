#' Coerce input to a \code{station} object
#' 
#' Transform an input object into the esd class \code{pca}. 
#' \code{as.pca} is an S3 method and will redirect to a fitting function depending on the output. 
#' The way in which the transformation is performed depends on the type of input data.
#' 
#' @seealso as.pca.ds as.pca.station PCA
#' 
#' @param x the input object
#' @param ... other arguments
#' 
#' @return a \code{pca} object
#' 
#' @export
as.pca <- function(x,verbose=FALSE,...) UseMethod("as.pca")

#' Coerce input to a \code{pca} object
#' 
#' Coerce a \code{ds} \code{pca} object into a \code{pca} object by simply replacing the class.
#'
#' @seealso as.pca DS
#' 
#' @param x the input object
#' @param verbose a boolean; if TRUE print information about progress
#' @param ... other arguments
#' 
#' @return a \code{pc} object
#' 
#' @export
as.pca.ds <- function(x,verbose=FALSE,...) {
  if(verbose) print("as.pca.ds")
  stopifnot(inherits(x,'pca'))
  class(x) <- class(x)[-1]
  return(x)
}

#' Coerce input to a \code{pca} object
#' 
#' Transform a \code{station} object into a \code{pca} object using the function \code{PCA}.
#'
#' @seealso as.pca PCA
#' 
#' @param x the input object
#' @param ... other arguments
#' 
#' @return a \code{pc} object
#' 
#' @export
as.pca.station <- function(x,verbose=FALSE,...) {
  if(verbose) print("as.pca.station")
  y <- PCA(x,verbose=verbose,...)
  return(y)
}
