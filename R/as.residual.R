#' Calculate residual
#' 
#' Caluculate the residual of a 'ds' object, i.e., the original data minus the fitted values.
#'
#' @aliases as.residual as.residual.ds as.residual.station as.residual.field
#
#' @param x a 'ds' object
#' @param verbose a boolean; if TRUE print information about progress
#' @param \dots additional arguments
#'
#' @export 
as.residual <- function(x,...) UseMethod("as.residual")

#' @export
as.residual.ds <- function(x,...,verbose=FALSE){
  if (verbose) print('as.residual.ds')
  if (is.ds(x)) {
    ## If the predictand was originally an EOF or PCA product, then
    ## the residual needs to inherits their attributes
    if (verbose) print('Re-construct and re-compute')
    ## Need to reconstruct the data matrix and re-calculate the EOFs/PCAs 
    if (is.eof(x)) {
      if (verbose) print('eof/field')
      z0 <- as.field(attr(x,'original_data'))
      z1 <- as.field(x)
      y <- z1 - z0
      y <- attrcp(z0,y); class(y) <- class(z0)
    } else if (is.pca(x)) {
      if (verbose) print('pca/station')
      z0 <- as.station(attr(x,'original_data'))
      z1 <- as.station(x)
      y <- z1 - z0
      y <- attrcp(z0,y); class(y) <- class(z0)
    } else if (is.station(x)) {
      if (verbose) print('station')
      z0 <- attr(x,'original_data')
      z1 <- x
      y <- z1 - z0
      y <- attrcp(z0,y); class(y) <- class(z0)
    }      
  }
  ## If the results are a field object, then the residuals are stored as EOFs.
  if (is.field(x)) {
    if (verbose) {print('x is a field object'); print(class(x))}
    y <- as.field(attr(x,'original_data')) -
         as.field(attr(x,'fitted_values'))
    y <- attrcp(attr(x,'original_data'),y)
    class(y) <- class(attr(x,'original_data'))
    y <- as.field(y)
    attr(y,'aspect') <- 'residual'
  }
  attr(y,'aspect') <- 'residual'
  attr(y,'history') <- history.stamp(x)
  invisible(y)
}

#' @export
as.residual.station <- function(x,...) {
  if (!is.null(attr(x,'calibration_data'))) {
    y <- as.residual.ds(x)
  } else {
    y <- NULL
  }
  invisible(y)
}

#' @export
as.residual.field <- function(x,...) {
  if (!is.null(attr(x,'calibration_data'))) {
    y <- as.residual.ds(x)
  } else {
    y <- NULL
  }
  invisible(y)
}

