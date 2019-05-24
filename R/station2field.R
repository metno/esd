#' Transform \code{station} to \code{field} 
#'
#' Function for converting station data to field data vie the computation of PCAs
#' gridding to EOFs and then transforming the EOFs to field object.
#'
#' @seealso pca2eof eof2field
#'
#' @param x a \code{station} object
#' @param verbose a boolean; if TRUE print information about progress
#'
#' @return a \code{field} object
#'
#' @export
station2field <- function(x,verbose=FALSE) {
    if (verbose) print('station2field')
    stopifnot(inherits(x,'station'))
    x <- pcafill(x,verbose=verbose)
    pca <- PCA(x)
    eof <- pca2eof(pca,verbose=verbose)
    X <- eof2field(eof,verbose=verbose)
    return(X)
}