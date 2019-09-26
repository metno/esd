#' Transform \code{pca} to \code{eof}
#'
#' \code{pca2eof} uses gridding (\code{\link{gridstation}}) to transform the
#' station map into a regular map making use of elevation.
#'
#' @param x a 'pca' object (see \code{link{PCA}})
#' @param verbose a boolean; if FALSE do not print information about progress (silent mode)
#'
#' @return an 'eof' object (see \code{link{EOF}})
#'
#' @seealso gridstation
#'
#' @export
pca2eof <- function(x,verbose=FALSE) {
  if (verbose) print('pca2eof')
  stopifnot(inherits(x,'pca'))
  y <- x
  if (!is.null(z <- attr(x,'pattern'))) {
    z <- attr(x,'pattern')
    attr(z,'longitude') <- lon(x)
    attr(z,'latitude') <- lat(x)
    attr(z,'altitude') <- alt(x)
  } else if (!is.null(x$pca)) z <- x$pca else
                              stop('Do not know how to handle this object!')
  d <- dim(z)
  Z <- list()
  if (verbose) print('Grid the modes')
  for (i in 1:d[2]) {
    Z[[i]] <- gridstation(z,i,verbose=verbose)
  }
  if (verbose) print('Grid the mean')
  zc <- attr(x,'mean'); dim(zc) <- c(length(zc),1)
  attr(zc,'longitude') <- lon(x)
  attr(zc,'latitude') <- lat(x)
  attr(zc,'altitude') <- alt(x)  
  clim <- gridstation(zc,1,verbose=verbose)
  z <- unlist(Z)
  dim(z) <- c(dim(Z[[1]]),d[2])
  z -> attr(y,'pattern')
  clim  -> attr(y,'mean')
  attr(y,'longitude') <- lon(Z[[1]])
  attr(y,'latitude') <- lat(Z[[1]])
  attr(y,'old_longitude') <- lon(zc)
  attr(y,'old_latitude') <- lat(zc)
  attr(y,'old_altitude') <- alt(zc)
  attr(y,'dimensions') <- c(dim(Z[[1]]),d[2])
  attr(y,'variable') <- varid(x)[1]
  attr(y,'unit') <- unit(x)[1]
  attr(y,'longname') <- attr(x,'longname')[1]
  attr(y,'greenwich') <- TRUE
  class(y) <- c('eof','field',class(x)[-c(1,2)])
  return(y)
}
