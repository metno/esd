#' Copy attributes
#'
#' Copy the attributes of x to y
#' 
#' @param x input object with attributes
#' @param y input object to receive attributes
#' @param ignore names of attributes that will not be copied from x to y
#'
#' @seealso softattr
#'
#' @return same as y but with attributes of x
#'
#' @export
attrcp <- function(x,y,ignore=c("name","model","n.apps","appendix","dimnames")) {
  nattr <- softattr(x,ignore=ignore)
  if (length(nattr)>0) {
    for (i in 1:length(nattr)) {
      attr(y,nattr[i]) <- attr(x,nattr[i])
    }
  }
  if (inherits(x,'field')) {
    if (!is.null(dim(x)) & !is.null(dim(y)) & (length(index(y))>0)) { 
      if (dim(y)[2]==dim(x)[2]) {
        attr(y,'dimensions') <- c(attr(x,'dimensions')[1:2],length(index(y)))
      }
    }
  }
  invisible(y)
}

#' Get names of attributes
#'
#' @param x input object with attributes
#' @param ignore names to ignore
#'
#' @return character string or vector with names of attributes
#'
#' @seealso attrcp
#'
#' @export
softattr <- function (x, ignore = NULL) {
  nattr <- names(attributes(x))
  if (sum(is.element(nattr, "names")) > 0) {
    nattr <- nattr[-grep("names", nattr)]
  }
  if (sum(is.element(nattr, "index")) > 0) {
    nattr <- nattr[-grep("index", nattr)]
  }
  if (sum(grep("dim",nattr)) > 0) {
    nattr <- nattr[-grep("dim",nattr)]
  }
  if (!is.null(ignore)) {
    for (i in 1:length(ignore)) {
      if (sum(is.element(nattr, ignore[i])) > 0) {
        nattr <- nattr[-grep(ignore[i], nattr)]
      }
    }
  }
  if (sum(is.element(nattr, "class")) > 0) {
    nattr <- nattr[-grep("class", nattr)]
  }
  return(nattr)
}
