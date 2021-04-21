#' @exportS3Method
#' @export
subset.mvcomb <- function(x, ..., it=NULL, is=NULL, ip=NULL, verbose=FALSE) {
  if(verbose) print("subset.mvcomb")  
  stopifnot(inherits(x,"mvcomb"))
  if(verbose) print(paste("ncomb=",length(lon(x))))
  if ((is.null(it)) & (is.null(is)) & is.null(ip)) return(x)
  ncomb <- length(lon(x))
  if(inherits(x, "eof")) {
    browser()
  } else {
    attrnames <- names(attributes(x))
    attrnames <- attrnames[!attrnames %in% c("dimnames")]
    y <- list()
    for(i in 1:ncomb) {
      id.i <- unique(attr(x, "id"))[i]
      x.i <- x[, attr(x, "id")==id.i]
      attr(x.i, "dimnames") <- list(index(x), paste0("x.",seq(ncol(x.i))))
      for(nm in attrnames) {
        if(is.list(attr(x, nm)) & length(attr(x, nm))==ncomb) {
          attr(x.i, nm) <- attr(x, nm)[[i]]
        }
      }
      class(x.i) <- class(x)[!class(x)=="mvcomb"]
      y[[i]] <- subset(x.i, it=it, is=is)
    }
    z <- mvcombine(y[[1]], y[[2]])
    if(ncomb>2) for(i in seq(3,ncomb)) z <- mvcombine(z, y[[i]])
  }
  attr(z,'history') <- history.stamp(x)
  invisible(z)
}
