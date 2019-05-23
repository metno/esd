#' Coerce input to a \code{comb} object
#' 
#' Transform an \code{eof} object into the esd class \code{comb} by applying \code{eof2field} to the input object and its appendices. 
#' 
#' @aliases as.comb as.comb.eof
#'
#' @seealso eof2field
#' 
#' @param x the input object
#' @param verbose a boolean; if TRUE print information about progress
#' @param ... other arguments
#' 
#' @return an \code{eof} object
#' 
#' @export
as.comb <- function(x,verbose=FALSE,...) UseMethod("as.comb")

#' @export
as.comb.eof <- function(x,verbose=FALSE,...) {
  if(verbose) print("as.comb.eof")
  stopifnot(inherits(x,'eof'))
  y <- eof2field(as.eof(x))
  n <- attr(x,'n.apps')
  if (!is.null(n)) {
    for ( i in 1:n ) {
      z <- as.eof(x,i)
      #plot(z)
      Z <- eof2field(z)
      eval(parse(text=paste("Z -> attr(y,'appendix.",i,"')",sep="")))
    }
  }
  n -> attr(y,'n.apps')
  attr(y,'history') <- history.stamp(x)
  attr(y,'quality') <- c(attr(x,'quality'),'EOF-filtered')
  return(y)
}
