#' Aggregate data - calculate 5 year average
#' 
#' @export
pentad <- function(x,l=5,it0=NULL,...) {
  if (!is.null(it0)) yr <- year(x) - it0 else yr <- year(x)
  yrl <- l*trunc(yr/l)
  if (!is.null(it0)) yrl <- yrl + it0  
  index(x) <- yrl
  
  xl <- aggregate(x,yrl,...)
  attr(xl,'dimnames') <- NULL
  return(xl)
}