## Provenance - show the history of an object for tracability.



#' Provenance %% ~~function to do ... ~~
#' 
#' A function that extracts the history/provenance of an object. Eg. the list
#' of steps in a chain of analysis. %% ~~ A concise (1-5 lines) description of
#' what the function does. ~~
#' 
#' 
#' @param x any esd-object %% ~~Describe \code{x} here~~
#' @param what What to return from the history attribute %% ~~Describe
#' \code{what} here~~
#' @seealso history.stamp
#' @examples
#' 
#' data(Oslo)
#' provenance(Oslo)
#' 
#' @export provenance
provenance <- function(x,what='call') return(unlist(attr(x,'history')[[what]]))
