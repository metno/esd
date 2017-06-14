## Provenance - show the history of an object for tracability.

provenance <- function(x,what='call') return(unlist(attr(x,'history')[[what]]))