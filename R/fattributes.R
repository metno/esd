## Author         A. Mezghani
## Create         11-12-2013
## Description    internal functions 
## Last update    08.01.2014

#' Shortcuts to attributes
#'
#' Fast access to attributes, e.g, \code{lon(x)} gives you the same output as \code{attr(x,"longitude")}.
#'
#' @aliases lon longitude lat latitude alt altitude stid qual quality calendar
#' cntr country loc location varid variable unit src aspect ref reference info information
#' history.esd pattern ele element err
#'
#' @param x input object
#'
#' @export
lon <- longitude <- function(x) {
  if (inherits(x,'trajectory')) return(x[,is.element(colnames(x),'lon')]) else 
    return(attr(x,"longitude"))
}

#' @export
lat <- latitude <- function(x) {
  if (inherits(x,'trajectory')) return(x[,is.element(colnames(x),'lat')]) else 
    return(attr(x,"latitude"))
}

#' @export
stid <- function(x) return(attr(x,"station_id")) 

# NOT EXPORTED
wmo <- function(x) return(attr(x,"wmo"))

#' @export
qual <- quality <-  function(x) return(attr(x,"quality"))

#' @export
alt <- altitude <- function(x) return(attr(x,"altitude"))

#' @export
calendar <- function(x) return(attr(x,"calendar"))

#' @export
cntr <- country <- function(x) return(attr(x,"country"))

#' @export
loc <- location <- function(x) { 
  if (!is.null(attr(x,"location"))) return(attr(x,"location"))
  if (inherits(x,'list')) {
    if (!is.null(x$pca)) return(attr(x$pca,"location")) else
    if (!is.null(names(x))) return(names(x))
  }
}

#' @export
varid <- variable <- function(x) return(attr(x,"variable"))

#' @export
unit <- function(x) return(attr(x, "unit"))

#' @export
src <-  function(x) return(attr(x,"source"))

#' @export
aspect <- function(x) return(attr(x, "aspect"))

#' @export
ref <- reference <-  function(x) (attr(x, "reference"))

#' @export
info <- information <- function(x) return(attr(x,"info"))

#' @export
history.esd <- function(x) return(attr(x, "history"))

#' @export
pattern <- function(x) return(attr(x,"pattern"))

#' @export
ele <- element <- function(x) return(attr(x,"element"))

#' @export
err <- function(x) return(attr(x,"standard.error"))

#' A function that extracts the history/provenance of an object. Eg. the list
#' of steps in a chain of analysis.
#' 
#' @param x any esd-object 
#' @param what What to return from the history attribute 
#' 
#' @seealso history.stamp
#' 
#' @examples
#' data(Oslo)
#' provenance(Oslo)
#' 
#' @export provenance
provenance <- function(x,what='call') return(unlist(attr(x,'history')[[what]]))

#' Create a label for plots
#'
#' Using the attributes of an object, put together a character string with 
#' the variable name and unit to be used as label in plots and maps.
#'
#' @param x an input object
#'
#' @return a character string
#'
#' @export
ylab <- function(x) {
  unit <- unit(x)[1]
  varnm <- varid(x)[1]
  if ( (is.na(unit) | is.null(unit)) ) unit <- " "
  if (is.character(unit)) {
    if ( (substr(unit,1,8)=='degree C') | (substr(unit,1,9)=='degrees C') |
         (unit=='deg C') | (unit=='degC') )
           unit <- 'degree*C'
    if (unit=='%') unit <- "'%'"
  }
  if (is.character(varnm)) if (varnm=='t2m') varnm <- 'T[2*m]'
  line <- paste("y <- expression(",varnm,"*phantom(0)*(",unit,"))")
  z <- try(eval(parse(text=line)),silent=TRUE)
  if (inherits(z,"try-error")) y <- unit(x)[1]
  return(y)
}
