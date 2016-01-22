## Author         A. Mezghani
## Create         11-12-2013
## Description    internal functions 
## Last update    08.01.2014


lon <- longitude <- function(x) return(attr(x,"longitude"))

lat <- latitude <- function(x) return(attr(x,"latitude"))

stid <- station_id <- function(x) return(attr(x,"station_id")) 

wmo <- function(x) return(attr(x,"wmo"))

qual <- quality <-  function(x) return(attr(x,"quality"))

alt <- altitude <- function(x) return(attr(x,"altitude"))

calendar <- function(x) return(attr(x,"calendar"))

cntr <- country <- function(x) return(attr(x,"country"))

loc <- location <- function(x) return(attr(x,"location"))

varid <- variable <- function(x) return(attr(x,"variable"))

unit <- function(x) return(attr(x, "unit"))

src <-  function(x) return(attr(x,"source"))

aspect <- function(x) return(attr(x, "aspect"))

ref <- reference <-  function(x) (attr(x, "reference"))

info <- function(x) return(attr(x,"info"))

history.esd <- function(x) return(attr(x, "history"))

pattern <- function(x) return(attr(x,"pattern"))

ele <- element <- function(x) return(attr(x,"element"))

err <- function(x) return(attr(x,"standard.error"))

ylab <- function(x) {
  unit <- unit(x)[1]
  varnm <- varid(x)[1]
  if ( (is.na(unit) | is.null(unit)) ) unit <- " "
    if ( (substr(unit,1,8)=='degree C') | (substr(unit,1,9)=='degrees C') |
         (unit=='deg C') | (unit=='degC') )
         unit <- 'degree*C'
  if (unit=='%') unit <- "'%'"
  if (varnm=='t2m') varnm <- 'T[2*m]'
  line <- paste("y <- expression(",varnm,"*phantom(0)*(",unit,"))")
  z <- try(eval(parse(text=line)),silent=TRUE)
  if (inherits(z,"try-error")) y <- unit(x)[1]
  return(y)
}
