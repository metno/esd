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

