#' Test object type
#' 
#' Test if an object is of a certain class or contains some variable
#' 
#' @aliases is.T is.precip is.field is.station is.eof is.pca is.cca
#' is.trajectory is.daily is.monthly is.seasonal is.annual is.model is.wind
#' is.direction is.pressure is.dates is.ds is.dsensemble is.events is.months
#' is.seasons is.url is.years
#'
#' @param x a data object
#' @return a boolean
#'
#' @keywords parameter,element
#'
#' @examples
#' data(ferder)
#' is.T(ferder)
#' is.precip(ferder)
#'
#' @export
is.T <- function(x) {
  return(sum(is.element(tolower(attr(x,'variable')),
                        c('t2m','tmax','tmin','tas','tasmax','tasmin','air','sst')))>0)
}

#' @export
is.precip <- function(x) {
  return( (sum(is.element(tolower(varid(x)),
                        c('pr','precip','rain','precipitation',
                          'mu','fw','f[w]','tp')))>0) |
          (sum(is.element(tolower(unit(x)),
                          c('mm/day','mm/month','mm/year','mm/season','mm')))))
}

#' @export
is.wind <- function(x) {
  return( (sum(is.element(tolower(varid(x)),
                        c('u','v','ug','vg','ff')))>0) |
          (sum(is.element(tolower(unit(x)),
                          c('m/s','km/k','knots')))))
}

#' @export
is.direction <- function(x) {
  return( (sum(is.element(tolower(varid(x)),c('dd','direction','heading')))>0) |
          (sum(is.element(tolower(unit(x)),c('degrees','radians')))))
}

#' @export
is.pressure <- function(x) {
  return( (sum(is.element(tolower(varid(x)),
                        c('p','slp','psl','pressure')))>0) |
          (sum(is.element(tolower(unit(x)),
                          c('hPa','Pa','millibar','bar','N/m^2')))))
}

#' @export
is.field <- function(x) (sum(inherits(x,'field')) > 0)

#' @export
is.ds <- function(x) (sum(inherits(x,'ds')) > 0)

#' @export
is.dsensemble <- function(x) (sum(inherits(x,'dsensemble')) > 0)

#' @export
is.station <- function(x) (sum(inherits(x,'station')) > 0)

#' @export
is.eof <- function(x) (sum(inherits(x,'eof')) > 0)

#' @export
is.pca <- function(x) (sum(inherits(x,'pca')) > 0)

#' @export
is.cca <- function(x) (sum(inherits(x,'cca')) > 0)

#' @export
is.trajectory <- function(x) (sum(inherits(x,'trajectory')) > 0)

#' @export
is.events <- function(x) (sum(inherits(x,'events')) > 0)

#' @export
is.daily <- function(x) (sum(inherits(x,'day')) > 0) | inherits(x,'Date')

#' @export
is.monthly <- function(x) (sum(inherits(x,'month')) > 0)

#' @export
is.seasonal <- function(x) (sum(inherits(x,'season')) > 0)

#' @export
is.annual <- function(x) (sum(inherits(x,'annual')) > 0)

#' @export
is.months <- function(x) all(sum(is.element(tolower(substr(x,1,3)),
                                            tolower(month.abb)))>0)

#' @export
is.seasons <- function(x) all(sum(is.element(tolower(substr(x,1,3)),
                                             names(season.abb())))>0)

#' @export
is.dates <- function(x) all(inherits(x,"Date") | !is.months(x) &
                            (levels(factor(nchar(x)))==10) |
                            (is.numeric(x) & levels(factor(nchar(x)))==8))

#' @export
is.years <- function(x) all(!is.months(x) & 
                            is.numeric(x) & levels(factor(nchar(x)))==4)

#' @export
is.model <- function(model,verbose=FALSE) {
  if (verbose) print(summary(model))
  return(inherits(model,c('lm','glm','mlm')))
}

#' @export
is.url <-function(x) {
  grepl("http:|https:|www.", x)
}