is.T <- function(x) {
  return(sum(is.element(tolower(attr(x,'variable')),
                        c('t2m','tmax','tmin','tas','tasmax','tasmin','air','sst')))>0)
}

is.precip <- function(x) {
  return( (sum(is.element(tolower(varid(x)),
                        c('pr','precip','rain','precipitation',
                          'mu','fw','f[w]','tp')))>0) |
          (sum(is.element(tolower(unit(x)),
                          c('mm/day','mm/month','mm/year','mm/season','mm')))))
}

is.wind <- function(x) {
  return( (sum(is.element(tolower(varid(x)),
                        c('u','v','ug','vg','ff')))>0) |
          (sum(is.element(tolower(unit(x)),
                          c('m/s','km/k','knots')))))
}

is.direction <- function(x) {
  return( (sum(is.element(tolower(varid(x)),c('dd','direction','heading')))>0) |
          (sum(is.element(tolower(unit(x)),c('degrees','radians')))))
}

is.pressure <- function(x) {
  return( (sum(is.element(tolower(varid(x)),
                        c('p','slp','psl','pressure')))>0) |
          (sum(is.element(tolower(unit(x)),
                          c('hPa','Pa','millibar','bar','N/m^2')))))
}


is.field <- function(x) (sum(inherits(x,'field')) > 0)
is.ds <- function(x) (sum(inherits(x,'ds')) > 0)
is.dsensemble <- function(x) (sum(inherits(x,'dsensemble')) > 0)
is.station <- function(x) (sum(inherits(x,'station')) > 0)
is.eof <- function(x) (sum(inherits(x,'eof')) > 0)
is.pca <- function(x) (sum(inherits(x,'pca')) > 0)
is.cca <- function(x) (sum(inherits(x,'cca')) > 0)
is.trajectory <- function(x) (sum(inherits(x,'trajectory')) > 0)
is.daily <- function(x) (sum(inherits(x,'day')) > 0) | inherits(x,'Date')
is.monthly <- function(x) (sum(inherits(x,'month')) > 0)
is.seasonal <- function(x) (sum(inherits(x,'season')) > 0)
is.annual <- function(x) (sum(inherits(x,'annual')) > 0)

is.months <- function(x) all(sum(is.element(tolower(substr(x,1,3)),
                                            tolower(month.abb)))>0)
is.seasons <- function(x) all(sum(is.element(tolower(substr(x,1,3)),
                                             names(season.abb())))>0)
## KMP 2017-06-07 Why wasn't inherits(x,"Date") included in is.dates? Changed to solve problem in subset.pc
is.dates <- function(x) all(inherits(x,"Date") | !is.months(x) &
                            (levels(factor(nchar(x)))==10) |
                            (is.numeric(x) & levels(factor(nchar(x)))==8))
is.years <- function(x) all(!is.months(x) & 
                            is.numeric(x) & levels(factor(nchar(x)))==4)

is.model <- function(model,verbose=FALSE) {
  if (verbose) print(summary(model))
  return(inherits(model,c('lm','glm','mlm')))
}
