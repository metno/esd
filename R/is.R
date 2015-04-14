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

is.field <- function(x) (sum(inherits(x,'field')) > 0)
is.station <- function(x) (sum(inherits(x,'station')) > 0)
is.eof <- function(x) (sum(inherits(x,'EOF')) > 0)
is.pca <- function(x) (sum(inherits(x,'pca')) > 0)
is.cca <- function(x) (sum(inherits(x,'cca')) > 0)
is.storm <- function(x) (sum(inherits(x,'storm')) > 0)
is.daily <- function(x) (sum(inherits(x,'day')) > 0)
is.monthly <- function(x) (sum(inherits(x,'month')) > 0)
is.seasonal <- function(x) (sum(inherits(x,'season')) > 0)
is.annual <- function(x) (sum(inherits(x,'annual')) > 0)

is.months <- function(x) all(sum(is.element(tolower(substr(x,1,3)),
                                            tolower(month.abb)))>0)
is.seasons <- function(x) all(sum(is.element(tolower(substr(x,1,3)),
                                             names(season.abb())))>0)
is.dates <- function(x) all(!is.months(x) &
                            (levels(factor(nchar(x)))==10) |
                            (is.numeric(x) & levels(factor(nchar(x)))==8))
is.years <- function(x) all(!is.months(x) & 
                            is.numeric(x) & levels(factor(nchar(x)))==4)
