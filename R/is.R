is.T <- function(x) {
  return(sum(is.element(tolower(attr(x,'variable')),
                        c('t2m','tmax','tmin','tas','tasmax','tasmin','air')))>0)
}

is.precip <- function(x) {
  return( (sum(is.element(tolower(varid(x)),
                        c('pr','precip','rain','precipitation',
                          'mu','fw','f[w]','tp')))>0) |
          (sum(is.element(tolower(unit(x)),
                          c('mm/day','mm/month','mm/year','mm/season','mm')))))
}
