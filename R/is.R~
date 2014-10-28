is.T <- function(x) {
  return(sum(is.element(tolower(attr(x,'variable')),
                        c('t2m','tmax','tmin','tas','tasmax','tasmin','air')))>0)
}

is.precip <- function(x) {
  return(sum(is.element(tolower(attr(x,'variable')),c('pr','precip','rain','precipitation')))>0)
}
