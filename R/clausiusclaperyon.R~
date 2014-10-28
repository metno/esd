# Clapeyron-Clausius equation (saturation evaporation pressure)

C.C.eq <- function(x) {
  stopifnot(!missing(x),(varid(x)=='t2m') | (varid(x)=='tas' |
                        (varid(x)=='sst')))
  unit <- attr(x,'unit')[1]

  # Check units:
  for (i in 1:length(unit)) {
    if ( (is.na(unit[i]) | is.null(unit[i])) ) unit[i] <- " "
    if ((unit[i]=='degree Celsius') | (unit[i]=='degC'))
         unit[i] <- 'deg C'
  }
  # If units are in deg C, convert to deg K: 
  if (unit[i]=='deg C') {
    x <- x + 273.15
    unit[i] <- 'deg K'
  }

  z <- (11.40 - 2353/x)
  e.s <- 10^z
  attr(e.s,'variable') <- 'p'
  attr(e.s,'unit') <- 'Pa'
  attr(e.s,'long_name') <- 'wapour saturation evaporation pressure'
  attr(e.s,'history') <- history.stamp(x)
  return(e.s)
}
