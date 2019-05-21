#' Various formulas, equations and transforms.
#' 
#' \code{C.C.eq}: Clapeyron-Clausius equation (saturation evaporation pressure)
#' where \code{x} is a data object holding the temperature.
#' 
#' \code{precip.vul}: and index for the vulerability to precipitation defined
#' as wetmean(x)/wetfreq(x). High when the mean intensity is high and/or the
#' frequency is low (it rains seldom, but when it rains, it really pours down).
#' 
#' \code{t2m.vul}: and index for the vulerability to temperature defined as the
#' mean spell length for heat waves with temperatures exceeding 30C (default).
#' 
#' \code{precip.rv}: a rough estimate of the return value for precipitation
#' under the assumption that it is exponentially distributed. Gives apprximate
#' answers for low return levels (less than 20 years). Advantage, can be
#' predicted given wet-day mean and frequency.
#' 
#' \code{nv}: number of valid data points.
#' 
#' \code{precip.Pr}: rough estimate of the probability of more than x0 of rain
#' based on an exponential distribution.
#' 
#' \code{t2m.Pr}: rough estimate of the probability of more than x0 of rain
#' based on a normal distribution.
#' 
#' \code{NE}: predicts the number of events given the probability Pr.
#' 
#' @aliases C.C.eq precip.vul t2m.vul precip.rv nv precip.Pr t2m.Pr NE
#'
#' @param x a data object
#' @param p a probability
#' @param x0 a threshold value
#' @param tau time scale (years)
#' @param is which of the spell results [1,2]
#' @param na.rm See \code{\link{mean}}.
#' @return The right hand side of the equation
#' @author R. Benestad
#'
#' @keywords parameter,element,Clausius-Clapeyron
#'
#' @examples
#' 
#' t2m <- t2m.DNMI(lon=c(-70,-10),lat=c(20,60))
#' es <- C.C.eq(t2m)
#' map(es)
#'
#' @references 'An Introduction to Atmospheric Physics' by RG. Fleagle, JA. Businger Academic Press, 9. jan. 1981, eq. 2.89, p. 72.
#'
#' @export
C.C.eq <- function(x) {
  stopifnot(!missing(x),(varid(x)=='t2m') | (varid(x)=='tas') |
                        (varid(x)=='air') | (varid(x)=='sst'))
  units <- unit(x)[1]
  if (is.null(units)) units <- 'deg C'
  
  # Check units:
  for (i in 1:length(units)) {
    if ( (is.na(units[i]) | is.null(units[i])) ) units[i] <- " "
    if ((units[i]=='degree Celsius') | (units[i]=='degC'))
         units[i] <- 'deg C'
  }
  # If units are in deg C, convert to deg K: 
  if ((units[i]=='deg C') | (units[i]=='degC')) {
    x <- x + 273.15
    units[i] <- 'deg K'
  }

  z <- (11.40 - 2353/x)
  e.s <- 10^z
  attr(e.s,'variable') <- 'p'
  attr(e.s,'unit') <- 'Pa'
  attr(e.s,'long_name') <- 'wapour saturation evaporation pressure'
  attr(e.s,'history') <- history.stamp(x)
  return(e.s)
}

#' @export
precip.vul <- function(x) {
  pv <- round(wetmean(x)/wetfreq(x))
  pv
}

#' @export
t2m.vul <- function(x,x0=30,is=1) {
  tv <- mean(subset(spell(x,threshold=x0),is=is))
  tv
}

#' @export
precip.rv <- function(x,tau=10) {
   rv <- -log( 1/(tau*365.25*wetfreq(x))) * wetmean(x)
   rv
}

#' @export
nv <- function(x) sum(is.finite(x))

#' @export
precip.Pr <- function(x,x0=10) {
  # Pr(X > x)
  mu <- wetmean(x)
  fw <- wetfreq(x)
  Pr <- fw*exp(-x0/mu)
  attr(Pr,'variable') <- paste('Pr(X>',x0,'mm/day)')
  attr(Pr,'unit') <- 'fraction'
  attr(Pr,'x0') <- x0
  attr(Pr,'size') <- nv(x)
  Pr
}

#' @export
t2m.Pr <- function(x,x0=10,na.rm=TRUE) {
    Pr <- 1-pnorm(x0,mean=mean(x,na.rm=na.rm),
                         sd=sd(x,na.rm=na.rm))
  # fix some house-keeping attributes:
   attr(Pr,'variable') <- paste('Pr(X>',x0,'deg C)')
   attr(Pr,'unit') <- 'fraction'
   attr(Pr,'x0') <- x0
   attr(Pr,'size') <- nv(x)
   Pr 
}

#' @export
NE <- function(p) {
  nel <- qbinom(p=0.05,size=attr(p,'size'),prob=p)
  nem <- qbinom(p=0.5,size=attr(p,'size'),prob=p)
  neh <- qbinom(p=0.95,size=attr(p,'size'),prob=p)
  ne <- c(nel,nem,neh)
  ## KMP 2018-11-12: x0 was not defined in the function.
  ## Should it be provided as input or as an attribute of p?
  x0 <- attr(p,"x0")
  attr(ne,'variable') <- paste('N(X>',x0,'mm/day)')
  attr(ne,'unit') <- 'days'
  return(ne)
}

