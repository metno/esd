#' Various formulas, equations and transforms.
#' 
#' \code{C.C.eq}: Clapeyron-Clausius equation (saturation evaporation pressure)
#' where \code{x} is a data object holding the temperature.
#' 
#' @seealso precip.vul t2m.vul precip.rv precip.Pr t2m.Pr NE
#'
#' @importFrom stats qbinom
#'
#' @param x an object containing air temperature data
#' @return the saturation vapor pressure
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
  attr(e.s,'long_name') <- 'saturation vapor pressure'
  attr(e.s,'history') <- history.stamp(x)
  return(e.s)
}

#' Various formulas, equations and transforms.
#' 
#' \code{precip.vul}: an index for the vulerability to precipitation defined
#' as wetmean(x)/wetfreq(x). High when the mean intensity is high and/or the
#' frequency is low (it rains seldom, but when it rains, it really pours down).
#' 
#' @seealso C.C.eq t2m.vul precip.rv precip.Pr t2m.Pr NE
#'
#' @importFrom stats qbinom
#'
#' @param x a data object
#' @return an index for the vulerability to precipitation
#' @author R. Benestad
#'
#' @keywords parameter,element,Clausius-Clapeyron
#'
#' @export
precip.vul <- function(x) {
  pv <- round(wetmean(x)/wetfreq(x))
  return(pv)
}

#' Various formulas, equations and transforms.
#' 
#' \code{t2m.vul}: an index for the vulnerability to temperature defined as the
#' mean spell length for heat waves with temperatures exceeding 30C (default).
#' 
#' @seealso t2m.vul precip.rv precip.Pr t2m.Pr NE
#'
#' @importFrom stats qbinom
#'
#' @param x a data object
#' @param x0 a threshold value
#' @param is which of the spell results [1,2]
#' @return an index for the vulnerability to temperature
#' @author R. Benestad
#'
#' @keywords parameter,element,Clausius-Clapeyron
#'
#' @export
t2m.vul <- function(x,x0=30,is=1) {
  tv <- mean(subset(spell(x,threshold=x0),is=is))
  return(tv)
}

#' Various formulas, equations and transforms.
#' 
#' \code{precip.rv}: a rough estimate of the return value for precipitation
#' under the assumption that it is exponentially distributed. Gives apprximate
#' answers for low return levels (less than 20 years). Advantage, can be
#' predicted given wet-day mean and frequency.
#'
#' @seealso C.C.eq precip.vul t2m.vul precip.Pr t2m.Pr NE
#'
#' @importFrom stats qbinom
#'
#' @param x a data object
#' @param tau time scale (years)
#' @return An estimate of the return value for the precipitation
#' @author R. Benestad
#'
#' @keywords parameter,element,Clausius-Clapeyron
#'
#' @export
precip.rv <- function(x,tau=10) {
   rv <- -log( 1/(tau*365.25*wetfreq(x))) * wetmean(x)
   return(rv)
}

#' Various formulas, equations and transforms.
#' 
#' \code{precip.Pr}: rough estimate of the probability of more than x0 of rain
#' based on an exponential distribution.
#' 
#' @seealso C.C.eq precip.vul t2m.vul precip.Pr t2m.Pr NE
#'
#' @param x a data object
#' @param x0 a threshold value
#' @return A probability
#' @author R. Benestad
#'
#' @keywords parameter,element,Clausius-Clapeyron
#'
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
  return(Pr)
}

#' Various formulas, equations and transforms.
#' 
#' \code{t2m.Pr}: rough estimate of the probability of more than x0 of rain
#' based on a normal distribution.
#' 
#' @seealso C.C.eq precip.vul t2m.vul precip.Pr t2m.Pr NE
#'
#' @importFrom stats pnorm
#'
#' @param x a data object
#' @param x0 a threshold value
#' @param na.rm See \code{\link{mean}}.
#' @return A probability
#' @author R. Benestad
#'
#' @keywords parameter,element,Clausius-Clapeyron
#'
#' @export
t2m.Pr <- function(x,x0=10,na.rm=TRUE) {
  Pr <- 1-pnorm(x0,mean=mean(x,na.rm=na.rm),
                sd=sd(x,na.rm=na.rm))
  attr(Pr,'variable') <- paste('Pr(X>',x0,'deg C)')
  attr(Pr,'unit') <- 'fraction'
  attr(Pr,'x0') <- x0
  attr(Pr,'size') <- nv(x)
  return(Pr) 
}

#' Various formulas, equations and transforms.
#' 
#' \code{NE}: predicts the number of events given the probability Pr.
#' 
#' @seealso C.C.eq precip.vul t2m.vul precip.Pr t2m.Pr
#'
#' @param p a probability?
#' @return The right hand side of the equation
#' @author R. Benestad
#'
#' @keywords parameter,element,Clausius-Clapeyron
#'
#' @export
NE <- function(p) {
  nel <- qbinom(p=0.05,size=attr(p,'size'),prob=p)
  nem <- qbinom(p=0.5,size=attr(p,'size'),prob=p)
  neh <- qbinom(p=0.95,size=attr(p,'size'),prob=p)
  ne <- c(nel,nem,neh)
  x0 <- attr(p,"x0")
  attr(ne,'variable') <- paste('N(X>',x0,'mm/day)')
  attr(ne,'unit') <- 'days'
  return(ne)
}

