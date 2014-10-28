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

precip.vul <- function(x) {
  pv <- round(wetmean(x)/wetfreq(x))
  pv
}

t2m.vul <- function(x,x0=30,is=1) {
  tv <- mean(subset(spell(x,threshold=x0),is=is))
  tv
}

precip.rv <- function(x,tau=10) {
   rv <- -log( 1/(tau*365.25*wetfreq(x))) * wetmean(x)
   rv
}

nv <- function(x) sum(is.finite(x))

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

NE <- function(p) {
  nel <- qbinom(p=0.05,
                size=attr(p,'size'),prob=p)
  nem <- qbinom(p=0.5,
                   size=attr(p,'size'),prob=p)
  neh <- qbinom(p=0.95,
                   size=attr(p,'size'),prob=p)
  ne <- c(nel,nem,neh)
  attr(ne,'variable') <- paste('N(X>',x0,'mm/day)')
  attr(ne,'unit') <- 'days'
  ne
}

