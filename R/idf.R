#' Calculate approximate intensity-duration-frequency (IDF) curves from daily data through
#' a crude downscaling in timescales. 
#' Calculate the residual of a 'ds' object, i.e., the original data minus the fitted values.
#'
#' @aliases IDF plot.IDF
#
#' @param x a daily station object for rainfall
#' @param L timescale in hours
#' @param tau return interval in years
#' @param exponent for the power-law dependency between return values and timescales
#' @param n0 number of observations per year
#' @param alpha best fit constant and slope for scaling factor between exponential distribution and empirical 
#' distribution taken from https://doi.org/10.1088/1748-9326/ab2bb2 
#' @param verbose a boolean; if TRUE print information about progress
#' @param cols provies the colour pallette: default is ` rev(heat.colors)`.
#'
#' @examples
#' y <- station(stid=18700,src='metnod.thredds',param='precip')
#' z <- IDF(y)
#' View(z)
#'
#' @export day2IDF
day2IDF <- function(x,L=c(1,2,3,6,12,24),tau=10,zeta=NULL,n0=365.25,
                    alpha=c(1.256,0.064),verbose=FALSE) {
  ## Formula to estimate the sub-daily return-values for return interval tau and timescales L (unit: hours)
  ## The parameter zeta was calibrated on Oslo-blindern May-November while alpha was
  ## derived from many stations in Benestad et al, 2019. DOI: 10.1088/1748-9326/ab2bb2.
  ## Return-value from Benestad, Rasmus E.; Parding, Kajsa M.; Erlandsen, Helene B.; Mezghani, Abdelkader, 
  ## “A simple equation to study changes in rainfall statistics", 
  ## Environ. Res. Lett. https://doi.org/10.1088/1748-9326/ab2bb2
  ## xzeta is the linear slope in the log-log relationships for IDF curves (Lutz et al, 2020).
  ## Here a set of estimates for zeta was estimated further down in the script: log(x) = zeta * log(L/24) + const. 
  ## Choose L/24 which equals 1 for 24-hr data so that the results are tied to the 24-hr data and that
  ## x = const = return value for 24-hr data provided by x_t = alpha mu ln(fw tau).
  ## Th e value for zeta depends on the return interval tau:
  
  stopifnot((is.precip(x)) & (class(x)[2]=='day'))
  mu <- wetmean(coredata(x))
  fw <- wetfreq(coredata(x))
  taus <- c(2, 5, 10, 20, 25, 50, 100, 200) 
  zetaestimates <- c(0.4251593, 0.4185929, 0.4161947, 0.4147515, 0.4144257, 0.4137387, 0.4134449, 0.4134594)
  if (is.null(zeta)) zeta <- approx(x=taus,y=zetaestimates,tau)$y
  alpha <- alpha[1] + alpha[2]*log(tau)
  scaleL <- (L/24)^zeta
  x.L <- alpha*mu*scaleL*log(fw*tau*n0)  ## use the 24-hr estimate as a starting point.
  if (verbose) print(round(c(zeta,mu,fw,tau,n0),3))
  attr(x.L,'mu') <- mu
  attr(x.L,'fw') <- fw
  attr(x.L,'tau') <- tau
  attr(x.L,'zeta') <- zeta
  attr(x.L,'alpha') <- alpha
  attr(x.L,'n0') <- n0
  attr(x.L,'history') <- history.stamp()
  attr(x.L,'unit') <- 'mm'
  class(x.L) <- c('return_value','numeric')
  return(x.L)
}

#' @export IDF
IDF <- function(x,plot=TRUE,L=c(0.25,0.5,1,2,3,6,12,24),tau=c(2,5,10,20,50,100),cols=NULL) {
  n <- length(L); m <- length(tau)
  X <- matrix(rep(NA,n*m),n,m)
  colnames(X) <- paste(tau,'years'); rownames(X) <- paste(L,'hours')
  for (i in 1:m) X[,i] <- day2IDF(x,L=L,tau=tau[i])
  attr(X,'L') <- L
  attr(X,'tau') <- tau
  attr(X,'original_data') <- x
  class(X) <- c('IDF','matrix')
  if (plot) plot(X,cols=cols)
  return(X)
}

#' @export plot.IDF
plot.IDF <- function(x,type='l',xlab='timescale (hrs)',ylab='return value (mm)',main=NULL,cols=NULL,...) {
  d <- dim(x)
  if (is.null(cols)) cols <- rev(heat.colors(d[2]))
  if (is.null(main)) main <- paste(loc(attr(x,'original_data')),stid(attr(x,'original_data')))
  plot(attr(x,'L'),x[,d[2]],type=type,xlab=xlab,ylab=ylab,main=main,...)
  grid()
  for (i in 1:d[2]) lines(attr(x,'L'),x[,i],col="grey30",lwd=2)
  for (i in 1:d[2]) lines(attr(x,'L'),x[,i],col=cols[i],lty=2,lwd=2)
  legend(0,max(x),colnames(x),cols,bty='n')
}
