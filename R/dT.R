#' Calculates the t-derivatives for a time series
#' 
#' \code{dX}, \code{dY}, and \code{dT} are functions to estimate derivatives for
#' gridded field objects based on a fit to truncated Fourier series.
#' The three functions give the x-, y- and time derivatives respectively.
#' See Benestad, R.E. (2005) 'A review of the solar cycle length estimates' 
#' GRL 32 L15714, doi:10.1029/2005GL023621, August 13
#'
#' @aliases dX dY dT regfit
#' @param y A zoo, station or field object 
#' @param m number of harmonics for fitting the Fourier series 
#' @param plot if TRUE show plot 
#' @param verbose show diagnostics of the progress 
#' 
#' @return a list with several components:
#' 
#' \item{Z}{original data} \item{a}{Fourier coefficients for cosine}
#' \item{b}{Fourier coeffieicnes for sine} \item{z0}{defunct?} \item{dZ}{The
#' component contains the first derivative.} \item{dZ2}{The component contains
#' the second derivative (quicker to do both in one go).} \item{lon}{longitude}
#' \item{lat}{latitude} \item{dx}{spatial resolution} \item{span}{spatial
#' extent}
#' 
#' @examples
#' Rz <- Sunspots()
#' ## Estimate the time derivative using only the gravest Fourier frequencies
#' d.Rz.dt <- dT(Rz,m=30)
#' ## The original series
#' plot(d.Rz.dt$y)
#' ## The fitted series
#' lines(d.Rz.dt$y.fit,col='red',lty=2)
#' ## first derivative
#' plot(d.Rz.dt$dy)
#' ## Second derivative
#' plot(d.Rz.dt$dy2)
#' co2 <- CO2()
#' dco2dt <- dT(co2)
#' plot(dco2dt$y)
#' lines(dco2dt$y.fit,col='red',lty=2)
#' @export dT
dT <- function(y,m=NULL,plot=FALSE,verbose=FALSE) {

  if (verbose) print('dT')
  Y <- y
  if (plot) plot(y)
  nt <- length(y)
  if (is.null(m)) m <- nt
  m <- min(nt,m)

  ## Generate frequencies for the siusoids: rows of matrix Wi[1..m,n..nt]
  ## contain the harmonics 1..m using the outer product '%o%'
  Wi <- 2*pi/nt*(1:m)
  Wii <- Wi %*% t(seq(1,nt))
  rownames(Wii) <- paste('harmonic',1:m)

  ## Generate the terms to include in the regression:
  terms <- paste('c',(1:m), ' + s',(1:m),sep='',collapse=' + ')

  ## generate data.frame containing the original data and the harmonics
  cal.dat <- eval(parse(text=paste('data.frame(',
                          paste('c',(1:m),'= cos(Wii[',(1:m),',]),',
                                's',(1:m),'= sin(Wii[',(1:m),',])',
                                sep='',collapse=', '),')')))

  ## Apply the regression fit:
  #good <- is.finite(y)
  #if (sum(good)>10) {
  beta <- regfit(as.numeric(y),cal.dat=cal.dat,terms=terms)
  coefs <- beta[,1]
  y0 <- coefs[1]
  a <- coefs[seq(2,2*m+1,by=2)]
  b <- coefs[seq(3,2*m+1,by=2)]

  # Calculate the regression fit
  y.fit <- y0 + a%*%cos(Wii) + b%*%sin(Wii)
  dim(y.fit) <- length(y)
  # Calculate the first derivative
  dy <- (-Wi*a)%*%sin(Wii) + (Wi*b)%*%cos(Wii)
  dim(dy) <- length(y)
  # Calculate the second derivative
  dy2 <- -(Wi^2*a)%*%cos(Wii) - (Wi^2*b)%*%sin(Wii)
  dim(dy2) <- length(y)

  # Convert fit to zoo object
  if(inherits(y,'zoo')) {
    y.fit <- zoo(y.fit,order.by=index(y))
    dy <- zoo(dy,order.by=index(y))
    dy2 <- zoo(dy2,order.by=index(y))
  }
  
  if (plot) lines(y.fit,lwd=2,col="grey")
  results <- list(y=Y,a=a,b=b,y0=y0,dy=dy,dy2=dy2,y.fit=y.fit)
  class(results) <- "dydt"
  attr(results,"long_name") <- "t-derivative"
  attr(results,"descr") <- "dT.R"
  invisible(results)
}
