#' Coherence spectrum - cross-spectrum analysis
#' 
#' Based on: \url{http://en.wikipedia.org/wiki/Wiener-Khinchin_theorem}; Press
#' et al. (1989) 'Numerical Recipes in Pascal', Cambridge, section 12.8
#' 'Maximum Entropy (All Poles) Method'; von Storch & Zwiers (1999)
#' 'Statistical Analysis in climate Research', Cambridge, section 11.4, eq
#' 11.67, p. 235;
#' 
#' A test with two identical series the original equation (eq 11.67) from von
#' Storch & Zwiers (1999) gave uniform values: 1. The denominator was changed
#' from \eqn{( \Gamma_{xx} * \Gamma_{yy} )}{(Gamxx * Gamyy )} to \eqn{(
#' }{sqrt(Gamxx * Gamyy )}\eqn{ \sqrt{\Gamma_{xx} * \Gamma_{yy} )}}{sqrt(Gamxx
#' * Gamyy )}.
#' 
#' 
#' @param x A vector (time series).
#' @param y A vector (time series).
#' @param dt time incremet - for plotting.
#' @param M Window length - default= half series length
#' @param plot Flag: plot the diagnostics.
#'
#' @return A complex vector .
#' 
#' @keywords manip
#' 
#' @export coherence
coherence <- function(x,y,dt=1,M=NULL,plot=TRUE) {
# Based on:  
# http://en.wikipedia.org/wiki/Wiener%E2%80%93Khinchin_theorem
# Press et al. (1989) 'Numerical Recipes in Pascal', Cambridge, section 12.8
#                     'Maximum Entropy (All Poles) Method'  
# von Storch & Zwiers (1999) 'Statistical Analysis in climate Research',
#                     Cambridge, section 11.4, eq 11.67, p. 235   
# Test with two identical series the original equation gave uniform values: 1
# The denominator was changed from
#  ( Gamxx * Gamyy ) to sqrt( Gamxx * Gamyy )

xy <- merge(zoo(x),zoo(y),all=TRUE)
x <- xy$x
y <- xy$y

if (is.null(M)) M <- length(x)/2
t <- ( 1:length(x) )/dt
n <- seq(-M,M,length=length(x)+1)  # Press et al (1989): eq.12.1.5
n <- c(sort(n[n>=0]),sort(n[n<0]))
#print(n)
fn <- n/(M * dt)
Phixy <- ccf(x, y, lag.max = M, plot=FALSE,type = "covariance")
Gamxy <- fft(Phixy$acf)

Phixx <- ccf(x, x, lag.max = M, plot=FALSE,type = "covariance")
Gamxx <- fft(Phixx$acf)

Phiyy <- ccf(y, y, lag.max = M, plot=FALSE,type = "covariance")
Gamyy <- fft(Phiyy$acf)

Axy <- abs(Gamxy)

coh <- Axy^2/sqrt( Gamxx * Gamyy )

attr(coh,'Reference') <- "von Storch & Zwiers (1999),  eq 11.67, p. 235"
attr(coh,'Method') <- "coh in clim.pact"
attr(coh,'Description') <- "Squared coherence"

n <- ceiling(length(x)/2)
tau <- 1/fn

plus <- tau>0
if (plot) {
  ylim <- c(min(log(abs(coh)),na.rm=TRUE),max(log(abs(coh)),na.rm=TRUE)*1.1)
  plot(tau[plus],log(abs(coh[plus])),type="l",lwd=2,ylim=ylim,
       main="Coherence",xlab="Periodicity",ylab="Power",log="x",
       sub="Maximum Entropy Method")
  grid()
  lines(range(1/t[plus]),rep(0.5*ylim[2],2),col="red",lty=2)
  lines(tau[plus],
        0.5*atan2(Im(coh[plus]),Re(coh[plus]))/pi *
        (ylim[2]-ylim[1]) + mean(ylim),
        col="red")
  axis(4,at=0.5*seq(-pi,pi,length=5)/pi * (ylim[2]-ylim[1]) + mean(ylim),
       labels=180*seq(-pi,pi,length=5)/pi,col="red")
  lines(tau[plus],log(abs(coh[plus])),lwd=2)
  mtext("Phase (degrees)",4,col="red")
}
invisible(coh)
}

