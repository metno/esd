#' Singular Spectrum Analysis
#' 
#' After von Storch & Zwiers (1999), Statistical Analysis in Climate Research,
#' p. 312
#' 
#' @aliases SSA
#'
#' @param x A station or eof object.
#' @param m Window length.
#' @param plot Flag: plot the diagnostics.
#' @param main main title (see \code{link{plot}}).
#' @param sub subtitle (see \code{link{plot}}).
#' @param anom TRUE if analysis on anomalies
#' @param ip If x is an eof-object, which PC to use.
#' @param verbose Print out diagnostics.
#'
#' @return A SSA object: An \code{link{svd}} object with additional parameters:
#' m (window length), nt (original length of series), Nm (effective length of
#' series= nt - m), anom (FLAG for use of anomaly), param (name of parameter,
#' typically 'precip' or 't2m'), station (the station object to which SSA is
#' applied).
#'
#' @keywords manip
#' 
#' @export
SSA <- function(x,m=12,plot=TRUE,main="SSA analysis",sub="",
                anom=TRUE,ip=1,verbose=FALSE) {
  if(verbose) print("SSA")
  # After von Storch & Zwiers (1999), Statistical Analysis in Climate Research, p. 312.  

  if ( (class(x)[1] != "station") & (class(x)[1] != "eof") ) {
    stop('SSA: need a station object or EOF object')
  }
  x.mean <- 0

  if (class(x)[2] == "station") {
    y <- "val"
  } else if (class(x)[1] == "eof") {
    y <- "PC[,ip]"
  }

  if (anom) x <- anomaly(y) 
  nt <- length(y)
  
  if (sum(!is.finite(y))>0) {
  # need to fill in missing data
    ok <- (1:nt)[is.finite(y)]
    y <- approx(ok,y[ok],xout=1:nt)$y
  }
  
  Nm <- nt - m + 1
  X <- matrix(rep(NA,Nm*m),Nm,m)
  #if(verbose) print(dim(X))
  for (i in 1:m) {
    ii <- i:(Nm-i+1)
    X[ii,i] <- coredata(y[ii])
  }
  
  if(verbose) print(str(c(X)))
  udv <- svd(X) 

  if (sub=="") sub <- paste("Window width=",m)

  ssa <- zoo(udv$v,order.by=index(x)[1:Nm])
  attr(ssa,'pattern') <- udv$u 
  attr(ssa,'eigenvalues') <- udv$d
  attr(ssa,'m') <- m; 
  attr(ssa,'Nm') <- Nm
  attr(ssa,'nt') <- nt
  attr(ssa,'anom') <- anom
  attr(ssa,'original') <- x
  attr(ssa,'history') <- history.stamp()
  class(ssa) <- c("ssa",class(x))
  if (plot) plot(ssa)
  invisible(ssa)
}


