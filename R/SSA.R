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


