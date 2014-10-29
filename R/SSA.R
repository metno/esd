SSA <- function(x,m=12,plot=TRUE,main="SSA analysis",sub="",
                anom=TRUE,i.eof=1) {
# After von Storch & Zwiers (1999), Statistical Analysis in Climate Research, p. 312.  

  if ( (class(x)[1] != "station") & (class(x)[1] != "eof") )
    stop('SSA: need a station object or EOF object')
  x.mean <- 0

  if (class(x)[2] == "monthly.station.record") param <- "val" else
  if (class(x)[1] == "eof") param <- "PC[,i.eof]"

  if (anom) x <- anomaly(x) 
  nt <- length(x)
  
  Nm <- nt - m + 1
  X <- matrix(rep(NA,Nm*m),Nm,m)
  #print(dim(X))
  for (i in 1:m) {
    ii <- i:(Nm-i+1)
    X[ii,i] <- coredata(x[ii])
  }

  browser() 
  ssa <- svd(X) 

  if (sub=="") sub <- paste("Window width=",m)

  ssa$m<- m; ssa$Nm <- Nm; ssa$nt <- nt
  ssa$anom <- anom
  ssa$param <- x$param
  ssa$x <- x
  class(ssa) <- c("ssa",class(x))
  if (plot) plotSSA(ssa)
  invisible(ssa)
}


