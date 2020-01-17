#' Bias correct
#'
#' Bias correction as described in  Imbert & Benestad (2005), Theor. Appl. Clim., DOI: 10.1007/s00704-005-0133-4
#'
#' @export
biasfix <- function(x,verbose=FALSE) {
  if(verbose) print("biasfix")
  stopifnot(!missing(x), inherits(x,"eof"),inherits(x,"comb"))
  if (!is.null(attr(x,'diagnose'))) return(x)
  diag <- diagnose(x)
  n <- attr(x,'n.apps')
  for ( i in 1:n ) {
    z <- NULL
    eval(parse(text=paste("z <- attr(x,'appendix.",i,"')",sep="")))
    Z <- coredata(z) 
    year.ox <- range(year(z)[is.element(year(z),year(x))])
    sd.o <- apply(coredata(subset(x,it=range(year.ox))),2,sd,na.rm=TRUE)
    sd.x <- apply(coredata(subset(z,it=range(year.ox))),2,sd,na.rm=TRUE)
    scalfac <- sd.o/sd.x
    me.o <- apply(coredata(subset(x,it=range(year.ox))),2,mean,na.rm=TRUE)
    me.x <- apply(coredata(subset(z,it=range(year.ox))),2,mean,na.rm=TRUE)
    mean.diff <- me.o - me.x
    
    ## diagnose: (1 + sd(z))/(1 + sd(x))
    ## x is reanalysis; z is gcm:
    for (j in 1:length(Z[1,]))
      Z[,j] <- scalfac[j]*Z[,j] + mean.diff[j]
    y <- zoo(Z,order.by=index(z))
    y <- attrcp(z,y)
    eval(parse(text=paste("y -> attr(x,'appendix.",i,"')",sep="")))
  }
  attr(x,'history') <- history.stamp(x)
  attr(x,'quality') <- "'bias' corrected -  ref (Imbert & Benestad (2005); Theor. Appl. Clim.; DOI: 10.1007/s00704-005-0133-4"
  attr(x,'baseline') <- year.ox
  attr(x,'diagnose') <- diag
  invisible(x)
}
