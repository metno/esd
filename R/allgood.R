# Rasmus Benestad
# A small function that removes stations with missing values from a group
# Useful before performing PCA.
allgood <- function(x,miss=.1,verbose=TRUE) {

  d <- dim(x)
  mfcol0 <- par()$mfcol
  if (verbose) {
    par(mfrow=c(2,1))
    diagram(x,main='Original data')
  }
  nok <- apply(coredata(x),1,nv)
  x <- subset(x,it=(1:d[1])[nok>=(1-miss)*d[2]])
  d <- dim(x)
  nok <- as.numeric(apply(coredata(x),2,nv))
  is <- (1:d[2])[is.element(nok,d[1])]
  y <- subset(x,is=is)
  if (verbose) {
    print("Removed stations ...")
    print(loc(x)[!is.element(loc(x),loc(y))])
    diagram(y,main='Filtered data')
  }
  par(mfcol=mfcol0)
  return(y)
}
