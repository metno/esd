# Rasmus Benestad
# A small function that removes stations with missing values from a group
# Useful before performing PCA.
allgood <- function(x) {
  d <- dim(x)
  nok <- apply(coredata(x),1,nv)
  x <- subset(x,it=(1:d[1])[nok>=0.75*d[2]])
  d <- dim(x)
  nok <- as.numeric(apply(coredata(x),2,nv))
  is <- (1:d[2])[is.element(nok,d[1])]
  y <- subset(x,is=is)
  return(y)
}
