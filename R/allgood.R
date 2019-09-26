#' A small function that removes stations with missing values from a group
#'
#' Useful before performing PCA.
#'
#' @param x a \code{station} object
#' @param miss fraction of data that may be missing, e.g., if miss=.1 then stations with than 10\% missing data are removed
#' @param verbose a boolean; if TRUE print information about progress
#'
#' @export
allgood <- function(x,miss=.1,verbose=FALSE) {

  d <- dim(x)
  if (verbose) {
    mfcol0 <- par()$mfcol
    par(mfrow=c(2,1))
    diagnose(x,main='Original data')
  }

  ## Remove stations with little data:
  nok <- as.numeric(apply(coredata(x),2,nv))


  is <- (1:d[2])[nok>=0.7*d[2]]
  if (length(is)>0) y <- subset(x,is=is) else
                    y <- x
  ## Remove time with missing data 
  nok <- apply(coredata(y),1,nv)
  it <- (1:d[1])[nok>=(1-miss)*d[2]]
  if (verbose) print(table(year(y)[it]))
  if (length(it)>0) z <- subset(y,it=it) else
                    z <- y
  d <- dim(z)
  nok <- as.numeric(apply(coredata(z),2,nv))
  is <- (1:d[2])[is.element(nok,d[1])]
  if (length(is)>0) w <- subset(z,is=is) else
                    w <- z
  if (verbose) {
    print("Removed stations ...")
    print(loc(x)[!is.element(loc(x),loc(w))])
    diagnose(w,main='Filtered data')
    par(mfcol=mfcol0)
  }
  return(w)
}