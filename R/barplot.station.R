#' create a barplot
#'
#' @param x a 'station' object
#' @param threshold threshold - midline of plot
#'
#' @importFrom graphics barplot title
#'
#' @export
barplot.station <- function(x,threshold=0,...) {
  stopifnot(inherits(x,'station'))
  x.above <- x.below <- x
  x.above[x < threshold] <- NA
  x.below[x > threshold] <- NA
  ylim <- range(pretty(coredata(x)),na.rm=TRUE)
  barplot(as.numeric(x),col='white',ylim=ylim,border=NA)
  barplot(as.numeric(x.above),col='red',names.arg=year(x),
          ylab=paste(varid(x),'[',attr(x,'unit'),']'),axes=FALSE,
          border=NA,add=TRUE)
  barplot(as.numeric(x.below),col='blue',axes=FALSE,border=NA,add=TRUE)
  title(toupper(loc(x)))
}