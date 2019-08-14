#' create a barplot
#'
#' @param height a 'station' object
#' @param threshold threshold - midline of plot
#' @param verbose a boolean; if TRUE print information about progress
#'
#' @importFrom graphics barplot title
#'
#' @export barplot.station
barplot.station <- function(height,...,threshold=0,verbose=FALSE) {
  if(verbose) print("barplot.station")
  stopifnot(inherits(height,'station'))
  height.above <- height.below <- height
  height.above[height < threshold] <- NA
  height.below[height > threshold] <- NA
  ylim <- range(pretty(coredata(height)),na.rm=TRUE)
  barplot(as.numeric(height),col='white',ylim=ylim,border=NA)
  barplot(as.numeric(height.above),col='red',names.arg=year(height),
          ylab=paste(varid(height),'[',attr(height,'unit'),']'),axes=FALSE,
          border=NA,add=TRUE)
  barplot(as.numeric(height.below),col='blue',axes=FALSE,border=NA,add=TRUE)
  title(toupper(loc(height)))
}