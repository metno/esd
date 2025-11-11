#' create a barplot
#'
#' @param height a 'station' object
#' @param threshold threshold - midline of plot
#' @param col.above color for values above threshold. Default: red.
#' @param col.above color for values below threshold. Default: blue.
#' @param verbose a boolean; if TRUE print information about progress
#' @param \dots additional arguments
#'
#' @importFrom graphics barplot title
#'
#' @exportS3Method
#' @export
barplot.station <- function(height,...,col.above='red',col.below='blue',threshold=0,verbose=FALSE) {
  if(verbose) print("barplot.station")
  stopifnot(inherits(height,'station'))
  height.above <- height.below <- height
  height.above[height < threshold] <- NA
  height.below[height > threshold] <- NA
  ylim <- range(pretty(coredata(c(height, threshold))),na.rm=TRUE)
  if(threshold!=0) {
    dummyaxis <- pretty(ylim, n = 5)
    realaxis <- dummyaxis - threshold
    barplot(as.numeric(height) - threshold, col='white', ylim=ylim - threshold, 
            border=NA, ylab=paste(varid(height),'[',attr(height,'unit'),']'), 
            yaxt="n")
    barplot(as.numeric(height.above)-threshold, col=col.above, names.arg=year(height),
            axes=FALSE, border=NA, add=TRUE)
    barplot(as.numeric(height.below) - threshold, col=col.below, 
            axes=FALSE, border=NA, add=TRUE)
    axis(side = 2, at=realaxis, labels = dummyaxis)
  } else {
    barplot(as.numeric(height),col='white',ylim=ylim,border=NA,
            ylab=paste(varid(height),'[',attr(height,'unit'),']'))
    barplot(as.numeric(height.above),col=col.above,names.arg=year(height),
            axes=FALSE,border=NA,add=TRUE)
    barplot(as.numeric(height.below),col=col.below,axes=FALSE,border=NA,add=TRUE)
  }
  title(toupper(loc(height)))
}