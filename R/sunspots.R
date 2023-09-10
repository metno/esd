#' sunspots
#' 
#' Reads up-to-date sunspot data from theInternet.
#' 
#' 
#' @param url 
#'
#' @return The call returns a zoo object
#'
#' @examples
#' 
#' ## S3 method for class 'station'
#' R <- sunspots()
#' plot(R)
#'
#' @export 
sunspots <- function(url='https://www.sidc.be/SILSO/DATA/SN_d_tot_V2.0.txt') {
  Raw <- readLines(url)
  for (i in 1:3) Raw <- gsub('  ',' ',Raw)
  R <- strsplit(Raw,split=' ')
  t <- as.Date(unlist(lapply(R,function(x) paste(x[1],x[2],x[3],sep='-'))))
  x <- as.numeric(unlist(lapply(R,function(x) x[5])))
  x[x < 0] <- NA
  ok <- is.finite(t)
  RN <- zoo(x[ok],order.by=t[ok])
  attr(RN,'url') <- url
  return(RN)
}