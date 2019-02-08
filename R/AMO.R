AMO <- function(url=NULL, verbose=FALSE) {
  if(verbose) print("AMO")
  if(is.null(url)) url <- 'http://www.esrl.noaa.gov/psd/data/correlation/amon.us.long.data'
  amo.test <- readLines(url)
  nrows <- sum(is.element(nchar(amo.test),max(nchar(amo.test))))
  skip <- min(which(is.element(nchar(amo.test), max(nchar(amo.test)))))-1
  amo <- read.table(url, skip=skip, nrows=nrows)
  amo[amo <= -99] <- NA
  if(ncol(amo)==13) {
    d <- as.Date(paste(sort(rep(amo$V1,12)),rep(1:12,length(amo$V1)),'01',sep='-'))
    amo <- zoo(c(t(as.matrix(amo[2:13]))), order.by=d)
  } else if(ncol(amo)==2) {
    
  }
  amo <- as.station(amo,loc=NA,param='AMO',unit='dimensionless',
                    lon=NA,lat=NA,alt=NA,
                    cntr=NA,longname='Atlantic Multi-decadal Oscillation unsmoothed from the Kaplan SST V2',
                    stid=NA,quality=NA,src='Calculated at NOAA/ESRL/PSD1', url=url,
                    reference=NA,info=NA, method= NA)
  return(amo)
}
