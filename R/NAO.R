NAO <- function(freq="monthly", url=NULL, header=FALSE) {
  if(verbose) print("NAO")
  if(is.null(url)) {
    if(freq=="daily") {
      url <- "ftp://ftp.cdc.noaa.gov/Public/gbates/teleconn/nao.reanalysis.t10trunc.1948-present.txt"
    } else {
      url <- 'http://www.cpc.ncep.noaa.gov/products/precip/CWlink/pna/norm.nao.monthly.b5001.current.ascii.table'
    }
  }
  nao <- read.table(url,header,fill=TRUE)
  if(ncol(nao)==4) {
    d <- as.Date(paste(nao[,1],nao[,2],nao[,3],sep="-"))
    naoi <- zoo(nao[,4], order.by=d)
  } else if(ncol(nao)==13) {
    d <- as.Date(paste(sort(rep(nao[,1],12)),1:12,'01',sep='-'))
    naoi <- zoo(c(unlist(t(nao[,2:13]))), order.by=d)
  } else {
    print("Warning! Unknown file structure of downloaded NAO data")
  }
  naoi <- as.station(naoi,loc='NAOI',param='NAOI',
                     unit='dimensionless')
  attr(naoi,'url') <- url
  return(naoi)
}
