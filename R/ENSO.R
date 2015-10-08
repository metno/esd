NINO3.4 <- function(url='ftp://ftp.cpc.ncep.noaa.gov/wd52dg/data/indices/ersst3b.nino.mth.ascii',header=TRUE) {
  enso <- read.table(url,header=header)
  nino3.4 <- zoo(enso[10],
                 order.by=as.Date(paste(enso$YR,enso$MON,'01',sep='-')))
  nino3.4 <- as.station(nino3.4,loc='Nino3.4',param='Nino3.4',
                        unit='dimensionless')
  attr(nino3.4,'url') <- url
  return(nino3.4)
}

NAO <- function(url='http://www.cpc.ncep.noaa.gov/products/precip/CWlink/pna/norm.nao.monthly.b5001.current.ascii.table',header=FALSE) {

  nao <- read.table(url,header,fill=TRUE)
  naoi <- zoo(c(unlist(t(nao[,2:13]))),
              order.by=as.Date(paste(sort(rep(nao[,1],12)),1:12,'01',sep='-')))
  naoi <- as.station(naoi,loc='NAOI',param='NAOI',
                     unit='dimensionless')
  attr(naoi,'url') <- url
  return(naoi)
}
