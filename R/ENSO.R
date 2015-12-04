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

## Southern Oscillation Index
SOI <- function(url='ftp://ftp.bom.gov.au/anon/home/ncc/www/sco/soi/soiplaintext.html',header=FALSE) {
  SOI <- readLines(url)
  SOI <- gsub('\t','  ',SOI)
  SOI <- gsub('<br>','',SOI)
  SOI <- gsub('<title>','',SOI)
  SOI <- gsub('<head>','',SOI)
  SOI <- gsub('<body>','',SOI)
  SOI <- gsub('</title>','',SOI)
  SOI <- gsub('</head>','',SOI)
  SOI <- gsub('</body>','',SOI)
  i1 <- grep('year',tolower(SOI))+1
  i2 <- grep('/pre',tolower(SOI))-1
  writeLines(SOI[i1:i2],con='SOI.txt')
  soi <- read.table('SOI.txt',na.string='*',skip=13)
  soi <- zoo(c(unlist(t(soi[,2:13]))),
              order.by=as.Date(paste(sort(rep(soi[,1],12)),1:12,'01',sep='-')))
  soi <- as.station(soi,loc='SOI',param='SOI',unit='dimensionless',url=url,
                    ref=paste(SOI[7:9],collapse=', '),longname=SOI[3])
  return(soi)
}
