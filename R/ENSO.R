NINO3.4 <- function(url='ftp://ftp.cpc.ncep.noaa.gov/wd52dg/data/indices/ersst3b.nino.mth.ascii',header=TRUE,
url2='http://www.cpc.ncep.noaa.gov/data/indices/sstoi.indices') {
  enso <- read.table(url,header=header)
  nino3.4 <- zoo(enso[10],
                 order.by=as.Date(paste(enso$YR,enso$MON,'01',sep='-')))
  nino3.4 <- as.station(nino3.4,loc='Nino3.4',param='Nino3.4',
                        unit='dimensionless')
  ## Combine with more updated data from url2 (which do not extend far back in time)                     
  if (!is.null(url2)) {
    y <- read.table(url2,header=TRUE)
    y <- zoo(y$ANOM.3,order.by=as.Date(paste(y$YR,y$MON,'01',sep='-')))
    y <- as.station(y,loc='Nino3.4',param='Nino3.4',lon=c(-170,-120),lat=c(-5,5),
                    unit='dimensionless')
    nino3.4 <- combine(nino3.4,y)                
  }
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

# CSIRO
# http://www.cmar.csiro.au/sealevel/GMSL_SG_2011_up.html
GSL <- function(url='https://www.epa.gov/sites/production/files/2016-08/sea-level_fig-1.csv') {

    sl <- read.csv(url,skip=6,header=TRUE)
    zsl <- zoo(sl[[2]]*2.54,order.by=sl[[1]])
    sl <- as.station(zsl,loc=NA,param='height',unit='cm',
                     lon=NA,lat=NA,alt=NA,
                     cntr=NA,longname='CSIRO - Adjusted global mean sea level',
                     stid=NA,quality=NA,src='CSIRO, 2015',url=url,
                     reference="EPA's Climate Change Indicators in the United States: www.epa.gov/climatechange",info=NA, method= NA)
    return(sl)
}
  
AMO <- function(url='http://www.esrl.noaa.gov/psd/data/correlation/amon.us.long.data') {
  amo.test <- readLines(url)
  nrows <- sum(is.element(nchar(amo.test),max(nchar(amo.test))))
  amo <- read.table(url,skip=1,nrows=nrows)
  amo[amo <= -99] <- NA
  amo <- zoo(c(t(as.matrix(amo[2:13]))),
             order.by=as.Date(paste(sort(rep(amo$V1,12)),rep(1:12,length(amo$V1)),'01',sep='-')))
  amo <- as.station(amo,loc=NA,param='AMO',unit='dimensionless',
                     lon=NA,lat=NA,alt=NA,
                     cntr=NA,longname='Atlantic Multi-decadal Oscillation unsmoothed from the Kaplan SST V2',
                     stid=NA,quality=NA,src='Calculated at NOAA/ESRL/PSD1',url=url,
                     reference=NA,info=NA, method= NA)
  return(amo)
}

QBO <- function(url='http://www.esrl.noaa.gov/psd/data/correlation/qbo.data') {
  qbo.test <- readLines(url)
  nrows <- length(qbo.test) - 6
  qbo <- read.table(url,skip=1,nrows=nrows)
  qbo[qbo <= -999] <- NA
  qbo <- zoo(c(t(as.matrix(qbo[2:13]))),
             order.by=as.Date(paste(sort(rep(qbo$V1,12)),
                               rep(1:12,length(qbo$V1)),'01',sep='-')))
  amo <- as.station(qbo,loc=NA,param='QBO',unit='dimensionless',
                     lon=NA,lat=NA,alt=NA,
                     cntr=NA,longname='Quasi-biennial oscillation',
                     stid=NA,quality=NA,src='Calculated at NOAA/ESRL PSD',
                     url=url,reference=NA, method= NA,
               info='http://www.esrl.noaa.gov/psd/data/climateindices/list/')
  return(qbo)
}

## Central England Temperature
CET <- function(url='http://hadobs.metoffice.com/hadcet/cetml1659on.dat') {
  cet <- read.table(url,skip=8)
  cet[cet <= -99] <- NA
  cet <- zoo(c(t(as.matrix(cet[2:13]))),
             order.by=as.Date(paste(sort(rep(cet$V1,12)),rep(1:12,length(cet$V1)),'01',sep='-')))
  cet <- as.station(cet,loc=NA,param='t2m',unit='degC',
                    lon=NA,lat=NA,alt=NA,
                    cntr=NA,longname='Central England Temperature',
                    stid=NA,quality=NA,src=NA,url=url,
                    reference='Parker, et. al. (1992), Int. J. Clim.',info=NA, method= NA)
  return(cet)
}

CO2 <- function(url='ftp://aftp.cmdl.noaa.gov/products/trends/co2/co2_mm_mlo.txt') {
  X <- read.table(url)
  X[X <= -99] <- NA
  co2 <- zoo(X$V4,order.by=as.Date(paste(X$V1,X$V2,'01',sep='-')))
  co2 <- as.station(co2,loc='Mauna Loa',param=expression(C*O[2]),unit='ppm',
                     lon=-155.5763,lat=19.5362,alt=3397,
                     cntr=NA,longname='Carbon dioxide',
                     stid=NA,quality=NA,src='NOAA/ESRL',
                     url=url,reference='C.D. Keeling, R.B. Bacastow, A.E. Bainbridge, C.A. Ekdahl, P.R. Guenther, and L.S. Waterman, (1976), Atmospheric carbon dioxide variations at Mauna Loa Observatory, Hawaii, Tellus, vol. 28, 538-551', method= NA,
               info='Use of these data implies an agreement to reciprocate.')
  return(co2)
}



