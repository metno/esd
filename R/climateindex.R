#' Functions to read climate indices from sources on the internet
#'
#' NAO: Daily North Atlantic Oscillation index from NCEP/NOAA 
#'
#' NINO3.4: Daily Nino3.4 index provided by NOAA/NCDC downloaded from KNMI Climate Explorer
#'
#' SOI: Souther Oscillation index from The Australian Bureau fo Meteorology (bom.gov.au)
#'
#' GSL: Global Average Absolute Sea level Change, 1880-2015,
#' from EPA's Climate Change Indicators in the United States: www.epa.gov/climate-indicators
#'
#' GSL.nasa: Global Average Sea level from NASA
#' GSL.aviso: Global Average Sea level from AVISO
#'
#' QBO: Quasi-Biennial Oscillation. Calculated at NOAA/ESRL/PSD from the zonal average
#' of the 30mb zonal wind at the equator as computed from the NCEP/NCAR Reanalysis.
#'
#' CET: Central England Temperature from the Hadley Center
#'
#' CO2: Carbon dioxide at Mauna Loa, Hawaii, from NOAA ESRL.
#' Reference: 'C.D. Keeling, R.B. Bacastow, A.E. Bainbridge, C.A. Ekdahl, P.R. Guenther, and L.S. Waterman, (1976),
#' Atmospheric carbon dioxide variations at Mauna Loa Observatory, Hawaii, Tellus, vol. 28, 538-551'
#'
#' AMO: Atlantic Multidecadal Oscillation, unsmoothed calculated from the Kaplan SST V2 at NOAA/ESRL/PSD1
#'
#' IOD: Indian Ocean Dipole index
#' 
#' Sunspots: updated monthly sunspot number
#' 
#' @aliases NAO NINO3.4 SOI GSL GSL.nasa QBO CET CO2 AMO IOD Sunspots
#'
#' @param freq frequency
#' @param url a URL or web address to location of data
#' @param header a boolean, indicating if file includes a header or not
#' @param verbose a boolean; if TRUE print information about progress
#'
#' @export
NAO <- function(freq="monthly", url=NULL, header=FALSE, verbose=FALSE) {
  if(verbose) print("NAO")
  if(is.null(url)) {
    if(freq=="daily") {
      url <- "ftp://ftp.cpc.ncep.noaa.gov/cwlinks/norm.daily.nao.index.b500101.current.ascii"
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

#' @export
NINO3.4 <- function(url=NULL, header=TRUE, freq="monthly", verbose=FALSE) {
  if(verbose) print("NINO3.4")
  if(is.null(url)) {
    if(freq=="daily") {
      url <- 'https://climexp.knmi.nl/data/inino34_daily.dat'
      header <- TRUE
    } else {
      url <- 'https://climexp.knmi.nl/data/inino5.dat'
      header <- FALSE
    }
  }
  enso <- read.table(url,header=header)
  if(ncol(enso)==2) {
    d <- as.Date(strptime(enso[,1],format="%Y%m%d"))
    nino3.4 <- zoo(enso[,2], order.by=d)
  } else if(ncol(enso)==13) {
    d <- as.Date(paste(sort(rep(enso[,1],12)),1:12,'01',sep='-'))
    nino3.4 <- zoo(c(unlist(t(enso[,2:13]))), order.by=d)
  } else if (all(c("YR","MON","ANOM.3") %in% colnames(enso))) {
    d <- as.Date(paste(enso$YR,enso$MON,'01',sep='-'))
    nino3.4 <- zoo(enso$ANOM.3, order.by=d)
  } else {
    print("Warning! Don't know how to read data from url",url)
  }
  nino3.4 <- as.station(nino3.4,loc='Nino3.4',param='Nino3.4',
                        unit='dimensionless')
  attr(nino3.4,'url') <- url
  return(nino3.4)
}

#' @export
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
  soi <- read.table('SOI.txt',na.strings='*',skip=13)
  soi <- zoo(c(unlist(t(soi[,2:13]))),
              order.by=as.Date(paste(sort(rep(soi[,1],12)),1:12,'01',sep='-')))
  soi <- as.station(soi,loc='SOI',param='SOI',unit='dimensionless',url=url,
                    ref=paste(SOI[7:9],collapse=', '),longname=SOI[3])
  return(soi)
}

# CSIRO
# http://www.cmar.csiro.au/sealevel/GMSL_SG_2011_up.html
#' @export
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

#' @export
GSL.nasa <- function(url='ftp://podaac.jpl.nasa.gov/allData/merged_alt/L2/TP_J1_OSTM/global_mean_sea_level/GMSL_TPJAOS_V4_199209_201708.txt',is=9) {
  
  sl <- read.table(url,skip=6,header=TRUE,comment.char='H')
  param <- switch(as.character(is),'1'='type','2'='cycle','3'='year+fraction of year',
                  '4'='number','5'='number','6'='height','7'='error','8'='error',
                  '9'='height','10'='error','11'='height','12'='height')
  longname=switch(as.character(is),'1'='altimeter type','2'='merged file cycle','3'='year+fraction of year',
                  '4'='number of observations','5'='number of weighted observations',
                  '6'='GMSL (Global Isostatic Adjustment (GIA) not applied)',
                  '7'='standard deviation of GMSL (GIA not applied)',
                  '8'='smoothed (60-day Gaussian type filter) GMSL (GIA not applied)',
                  '9'='GMSL (Global Isostatic Adjustment (GIA) applied)',
                  '10'='standard deviation of GMSL (GIA applied)',
                  '11'=' smoothed (60-day Gaussian type filter) GMSL (GIA applied)',
                  '12'='smoothed (60-day Gaussian type filter) GMSL (GIA applied) anomaly')
  unit=switch(as.character(is),'1'='rtype','2'='cycle','3'='year',
              '4'='number','5'='number','6'='mm','7'='mm','8'='mm',
              '9'='mm','10'='mm','11'='mm','12'='mm')
  hdr <- readLines(url,n=42)
  yr <- trunc(sl[[3]]); mo <- trunc(12*(sl[[3]] - trunc(sl[[3]]))) %% 12 + 1;
  dy <- trunc(365.25*(sl[[3]] - yr - (mo-1)/12))
  time <- as.Date(paste(yr,mo,dy,sep='-'))
  bad0 <- is.na(time) & (dy ==0)
  bad29 <- is.na(time) & (dy >=29)
  ## Fudge - quick fix for bad dates:
  #print(time); browser()
  time[bad0] <-  as.Date(paste(yr[bad0],mo[bad0],1,sep='-'))
  time[bad29] <-  as.Date(paste(yr[bad29],mo[bad29],dy[bad29]-2,sep='-'))

  zsl <- zoo(sl[[is]],order.by=time)
  sl <- as.station(zsl,loc=NA,param=param,unit=unit,
                   lon=NA,lat=NA,alt=NA,
                   cntr=NA,longname=longname,
                   stid=NA,quality=NA,src='NASA',url=url,
                   reference="NASA's Climate Change Indicators",info=hdr, method= NA)
  return(sl)
}

GSL.aviso <- function(url='ftp://ftp.aviso.altimetry.fr/pub/oceano/AVISO/indicators/msl/MSL_Serie_MERGED_Global_AVISO_GIA_Adjust_Filter2m.txt') {
  gsl <- read.table(url)
  yr <- trunc(gsl[,1])
  t <- as.Date(julian(as.Date(paste0(yr,'-01-01'))) + 365.25*(gsl[,1]- yr))
  z <- zoo(gsl[,2],order.by=t)
  return(z)
}

#' @export
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
#' @export
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

#' @export
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

#' @export
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

#' @export
IOD <- function(url='https://www.esrl.noaa.gov/psd/gcos_wgsp/Timeseries/Data/dmi.long.data') {
  ## Indian Ocean Dipole
  test <- readLines(url)
  n <- length(test)
  X <- as.matrix(read.table(url,skip=1,nrow = n-8))
  X[X <= -999] <- NA
  yr <- sort(rep(X[,1],12))
  mo <- rep(1:12,length(X[,1]))
  iod <- c(t(X[,2:13]))
  y <- zoo(iod,order.by=as.Date(paste(yr,mo,'01',sep='-')))
  return(y)
}

#'@export
Sunspots <- function(url='http://sidc.oma.be/silso/DATA/SN_m_tot_V2.0.txt') {
  S <- read.table(url,comment.char = "*")
  s <- zoo(S$V4,order.by=as.Date(paste(S$V1,S$V2,'01',sep='-')))
  attr(s,'url') <- url
  attr(s,'variable') <- 'Sunspots'
  attr(s,'unit') <- 'monthly sum'
  return(s)
}
  


