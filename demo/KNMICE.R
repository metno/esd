## Read daily station data from the KNMI Climate Explorer
## This script works for the KNMI ClimateExplorer https://climexp.knmi.nl/
## - go to the webpage, select daily stations and then a location
## - then select 'raw data' and use the URL for this page as the string for 
##   the url argument in this function
## REB 2025-06-25

KNMICE <- function(url) {
  x <- read.table(url)
  header <- sub('# ','',readLines(url,n=25))
  x <- zoo(x$V4,order.by=as.Date(paste(x$V1,x$V2,x$V3,sep='-')))
  param <- try(tolower(strsplit(header[grep('\\[',header)],split=' ')[[1]][1]))
  if (inherits(param,'try-error')) browser()
  lon <- as.numeric(sub(' degrees_east','',sub('longitude :: ','',
                                               header[grep('longitude',header)])))
  lat <- as.numeric(sub(' degrees_north','',sub('latitude :: ','',
                                               header[grep('latitude',header)])))
  alt <- as.numeric(sub('elevation :: ','',header[grep('elevation',header)]))
  stid <- sub('station_code :: ','',header[grep('station_code',header)])
  loc <- sub('station_name :: ','',header[grep('station_name',header)])
  cntr <- sub('station_country ::','',header[grep('station_country',header)])
  src <- header[grep('source_doi',header)]
  unit <- sub('\\]','',sub('\\[','',strsplit(header[grep('\\[',header)],split=' ')[[1]][2]))
  x <- as.station(x,param=variable,unit=unit,lon=lon,lat=lat,alt=alt,stid=stid,
                  loc=loc,cntr=cntr,src=src,history=paste(header[18],header[17]),
                  ref=header[7],url=url)
  print(paste(param,unit,loc,cntr,lon,lat,alt,stid))
  class(x) <- c('station','day','zoo')
  return(x)
}

seattle <- KNMICE('https://climexp.knmi.nl/data/xgdcnUSW00024281.dat')
sfo <- KNMICE('https://climexp.knmi.nl/data/xgdcnUSC00047767.dat')
la <- KNMICE('https://climexp.knmi.nl/data/xgdcnUSW00093134.dat')
kansas <- KNMICE('https://climexp.knmi.nl/data/xgdcnUSC00234379.dat')
dallas <- KNMICE('https://climexp.knmi.nl/data/xgdcnUSC00412243.dat')
atlanta <- KNMICE('https://climexp.knmi.nl/data/xgdcnUSW00053819.dat')
houston <- KNMICE('https://climexp.knmi.nl/data/xgdcnUSC00414326.dat')
boston <- KNMICE('https://climexp.knmi.nl/data/xgdcnUSW00014739.dat')
philadelphia <- KNMICE('https://climexp.knmi.nl/data/xgdcnUSW00013779.dat')
miami <- KNMICE('https://climexp.knmi.nl/data/xgdcnUSC00088396.dat')
ny <- KNMICE('https://climexp.knmi.nl/data/xgdcnUSC00308721.dat')
guadalajara <- KNMICE('https://climexp.knmi.nl/data/xgdcnMXM00076612.dat')
mexicomity <- KNMICE('https://climexp.knmi.nl/data/xgdcnMXM00076680.dat')
monterrey <- KNMICE('https://climexp.knmi.nl/data/xgdcnMXM00076393.dat')
vancouver <- KNMICE('https://climexp.knmi.nl/data/xgdcnCA001101200.dat')
toronto <- KNMICE('https://climexp.knmi.nl/data/xgdcnCA006158355.dat')

Y <- combine.stations(seattle,sfo,la,kansas,dallas,atlanta,houston,boston,philadelphia,
                      miami,ny, guadalajara,mexicomity,monterrey,vancouver,toronto)
map(Y,FUN='max',add.text = TRUE,colbar=list(show=TRUE,pal='t2m.ipcc'))
