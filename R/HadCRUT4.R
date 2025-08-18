#' Download HadCRUT4 temperature data from UK MetOffice
#'
#' @param url url
#' @param plot a boolean; if TRUE show resutls in plot
#'
#' @export
HadCRUT4 <- function(url="http://www.metoffice.gov.uk/hadobs/hadcrut4/data/current/time_series/HadCRUT.4.6.0.0.monthly_ns_avg.txt",
                     plot=FALSE) {
  ## REB New URL for new version 2017-10-08
  X <- read.table(url)
  year <- as.numeric(substr(X$V1,1,4))
  month <-  as.numeric(substr(X$V1,6,7))
  T2m <- zoo(X$V2,as.Date(paste(year,month,"15",sep="-")))
  T2m <- as.station(T2m,param='t2m',unit='deg C',loc='global',
                    lon=NA,lat=NA,longname='global mean temperature',
                    ref='HadCrut4, UK Met Office',
                    url=url)
  attr(T2m,'history') <- history.stamp()
  if (plot) plot(T2m)
  return(T2m)
}

#' @export
HadCRUT5 <- function(url='https://crudata.uea.ac.uk/cru/data/temperature/HadCRUT5.0Analysis_gl.txt', plot=FALSE, sep=" ") {
## Different format to the HadCRUT4 data - annoying...
  X <- readLines(url)
  i1 <- seq(1,length(X),by=2); n <- length(i1)
  yr <- sort(rep(as.numeric(substr(X[i1],1,5)),12))
  mo <- rep(1:12,n)
  ## temperature - remove the last with the annual mean
  for (ic in seq(10,2,by=-1)) X <- gsub(paste0(rep(' ',ic),collapse=''),' ',X)
  t2m <- unlist(lapply(X[i1],function(x) as.numeric(strsplit(substr(x,7,nchar(x)),sep)[[1]][-13])))
  # print(c(length(yr),length(mo),length(t2m)))
  # print(summary(t2m))
  ## Area coverage- remove the last with the annual mean
  A <- unlist(lapply(X[i1+1],function(x) as.numeric(strsplit(substr(x,7,nchar(x)),sep)[[1]][-13])))
  #print(summary(A))
  t <- as.Date(paste(yr,mo,15,sep='-'))
  t2m[t2m <= -9.99] <- NA
  T2m <- zoo(x=t2m,order.by=t)
  A <- zoo(x=A,order.by=t)
  y <- as.station(T2m,param='t2m',unit='deg C',loc='global',
                    lon=NA,lat=NA,longname='global mean temperature',
                    ref='HadCrut4, UK Met Office',
                    url=url)
  attr(y,'history') <- history.stamp()
  attr(y,'area-cover') <- A
  if (plot) plot(T2m)
  return(y)
} 

#' Download GISS Sea Surface Temperature data from NASA
#'
#' @param url url
#' @param plot a boolean; if TRUE show resutls in plot
#'
#' @export
NASAgiss <- function(url='http://data.giss.nasa.gov/gistemp/tabledata_v3/GLB.Ts+dSST.txt',plot=FALSE) {
  olines <- readLines(url)
  ll <- nchar(olines)
  lines <- olines[ll==max(ll)]
  lh <- grep('Year',lines)
  lines <- lines[-lh]
  for (i in 1:3) lines <- gsub('  ',' ',lines)
  x <- as.numeric(unlist(strsplit(lines,' ')))
  ny <- length(lines)
  dim(x) <- c(length(gregexpr(' ',lines[1])[[1]])+1,ny)
  y <- zoo(0.01*c(x[2:13,]),order.by=as.Date(paste(sort(rep(x[1,],12)),rep(1:12,ny),'01',sep='-')))
  y <- as.station(y,param='t2m',unit='deg C',loc='global',
                    lon=NA,lat=NA,longname='global mean temperature',
                    ref='NASA/GISS',url=url)
  attr(y,'history') <- history.stamp()
  if (plot) plot(y)
  y
}

#' @export
RSS <- function(url='https://data.remss.com/msu/monthly_time_series/RSS_Monthly_MSU_AMSU_Channel_TLT_Anomalies_Land_and_Ocean_v04_0.txt') {
  x <- read.table(url,skip=3)
  rss <- zoo(x$V3,order.by=as.Date(paste(x$V1,x$V2,'01',sep='-')))
  rss <- as.station(rss,param='temperature',unit='degC',descr='-70.0/82.5',
                    ref='https://www.remss.com/research/climate/',
                    url=url,loc='pseudo-global',longname='Tropospheric Temperature',
                    src='Remote Sensing System (RSS)')
  return(rss)
}

#' @export
UAH <- function(url='https://www.nsstc.uah.edu/data/msu/v6.1/tlt/uahncdc_lt_6.1.txt') {
  test <- readLines(url)
  nrows <- grep('Trend',test) - 4
  x <- read.table(url,header=TRUE,nrows=nrows)
  uah <- zoo(x$Globe,order.by=as.Date(paste(x$Year,x$Mo,'01',sep='-')))
  uah <- as.station(uah,param='temperature',unit='degC',descr='-70.0/82.5',
                    ref='https://www.nsstc.uah.edu/climate/',
                    url=url,loc='global',longname='Tropospheric Temperature',
                    src='University of Alabama in Huntsville (UAH)')
  return(uah)
}
