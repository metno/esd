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
  T2m
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
