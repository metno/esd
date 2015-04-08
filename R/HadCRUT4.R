HadCRUT4 <- function(url="http://www.metoffice.gov.uk/hadobs/hadcrut4/data/current/time_series/HadCRUT.4.3.0.0.monthly_ns_avg.txt",plot=FALSE) {

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
