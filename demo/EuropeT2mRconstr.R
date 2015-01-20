# R.E. Benestad, R-script plotting the Luterbacher et al. (2004) temperature
# reconstruction for Europe.
# http://www.ncdc.noaa.gov/paleo/pubs/luterbacher2004/luterbacher2004.html

url <- "ftp://ftp.ncdc.noaa.gov/pub/data/paleo/historical/europe-seasonal.txt"
Luterbacher <- read.table(url,skip=118,header=TRUE)
info <- readLines(url,n=117)

t2m.zoo <- zoo(x=Luterbacher[,2:6],order.by=Luterbacher$Year)
t2m <- as.station(t2m.zoo,loc='Europe',param='t2m',unit="deg C",
                  longname='Luterbacher et al. (2004) temperature reconstruction for Europe',
                  url=url,info=paste(info,collapse='; '),aspect='anomaly')

plot.zoo(t2m)
