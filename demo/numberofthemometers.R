## Thermometer statistics
## Rasmus.Benestad@met.no

library(esd)
ss <- select.station(param='tmin',src=c('ghcnd','ecad'))
n <- rep(NA,217)
i <- 1
for (y in 1801:2017) {n[i] <- sum( (as.numeric(ss$start) <= y) & (as.numeric(ss$end) >=y),na.rm=TRUE ); i <- i + 1}
plot(1801:2017,n,main='Number of themometers (daily minimum, tmin)',xlab='')
grid()

duration <- as.numeric(ss$end) - as.numeric(ss$start) + 1

hist(duration, breaks= seq(0,350,by=10),xlim=c(0,300),
     main='Number of long records',xlab='years',col='blue',border='grey',lwd=3)
grid()

ok <- is.finite(duration)
col <- rgb(0.6,0,0.1,duration[ok]/max(duration,na.rm=TRUE))
plot(ss$longitude[ok],ss$latitude[ok],pch=19,col=col,cex=0.3)
data(geoborders)
lines(geoborders)