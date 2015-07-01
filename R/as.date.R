vals <- nc$dim$time$vals
yrs <- time$vals%/%365 + 1850
mns <-

   table( (vals - 4015) %/% round(median(diff(vals))) %/% 12 + 1)

mndays <- rep(c(31,28,31,30,31,30,31,31,30,31,30,31),240)
d.lag <- lag(mndays,1)
daysfromyr <- cumsum(d.lag+mndays)/2-mndays/2 + 4015

(vals-4015)%/%365+1
months <- (cumsum(d.lag+mndays)/2-mndays/2) %/% 30 + 1
nm <- rep(1,12)

ndays <- (vals/365-vals%/%365) * 365

((vals-4015) %/% 365 ) * 365

rep(mndays,240)

## generate months
vals <- nc$dim$time$vals
l <- length(vals)/12
mo <- 1
my1 <- seq(mo,by=1,length=12-mo+1) ## months in year 1
mye <- seq(1,by=1,length= l - ((l-1)-length(my1)+1))
months <- c(my1,rep(rep(1:12),l-2),mye)
yrs <- vals%/%365 + 1850
dates <- as.Date(paste(yrs,months,"01",sep="-"))
