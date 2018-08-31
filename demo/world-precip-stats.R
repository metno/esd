## Draw map of precipitation statistics:

colscheme <- function(z,d=0.1,w=0.3) {
  z[!is.finite(z)] <- 0
  z <- (z - min(z,na.rm=TRUE))/diff(range(z,na.rm=TRUE))
  wet <- z > w; mod <- (z >= d) & (z <= w); dry <- z < d
  colwet <- rgb(1-z[wet],1-z[wet],1,0.5)
  coldry <- rgb(2*z[dry]+0.5*(1-d),0.5*(2*z[dry]+0.9*(1-d)),0,0.5)
  colmod <- rgb(0.3,1.5*z[mod]+0.3,0.3,0.5)
  col <- rgb(z,z,z,0.5)
  col[wet] <- colwet; col[mod] <- colmod; col[dry] <- coldry
  invisible(col)
}

if (!file.exists('qqt5yr.world.rda'))
  download.file(url = 'https://ndownloader.figshare.com/files/5668374',destfile = 'qqt5yr.world.rda')
load('qqt5yr.world.rda')

dev.new(fig.width=8,fig.height=12)
par(bty='n',mfcol=c(3,1),xaxt='n',yaxt='n',mar=c(1,1,1,1),xpd=NA)

## Mean precipitation
z <- apply(qqt5yr.world$XBAR,1,mean,na.rm=TRUE)*365.25
z[z > 5000] <- NA 

col <- colscheme(z)
plot(qqt5yr.world$lon,qqt5yr.world$lat,col=col,pch=10,cex=0.5,
     xlab='',ylab='')
data(geoborders)
lines(geoborders)
text(-180,85,expression(bar(x)),cex=1.5)

## Wet-day frequency
z <- apply(qqt5yr.world$FW,1,mean,na.rm=TRUE)
col <- colscheme(z)
plot(qqt5yr.world$lon,qqt5yr.world$lat,col=col,pch=10,cex=0.5,
     xlab='',ylab='')
data(geoborders)
lines(geoborders)
text(-180,85,expression(f[w]),cex=1.5)

## Wet-day mean precipitation
z <- apply(qqt5yr.world$MU,1,mean,na.rm=TRUE)
z[z > 50] <- NA
col <- colscheme(z)
plot(qqt5yr.world$lon,qqt5yr.world$lat,col=col,pch=10,cex=0.5,
     xlab='',ylab='')
data(geoborders)
lines(geoborders)
text(-180,85,expression(mu),cex=1.5)
