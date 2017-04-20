## Rasmus.Benestad@met.no
## climmodel - Model of the meanseasonal variations

data(ferder)
yrs <- rownames(table(year(ferder)))
ny <- length(yrs)
ac <- rep(0,365)
par(bty='n')
plot(c(1,365),range(ferder,na.rm=TRUE),type='n',ylab='C',xlab='')

for (i in 1:ny) {
  z <- coredata(window(ferder,start=as.Date(paste(yrs[i],'-01-01',sep='')),
                       end=as.Date(paste(yrs[i],'-12-31',sep=''))))
  y <- approx(1:length(z),z,xout=1:365)$y
  lines(1:365,y,col=rgb(0,0,0,0.01))
  ac <- ac + y
}
ac <- ac/ny
lines(1:365,ac,lwd=4,col=rgb(1,0.5,0.5,0.5))
lines(1:365,y)

legend(150,-10,c('modell','observasjon'),col=c(rgb(1,0.5,0.5,0.5),'black'),
       lwd=c(4,1),bty='n')