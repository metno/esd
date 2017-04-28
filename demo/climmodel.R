## Rasmus.Benestad@met.no
## climmodel - Model of the meanseasonal variations

x <- station(stid=18700,src='metnod',user='metno')
yrs <- rownames(table(year(x)))
ny <- length(yrs)
ac <- rep(0,365)
par(bty='n')
plot(c(1,365),range(x,na.rm=TRUE),type='n',ylab='C',xlab='')

n <- 0
for (i in 1:ny) {
  z <- coredata(window(x,start=as.Date(paste(yrs[i],'-01-01',sep='')),
                       end=as.Date(paste(yrs[i],'-12-31',sep=''))))
  if (sum(is.finite(z)) >= 365) {
    y <- approx(1:length(z),z,xout=1:365)$y
    lines(1:365,y,col=rgb(0,0,0,0.01))
    ac <- ac + y; n <- n + 1
  } else if (i==ny) {
    y <- coredata(window(x,start=as.Date(paste(yrs[i],'-01-01',sep='')),end=end(x)))
                      #end=as.Date(paste(yrs[i],'-03-01',sep=''))))
    
  }
}
ac <- ac/n
cal <- data.frame(y=ac,x1=cos(2*pi*seq(0,1,length=365)),x2=sin(2*pi*seq(0,1,length=365)),
                       x3=cos(4*pi*seq(0,1,length=365)),x4=sin(4*pi*seq(0,1,length=365)))
model <- lm(y ~x1 + x2 + x3 + x4, data=cal)
lines(1:365,predict(model),lwd=4,col=rgb(1,0.5,0.5,0.5))
lines(1:length(y),y,col=rgb(0,0,0,0.5),lwd=2)
lines(seq(as.Date(paste(yrs[i],'-03-25',sep='')),end(x),by='1 day') - as.Date(paste(yrs[i],'-01-01',sep=''))+1,
      coredata(window(x,start=as.Date(paste(yrs[i],'-03-25',sep='')),end=end(x))),lwd=2,col=rgb(0,0,0,0.5))
lines(seq(as.Date(paste(yrs[i],'-03-25',sep='')),end(x),by='1 day') - as.Date(paste(yrs[i],'-01-01',sep=''))+1,
      coredata(trend(window(x,start=as.Date(paste(yrs[i],'-03-25',sep='')),end=end(x)))),lwd=2,col='red')

legend(150,-10,c('modell','observasjon'),col=c(rgb(1,0.5,0.5,0.5),'black'),
       lwd=c(4,1),bty='n')