## R-script to illustrate downscaling climate information: the parameters of pdfs.

library(esd)
data(ferder)
data(bjornholt)
mx <- mean(subset(ferder,it='jja'),na.rm=TRUE)
sx <- sd(subset(ferder,it='jja'),na.rm=TRUE)
mu <- wetmean(bjornholt)

X <- seq(0,40,by=0.1)
pdf.t1 <- dnorm(X,mean=mx,sd=sx)
pdf.t2 <- dnorm(X,mean=mx+4,sd=sx+1)

pdf.p1 <- dexp(X,rate=1/mu)
pdf.p2 <- dexp(X,rate=1/(mu+0.4*mu))

par(mfcol=c(1,2),bty='n',yaxt='n',xpd=NA)
plot(X,pdf.t1,xlab=expression(deg*C),ylab='',lwd=2,type='l',
     main='Temperature statistics')
lines(X,pdf.t2)
polygon(c(X,max(X),min(X)),c(pdf.t1,0,0),col=rgb(0.5,0.5,0.5,0.5),lwd=3,border=rgb(0.5,0.5,0.5,0.5))
polygon(c(X,max(X),min(X)),c(pdf.t2,0,0),col=rgb(0.75,0.1,0.1,0.5),lwd=3,border=rgb(0.75,0.1,0.1,0.5))

text(mx,max(pdf.t1),expression(mu[1]),col=rgb(0.5,0.5,0.5),pos=3)
text(mx+4,max(pdf.t2),expression(mu[2]),col=rgb(0.75,0.1,0.1),pos=3)
arrows(mx,max(pdf.t1),mx,0,col=rgb(0.5,0.5,0.5),length=0.1)
arrows(mx+4,max(pdf.t2),mx+4,0,col=rgb(0.75,0.1,0.1),length=0.1)
lines(mx+c(-sx,sx),approx(X,pdf.t1,mx+c(-sx,sx))$y,lty=2)
lines(mx+4+c(-sx-1,sx+1),approx(X,pdf.t2,mx+4+c(-sx-1,sx+1))$y,lty=2)
text(mx-sx,approx(X,pdf.t1,mx-sx)$y,expression(sigma[1]),pos=2)
text(mx+5+sx,approx(X,pdf.t2,mx+3-sx)$y,expression(sigma[2]),pos=4)

plot(X,pdf.p1,xlab='mm/day',ylab='',lwd=2,type='l',
     main='Precipitation statistics')
lines(X,pdf.p2)
polygon(c(X,max(X),min(X)),c(pdf.p1,0,0),col=rgb(0.5,0.5,0.5,0.5),lwd=3,border=rgb(0.5,0.5,0.5,0.5))
polygon(c(X,max(X),min(X)),c(pdf.p2,0,0),col=rgb(0.75,0.1,0.1,0.5),lwd=3,border=rgb(0.75,0.1,0.1,0.5))

text(mu,approx(X,pdf.p1,mu)$y,expression(mu[1]),col=rgb(0.5,0.5,0.5),pos=4)
text(mu+0.4*mu,approx(X,pdf.p2,mu+0.4*mu)$y,expression(mu[2]),col=rgb(0.75,0.1,0.1),pos=4)
arrows(mu,approx(X,pdf.p1,mu)$y,mu,0,col=rgb(0.5,0.5,0.5),length=0.1)
arrows(mu+0.4*mu,approx(X,pdf.p2,mu+0.4*mu)$y,mu+0.4*mu,0,col=rgb(0.75,0.1,0.1),length=0.1)