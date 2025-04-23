## Illustrate downscaling with common EOFs using sample EOFs from the esd-package

warrow <- function(col='grey',border='darkgrey',x=0.5,y=0.1,n=100) {
  par(mar=c(1,1,1,1),bty='n',xaxt='n',yaxt='n')
  plot(c(0,1.5),c(0,1),type='n',xlab='',ylab='')
  s <- 0.05*sin(seq(0,2*pi*3,length=n))
  x <- c(0,0.1,seq(0.1,0.9,length=n),0.9,1,
           0.9,seq(0.9,0.1,length=n),0.1,0)

  y <- c(0.5,1,rep(0.7,n)+s,1,0.5,
             0,rep(0.3,n)+s,0,0.5)
  polygon(x,y,col=col,border=border)
  #rect(0,0,1,1,lty=2)
}

frame <- function(x,ip,col='black') {
  data(geoborders)
  d <- dim(x)
  nx <- d[1]+2*d[3]; ny <- d[2]+2*d[3]
  X <- matrix(rep(NA,nx*ny),nx,ny)
  X[2*ip + (1:d[1]) -1, 2*ip + (1:d[2])-1] <- x[,,ip]
  #print(c(nx,ny,range(2*ip + (1:d[1]) -1),range(2*ip + (1:d[2])-1)))
  image(X,add=TRUE)
  xcoords <- approx(1:nx,seq(0,1,length=nx),xout=c(2*ip,2*ip+d[1]-1))$y
  ycoords <- approx(1:ny,seq(0,1,length=ny),xout=c(2*ip+1,2*ip+d[2]-1))$y
  rect(xcoords[1],ycoords[1],xcoords[2],ycoords[2],border=col)
  lines(geoborders,col='grey')
}

## Function that draws a wavy arrow
similar <- function(x) {
  polygon()
}

## Get some sample data for illustration purposes  
library(esd)
print('NCEP')
t2m.1 <- annual(subset(t2m.NCEP(),it=c(1961,1990),is=list(lon=c(-180,360),lat=c(-90,90))))
print('NorESM.M')
t2m.2 <- annual(subset(t2m.NorESM.M(),it=c(1961,1990),is=list(lon=c(-180,360),lat=c(-90,90))))
print('combine')
t2m <- combine(t2m.1,t2m.2)
print('EOF')
eof.t2m <- EOF(t2m,n=3)
x <- t(rbind(coredata(t2m),coredata(attr(t2m,'appendix.1'))))
dim(x) <- c(attr(t2m,'dimensions')[1:2],60)
## Dimensions of the spatial patterns
d <- dim(x)
print(d)
nx <- d[1]; ny <- d[2]; nt <- d[3]
col <- c(rep('black',30),rep('grey',30))
par(xaxt='n',yaxt='n',bty='n',xpd=NA)
X <- matrix(rep(NA,nx*ny),nx,ny)

## First figure: Illustrate the EOFs.
print('Fig. A1')
layout(matrix(c(c(1,1,0,2,2),c(1,1,0,2,2),c(1,1,3,4,4),c(1,1,0,4,4),c(1,1,0,5,0)),5,5,byrow = TRUE))
image(X,xlab='',ylab='',main=paste('Maps: ',d[1],'x',d[2],'x',d[3]))

for (it in seq(nt,1,by=-1)) frame(x,it,col=col[it])

image(X,xlab='',ylab='',main=paste('Modes: ',d[1],'x',d[2],'x 3 + (3 x',d[3],'+ 3)'))
z <- attr(eof.t2m,'pattern')
np <- dim(z)[3]
for (ip in seq(np,1,by=-1)) frame(z,ip)

warrow()

par(mar=c(1,1,1,1),xaxt='s',yaxt='s',bty='n',xpd=NA)
PCs <- zoo(rbind(coredata(eof.t2m),coredata(attr(eof.t2m,"appendix.1"))))
plot(PCs,plot.type='single',xlab='')
rect(1,-0.3,30.5,0.3,col=rgb(1,0.8,0.8),border=rgb(1,0.9,0.9))
rect(30.5,-0.3,60,0.3,col=rgb(0.8,0.8,1),border=rgb(0.9,0.9,1))
for (ip in 1:np) {
  lines(PCs[1:30,ip],col='black',lwd=2,lty=ip)
  lines(PCs[31:60,ip],col='grey',lwd=2,lty=ip)
}

plot(100*attr(eof.t2m,'eigenvalues')^2/attr(eof.t2m,'tot.var'),
     type='b',ylab='Variance (%)',xlab='Mode',ylim=c(0,100),lwd=2)

## Second figure: illustrate the downscaling
print('Fig. A2')
par(xaxt='n',yaxt='n',bty='n',xpd=NA)
layout(matrix(1:4,2,2))
## Get a time series that represents the predictand (local information)

nino3.4 <- NINO3.4()

par(xaxt='s',yaxt='s',bty='n',xpd=NA)
plot(PCs,plot.type='single',xlab='')
rect(1,-0.3,30.5,0.3,col=rgb(1,0.8,0.8),border=rgb(1,0.9,0.9))
rect(30.5,-0.3,60,0.3,col=rgb(0.8,0.8,1),border=rgb(0.9,0.9,1))
for (ip in 1:np) {
  lines(PCs[1:30,ip],col='black',lwd=2,lty=ip)
  lines(PCs[31:60,ip],col='grey',lwd=2,lty=ip)
}

plot(nino3.4)

plot(zoo(DS(nino3.4,eof.t2m)))


