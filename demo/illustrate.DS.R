## Illustrate downscaling with common EOFs using sample EOFs from the esd-package


frame <- function(x,ip,col='black') {
  d <- dim(x)
  nx <- d[1]+2*d[3]; ny <- d[2]+2*d[3]
  X <- matrix(rep(NA,nx*ny),nx,ny)
  X[2*ip + 1:d[1],2*ip + 1:d[2]] <- x[,,ip]
  image(X,add=TRUE)
  rect(2*ip,2*ip,2*ip+d[1],2*ip+d[2],border=col)
}

## Function that draws a wavy arrow
similar <- function(x) {
  polygon()
}

## Get some sample data for illustration purposes  
library(esd)
print('NCEP')
t2m.1 <- annual(subset(t2m.NCEP(),it=c(1991,2020),is=list(lon=c(-180,360),lat=c(-90,90))))
print('NorESM.M')
t2m.2 <- annual(subset(t2m.NorESM.M(),it=c(1991,2020),is=list(lon=c(-180,360),lat=c(-90,90))))
print('combine')
t2m <- combine(t2m.1,t2m.2)
print('EOF')
eof.t2m.NCEP <- EOF(t2m,n=5)
x <- attr(eof.t2m.NCEP,'pattern')
## Dimensions of the spatial patterns
d <- dim(x)
print(d)
nx <- d[1]; ny <- d[2]; np <- d[3]
col <- c(rep('black',60),rep('grey',60))
par(xaxt='n',yaxt='n',bty='n',xpd=NA)
X <- matrix(rep(NA,nx*ny),nx,ny)

## First figure: Illustrate the EOFs.
print('Fig. A1')
layout(matrix(c(rep(0,5),c(1,1,0,2,2),c(1,1,3,2,2),c(1,1,0,4,4),c(0,0,5,4,4)),5,5,byrow = TRUE))
image(X)
z <- t(as.matrix(t2m))
dim(z) <- attr(t2m,'dimensions')
nt <- dim(z)[3]
for (it in seq(nt,1,by=-1)) frame(z,it,col=col[it])

## Second figure: illustrate the downscaling
print('Fig. A2')
layout(matrix(1:4,2,2))
## Get spatial patterns 

image(X)
for (ip in seq(np,1,by=-1)) frame(x,ip)
