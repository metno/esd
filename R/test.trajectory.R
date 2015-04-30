# K. Parding
# Test functions for handling IMILAST stormtrack files.

library(esd)

## source('~/git/esd/R/as.trajectory.R')
## source('~/git/esd/R/subset.trajectory.R')
## source('~/git/esd/R/annual.R')
## source('~/git/esd/R/rtools.trajectory.R')
## source('~/git/esd/R/map.trajectory.R')
## source('~/git/esd/R/plot.trajectory.R')
## source('~/git/esd/R/pca.trajectory.R')
data(imilast.M03)
y <- imilast.M03

fname="/vol/fou/klima/IMILAST/ERAinterim_1.5_NH_M03_19890101_20090331_ST.txt"
m03 <- read.fwf(fname,width=c(2,7,4,11,5,5,5,4,7,7,10),
        col.names=c("Code99","Trajectory","timeStep","Date","Year",
        "Month","Day","Time","Lon","Lat","slp"))
x <- as.trajectory(m03,src="IMILAST",method="M03",
      param="storm tracks",longname="mid-latitude storm tracks",
      url="http://journals.ametsoc.org/doi/abs/10.1175/BAMS-D-11-00154.1",
      reference="Neu, et al. , 2013: IMILAST: A Community Effort to Intercompare Extratropical Cyclone Detection and Tracking Algorithms. Bull. Amer. Meteor. Soc., 94, 529–547.",n=12,verbose=T)
it <- param.trajectory(x,param='slp',FUN=min)<980 &
      param.trajectory(x,param='n')>8
is <- list(lon=c(-60,60),lat=c(40,90))
y <- subset.trajectory(x,it=it,is=is,verbose=FALSE)
#imilast.M03 <- y
#save(imilast.M03,file="imilast.M03.rda")
#load("imilast.M03.rda")

m03.n <- tdf2numeric(m03)
m03.r <- matrix2tdf(x)

# Example of as.trajectory:
id = c(rep(1,9),rep(2,5),rep(3,12))
lon = c(1:9,8:12,12:23)
lat = c(seq(45,49,0.5)+sample(seq(0,1,0.1),9),
  seq(47,50,0.7)+sample(seq(0,1,0.1),5),
  seq(50,60,0.9)+sample(seq(0,1,0.1),12))
year = rep(2009,9+5+12)
month = rep(12,9+5+12)
day = c(rep(12,9),rep(14,5),rep(24,12))
time = c(2:10,7:11,11:22)
x <- data.frame("Trajectory"=id,"Lat"=lat,"Lon"=lon,
                "Year"=year,"Month"=month,"Day"=day,"Time"=time)

st <- season.trajectory(y)
slp <- param.trajectory(y,param='slp',FUN=min)
lat <- param.trajectory(y,param='lat',FUN=max)
lon <- param.trajectory(y,param='lon',FUN=mean)
n <- param.trajectory(y,param='n')
start <- param.trajectory(y,param='start')

y.sorted <- sort(y)
y.anomaly <- anomaly(y,param=c('lon','lat'),type='first',verbose=T)
y.anomaly <- anomaly(y.anomaly,param=c('slp'),type='mean',verbose=T)
yr <- anomaly2trajectory(y.anomaly,verbose=T)

y.pfit <- pfit.trajectory(y.anomaly)
coefs <- attr(y.pfit,'coefs')
map(y.pfit,alpha=0.1)
for(i in which(coefs[2,]>q95(coefs[2,]))) {
  print(i)
  print(coefs[,i])
  lines(y.pfit[i,colnames(y)=='lon'],y.pfit[i,colnames(y)=='lat'],lwd=1,col='black')
}
for(i in which(coefs[2,]<q5(coefs[2,]))) {
  print(i)
  print(coefs[,i])
  lines(y.pfit[i,colnames(y)=='lon'],y.pfit[i,colnames(y)=='lat'],lwd=1,col='green')
}

plot.trajectory(y)
map(y,projection='latlon',xlim=c(-90,90))
map(y,projection='sphere')
map.hexbin.trajectory(y,xlim=c(-90,90))
map.sunflower.trajectory(y,xlim=c(-90,90))

# PCA
y.anomaly <- anomaly.trajectory(y,param=c('lon','lat'))
yp <- PCA.trajectory(y,neofs=20,anomaly=T,param=c('lat','lon'))
plot.pca.storm(yp)

# Compare histogram to Poisson distribution

# by month
ylim <- c(0,0.15)
x <- seq(0,50)
dx <- 1
n <- count.storm(y,by='month')
h <- hist(n,xlim=c(min(x),max(x)),breaks=seq(min(x),max(x),dx),ylim=ylim,
     probability=TRUE,col="gray",border="white",main="Storm count")
p <- dpois(x,lambda=mean(n))
lines(x,p,col='red',lty=1)
legend('topleft',legend=c('histogram','Poisson distribution'),
       col=c("gray","red"),lty=c(1,1),lwd=c(7,1))

# by year
ylim <- c(0,0.05)
x <- seq(40,150)
dx <- 5
n <- count.storm(y,by='year')
h <- hist(n,xlim=c(min(x),max(x)),breaks=seq(min(x),max(x),dx),ylim=ylim,
     probability=TRUE,col="gray",border="white",main="Storm count")
p <- dpois(x,lambda=mean(n))
lines(x,p,col='red',lty=1)
legend('topleft',legend=c('histogram','Poisson distribution'),
       col=c("gray","red"),lty=c(1,1),lwd=c(7,1))



# reconstruct stormmatrix
yr <- pca2storm(yp)
map.storm(y,projection='sphere',alpha=0.1)
map.storm(yr,projection='sphere',col='blue',alpha=0.1)

# map.pca.storm
map.pca.storm(yp)
map.pca.storm(yp,projection='latlon',xlim=c(-40,120))


# Plot storm tracks
# png("/home/kajsamp/Desktop/stormtracks_djf_M07.png",
#     height=800,width=800,res=120,pointsize=14)
#map.stormmatrix(x,ylim=c(30,90),xlim=c(-60,60),
#    it='djf',col='red',alpha=0.1,projection='latlon',
#    pfit=F,new=T)
# dev.off()

# Color = year (should add scale)
#map.stormmatrix(x,ylim=c(30,90),xlim=c(-60,60),
#    it='djf',col='red',alpha=0.5,projection='latlon',
#    FUN=year.stormmatrix,pfit=F,new=T)

## # Show patterns in units of lon-lat. 
## U <- attr(x2.pca,'pattern')
## V <- coredata(x2.pca)
## W <- attr(x2.pca,'eigenvalues')

## x2.r <- pca2stormmatrix(x2.pca) # reconstruct x2 from PCA

## PC1 <- matrix(NA*dim(V)[1]*dim(U)[2],dim(V)[1],dim(U)[1])
## PC2 <- PC1
## PC3 <- PC1
## for (i in 1:dim(V)[1]) {
##   PC1[i,] <- V[i,1]*(U[,1]*W[1])
##   PC2[i,] <- V[i,2]*(U[,2]*W[2])
##   PC3[i,] <- V[i,3]*(U[,3]*W[3])   
## }
#colnames(PC1) <- colnames(x2.r)
#colnames(PC2) <- colnames(x2.r)
#colnames(PC3) <- colnames(x2.r)
#map.stormmatrix(PC1,projection='latlon',col=adjustcolor('red',alpha.f=0.2))
#map.stormmatrix(PC2,projection='latlon',col=adjustcolor('blue',alpha.f=0.2))
#map.stormmatrix(PC3,projection='latlon',col=adjustcolor('green',alpha.f=0.2))

## dev.new()
## plot(t(PC1[,1:10]),t(PC1[,11:20]),type='n')
## matlines(t(PC1[,1:10]),t(PC1[,11:20]),type='l',lty=1,
##          col=adjustcolor('red',alpha.f=0.2))
## points(PC1[,1],PC1[,11],pch=19,
##          col=adjustcolor('black',alpha.f=0.2))

## dev.new()
## plot(t(PC2[,1:10]),t(PC2[,11:20]),type='n')
## matlines(t(PC2[,1:10]),t(PC2[,11:20]),type='l',lty=1,
##          col=adjustcolor('blue',alpha.f=0.2))
## points(PC2[,1],PC2[,11],pch=19,
##          col=adjustcolor('black',alpha.f=0.2))


## dev.new()
## plot(t(PC3[,1:10]),t(PC3[,11:20]),type='n')
## matlines(t(PC3[,1:10]),t(PC3[,11:20]),type='l',lty=1,
##          col=adjustcolor('green',alpha.f=0.2))
## points(PC3[,1],PC3[,11],pch=19,
##          col=adjustcolor('black',alpha.f=0.2))

## PC123 <- PC1+PC2+PC3
## dev.new()
## plot(t(PC123[,1:10]),t(PC123[,11:20]),type='n')
## matlines(t(PC123[,1:10]),t(PC123[,11:20]),type='l',lty=1,
##          col=adjustcolor('grey',alpha.f=0.2))
## points(PC123[,1],PC123[,11],pch=19,
##          col=adjustcolor('black',alpha.f=0.2))



#dev.new()
#plot(0,0,type='n',xlim=c(-90,90),ylim=c(-90,90))
#for (i in 1:m) {
#  PCi <- median(V[,i])*(U[,i]*W[i])
#  lines(PCi[1:10],PCi[11:20],lty=1,col=colvec[i])
#  points(PCi[1],PCi[11],col=colvec[i],pch=19)
#  PCi <- max(V[,i])*(U[,i]*W[i])
#  lines(PCi[1:10],PCi[11:20],lty=3,col=colvec[i])
#  points(PCi[1],PCi[11],col=colvec[i],pch=20)
#  PCi <- min(V[,i])*(U[,i]*W[i])
#  lines(PCi[1:10],PCi[11:20],lty=4,col=colvec[i])
#  points(PCi[1],PCi[11],col=colvec[i],pch=20)
#}




## png("/home/kajsamp/Desktop/stormtracks_djf_M07.png",
##     height=800,width=800,res=120,pointsize=14)
## map.stormmatrix(x,ylim=c(30,90),xlim=c(-60,60),
##    it='djf',col='red',alpha=0.1,projection='latlon',
##    pfit=F,new=F)
## dev.off()

## png("/home/kajsamp/Desktop/stormtracks_djf_M07_sphere.png",
##     height=800,width=800,res=120,pointsize=12,type='cairo')
## map.stormmatrix(x,it='djf',col='red',alpha=0.1,
##     projection='sphere',pfit=F,new=F)
## dev.off()

## png("/home/kajsamp/Desktop/stormtracks_djf_M07_years.png",
##     height=800,width=800,res=120,pointsize=14)
## map.stormmatrix(x,ylim=c(30,90),xlim=c(-60,60),
##    it='dec',FUN=year.stormmatrix,col='red',alpha=0.5,projection='latlon',
##    pfit=F,new=F)
## dev.off()

## png("/home/kajsamp/Desktop/stormtracks_djf_M07_sunflower50.png",
##     height=800,width=800,res=120,pointsize=14)
## sunflower.stormmatrix(x,it='djf',col='blue',alpha=0.7,
##     xlim=c(-60,60),ylim=c(30,90),projection='latlon',
##     petalsize=50,dlon=4,dlat=2)
## dev.off()

## png("/home/kajsamp/Desktop/stormtracks_djf_M07_sunflower25.png",
##     height=800,width=800,res=120,pointsize=14)
## sunflower.stormmatrix(x,it='djf',col='blue',alpha=0.7,
##     xlim=c(-60,60),ylim=c(30,90),projection='latlon',
##     petalsize=25,dlon=4,dlat=2)
## dev.off()

# test subset
y <- subset.storm(x,it=c("1988-01-01","2000-12-31"))
y <- subset.storm(x,it='mam')
y <- subset.storm(x,it=c(1988,2000),is=list(lat=c(45,90),lon=c(-40,60)))

djf <- subset.storm(x,it='djf')
mam <- subset.storm(x,it='mam')
jja <- subset.storm(x,it='jja')
son <- subset.storm(x,it='son')

arctic <- subset.storm(djf,is=list(lat=c(75,90),lon=c(-30,30)))
north <- subset.storm(djf,is=list(lat=c(60,75),lon=c(-30,30)))
central <- subset.storm(djf,is=list(lat=c(45,60),lon=c(-30,30)))
south <- subset.storm(djf,is=list(lat=c(30,45),lon=c(-30,30)))

# test map
map.storm(x)
map.storm(djf)
map.storm(arctic)
map.sunflower.storm(djf)

# test stormtools
# stormcount
n <- count.storm(x,by='year')
n1 <- count.storm(djf)
n2 <- count.storm(mam)
n3 <- count.storm(jja)
n4 <- count.storm(son)
nmin <- min(sapply(list(n1,n2,n3,n4),min)) 
nmax <- max(sapply(list(n1,n2,n3,n4),max))
plot(n1,type='lines',col='blue',ylab='Storm frequency',
     ylim=c(signif(nmin,2)-10,signif(nmax,2)+10))
lines(n2,lty=2,col='green')
lines(n3,lty=3,col='red')
lines(n4,lty=4,col='purple')

n1 <- count.storm(arctic)
n2 <- count.storm(north)
n3 <- count.storm(central)
n4 <- count.storm(south)
nmin <- min(sapply(list(n1,n2,n3,n4),min)) 
nmax <- max(sapply(list(n1,n2,n3,n4),max))
plot(n1,type='lines',col='blue',ylab='Storm frequency',
     ylim=c(signif(nmin,2)-10,signif(nmax,2)+10))
lines(n2,lty=2,col='green')
lines(n3,lty=3,col='red')
lines(n4,lty=4,col='purple')

r1 <- subset.storm(djf,is=list(lat=c(30,60),lon=c(-60,-10)))
r2 <- subset.storm(djf,is=list(lat=c(60,90),lon=c(-60,-10)))
r3 <- subset.storm(djf,is=list(lat=c(30,60),lon=c(-10,40)))
r4 <- subset.storm(djf,is=list(lat=c(60,90),lon=c(-10,40)))
n1 <- count.storm(r1)
n2 <- count.storm(r2)
n3 <- count.storm(r3)
n4 <- count.storm(r4)
nmin <- min(sapply(list(n1,n2,n3,n4),min)) 
nmax <- max(sapply(list(n1,n2,n3,n4),max))
plot(n1,type='lines',col='red',ylab='Storm frequency',
     ylim=c(signif(nmin,2)-10,signif(nmax,2)+10))
lines(n2,lty=1,col='blue')
lines(n3,lty=2,col='red')
lines(n4,lty=2,col='blue')

# polyfit
#p <- polyfit.stormmatrix(x)
#y <- x
#y[,11:20] <- t(p)
#fn <- function(x) (x[33]>5)
#z <- subset.stormmatrix(y,is=list(FUN=fn))
#map.stormmatrix(z)
