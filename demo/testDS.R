library(esd)
lon <- c(-12,37)
lat <- c(52,72)
ylim <- c(-6,6)
t2m <- t2m.NCEP(lon=lon,lat=lat)
T2m <- t2m.NorESM.M(lon=lon,lat=lat)
data(Oslo)
X <- combine(t2m,T2m)

eof <- EOF(X,it=7)
ds <- DS(Oslo,eof)
plot(ds)

DS(Oslo,X,station=FALSE) -> y
Y <- combine.ds.comb(y)
plot(Y)


source("esd/R/DS.R",local=environment())
data(ferder)
t2m <- t2m.NCEP(lon=c(-30,50),lat=c(40,70))
slp <- slp.NCEP(lon=c(-30,50),lat=c(40,70))
T2m <- as.4seasons(t2m)
SLP <- as.4seasons(slp)
X <- EOF(T2m,it=1)
Z <- EOF(SLP,it=1)
y <- ferder
sametimescale(y,X) -> z
ym <- as.4seasons(y,FUN="mean")
ys <- as.4seasons(y,FUN="sd")
dsm <- DS(ym,X)
dss <- DS(ys,Z)

