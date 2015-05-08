# Script for setting up and running CORDEX ESD experiment 1 for precip
# R.E. Benestad
# The experiment protocol: http://wcrp-cordex.ipsl.jussieu.fr/images/pdf/guidelines/CORDEX_ESD_Experiment1.pdf
# Tier 1: calibration: 1979-1995; validation: 1996-2006 
# Tier 2: 5-fold cross-validation: [1979,1983], [1984,1988],[1989,1993],[1994,1998],[1999,2003] 

rm(list=ls()); gc(reset=TRUE)
print('CORDEX.ESD.exp.1.pr')

# Number of PCAs used to represent the predictand data - determines the degree
# of detail but also the robustness of the results
# Trials:
# Multi-predictors:
#npca <- 5; lon <- c(-90,-30); lat <- c(-35,-15); eofs <- 1:10 # mu: r=0.786; fw: r=0.848
#npca <- 10; lon <- c(-90,-30); lat <- c(-35,-15); eofs <- 1:10 # mu: r=0.768; fw: r=0.831
#npca <- 10; lon <- c(-90,-30); lat <- c(-35,-15); eofs <- 1:5 # mu: r=0.811; fw: r=0.794
#npca <- 3; lon <- c(-90,-30); lat <- c(-35,-15); eofs <- 1:15 # mu: r=0.667; fw: r=0.824
#npca <- 15; lon <- c(-120,-20); lat <- c(-45,-5); eofs <- 1:15 # mu: r=0.588; fw: r=0.692
#npca <- 7; lon <- c(-70,-15); lat <- c(-37,-17); eofs <- 1:5 # mu: r=0.813; fw: r=0.799
#npca <- 3; lon <- c(-120,-20); lat <- c(-45,-5); eofs <- 1:15 # mu: r=663; fw: r=0.785
#npca <- 15; lon <- c(-70,-15); lat <- c(-37,-17); eofs <- 1:5 # mu: r=0.798; fw: r=0.79
## Test with a predictor region presumed to have no skill - correlation due to climatology:
#npca <- 10; lon <- c(50,120); lat <- c(35,65); eofs <- 1:15 # mu: r=687; fw: r=677
## ----------------------------
#npca <- 3; lon <- c(-70,-15); lat <- c(-37,-17); eofs <- 1:5  # mu: r= 0.83; fw: r= 0.8
# Single predictors:
#npca <- 3; lon <- c(-70,-15); lat <- c(-37,-17); eofs <- 1:5  #
npca <- 17; lon <- c(-70,-15); lat <- c(-37,-17); eofs <- 1:5  #

# load the predictands: CLARIS precip
print('Predictand')
load('claris.Pr.rda')
attr(Pr,'location')[77:81] <- c("Aerodromo de Pedro Juan Caballero","Aerodromo de Concepcion",
                                "Villarrica del Espedritu Santo","Aerodromo de Pilar",
                                "Encarnacion")
# Limit to the prescribed interval
Pr <- subset(Pr,it=c(1979,2006))

# retrieve the predictors
# Out-going long-wave radiation
print('predictors')
olr <- annual(retrieve('data/ERAINT/ERAINT_olr_mon.nc',
                       lon=lon,lat=lat),FUN='mean')
attr(olr,'unit') <- "W * m**-2"
eof.olr <- EOF(olr)

# Temperature
t2m <- retrieve('data/ERAINT/ERAINT_t2m_mon.nc',lon=lon,lat=lat)
es <- annual(C.C.eq(t2m),FUN='mean')
eof.es <- EOF(es)

# Mean sea level pressure
slp <- annual(retrieve('data/ERAINT/ERAINT_slp_mon.nc',
                       lon=lon,lat=lat),FUN='mean')
eof.slp <- EOF(slp)

t2m.c1 <- subset(t2m,it=c(1979,1995));  t2m.v1 <- subset(t2m,it=c(1996,2006));
slp.c1 <- subset(t2m,it=c(1979,1995));  slp.v1 <- subset(t2m,it=c(1996,2006));
olr.c1 <- subset(t2m,it=c(1979,1995));  olr.v1 <- subset(t2m,it=c(1996,2006));

# Process the precipitation - predictand:
print('prepare the predictand')
print('estimate annual wet-day mean')

## Extract the predictands for the calibration period:
mu <- annual(Pr,FUN='wetmean',nmin=100)
mu1 <- subset(mu,it=c(1979,1995))

print('estimate annual wet-day frequency')
fw <- annual(Pr,FUN='wetfreq',nmin=100)
fw1 <- subset(fw,it=c(1979,1995))

## Missing values are a problem for PCA - replace missing anomalies with zero

miss.m <- !is.finite(coredata(mu))
coredata(mu)[miss.m] <- 0
miss.m1 <- !is.finite(coredata(mu1))
coredata(mu1)[miss.m1] <- 0
miss.f <- !is.finite(coredata(fw))
coredata(fw)[miss.f] <- 0
miss.f1 <- !is.finite(coredata(fw1))
coredata(fw1)[miss.f1] <- 0

pca.mu <- PCA(mu,neofs=npca)
pca.mu1 <- PCA(mu1,neofs=npca)
pca.fw <- PCA(fw,neofs=npca)
pca.fw1 <- PCA(fw1,neofs=npca)

print('Tier 1')
print('the actual downsscaling')
#z.mu1 <- DS(pca.mu1,list(olr=eof.olr,es=eof.es),
#            rmtrend=FALSE,m=NULL,eofs=eofs,verbose=FALSE)
z.mu1 <- DS(pca.mu1,eof.es,rmtrend=FALSE,m=NULL,eofs=eofs,verbose=FALSE)
#z.fw1 <- DS(pca.fw1,list(olr=eof.olr,es=eof.es,slp=eof.slp),
#            rmtrend=FALSE,m=NULL,eofs=eofs,verbose=FALSE)
z.fw1 <- DS(pca.fw1,eof.slp,rmtrend=FALSE,m=NULL,eofs=eofs,verbose=FALSE)
#plot(z.mu1)

## Predict the results.
print('predict the results')
#mu.tier1 <-as.station(predict(z.mu1,newdata=list(olr=eof.olr,es=eof.es),verbose=FALSE))
mu.tier1 <-as.station(predict(z.mu1,newdata=eof.es,verbose=FALSE))
#fw.tier1 <-as.station(predict(z.fw1,newdata=list(olr=eof.olr,es=eof.es,slp=eof.slp)))
fw.tier1 <-as.station(predict(z.fw1,newdata=eof.slp,verbose=FALSE))

print('Tier 2')
print('the actual downsscaling')
# All post-processing is finished. Proceed with the actual downscaling:
print('wet-day mean')
#z.mu <- DS(pca.mu,list(olr=eof.olr,es=eof.es),
#           m='cordex-esd-exp1',eofs=eofs,rmtrend=FALSE,verbose=FALSE)
z.mu2 <- DS(pca.mu,eof.es,m='cordex-esd-exp1',eofs=eofs,rmtrend=FALSE,verbose=FALSE)
print('wet-day frequency')
#z.fw <- DS(pca.fw,list(olr=eof.olr,es=eof.es,slp=eof.slp),
#           m='cordex-esd-exp1',eofs=eofs,rmtrend=FALSE,verbose=FALSE)
z.fw2 <- DS(pca.fw,eof.slp,m='cordex-esd-exp1',eofs=eofs,rmtrend=FALSE,verbose=FALSE)

# Post-processing - get the data into the right format:
print('post-processing of the results')
#mu.ds <- attr(z.mu,'evaluation')
#fw.ds <- attr(z.fw,'evaluation')
mu.tier2 <- pca2station(z.mu2,what='xval')
fw.tier2 <- pca2station(z.fw2,what='xval')

## Reset missing values to missing values
coredata(mu)[miss.m] <- NA
coredata(fw)[miss.f] <- NA
coredata(mu1)[miss.m1] <- NA
coredata(fw1)[miss.f1] <- NA

y <- pca2station(z.mu2)
x <- matchdate(mu,y)

# Check the results:
print('check the results')
dev.new(width=5,height=9)
par(bty='n',las=1,oma=rep(0.25,4),mfcol=c(2,1),cex=0.5)
plot(coredata(x),coredata(y),
     pch=19,col=rgb(0.5,0.5,1,0.25),
     xlab=expression(paste('Observed ',mu,(mm/day))),
     ylab=expression(paste('Downscaled ',mu,(mm/day))),
     main='Wet-day mean precipitation',
     sub=paste('predictand: CLARIS; #PCA=',npca,' #EOFs=',length(eofs)))
grid()
x <- c(coredata(matchdate(mu,mu.tier2)))
y <- c(coredata(mu.tier2))
points(x,y,col=rgb(0,0,1,0.5),lwd=2)
abline(lm(y ~ x),col=rgb(0,0,1),lty=2)
ok <-  is.finite(x) & is.finite(y)
mtext(side=4,paste('r=',round(cor(x[ok],y[ok]),3)),las=3)

# Wet-day frequency:

y <- pca2station(z.fw2)
x <- matchdate(fw,y)

# Check: Figure: scatter plot
plot(coredata(x),coredata(y),
     pch=19,col=rgb(0.4,0.4,0.8,0.25),
     xlab=expression(paste('Observed ',f[w])),
     ylab=expression(paste('Downscaled ',f[w])),
     main='Wet-day frequency',
     sub=paste('predictand: CLARIS; #PCA=',npca,' #EOFs=',length(eofs)))
grid()
x <- c(coredata(matchdate(fw,fw.tier2)))
y <- c(coredata(fw.tier2))
points(x,y,col=rgb(0,0,1,0.5),lwd=2)
abline(lm(y ~ x),col=rgb(0,0,1),lty=2)
ok <-  is.finite(x) & is.finite(y)
mtext(side=4,paste('r=',round(cor(x[ok],y[ok]),3)),las=3)

dev.copy2pdf(file='CORDEX.ESD.exp1.pr.esd.pdf')

mu.tier1 <- subset(mu.tier1,it=c(1996,2006))
fw.tier1 <- subset(fw.tier1,it=c(1996,2006))
mu.tier2 <- subset(mu.tier2,it=c(1979,2003))
fw.tier2 <- subset(fw.tier2,it=c(1979,2003))

X <- list(mu.tier1=mu.tier1,fw.tier1=fw.tier1,
          mu.tier2=mu.tier2,fw.tier2=fw.tier2,
          z.mu=z.mu2,z.fw=z.fw2,mu0=mu,fw0=fw)
attr(X,'description') <- 'cross-validation'
attr(X,'experiment') <- 'CORDEX ESD experiment 1'
attr(X,'predictand_file') <- 'claris.Pr.rda'
attr(X,'predictor_file') <- c('ERAINT_olr_mon.nc','ERAINT_t2m_mon.nc','ERAINT_slp_mon.nc')
attr(X,'predictor_domain') <- 'lon=c(-90,-30),lat=c(-35,-15)'
attr(X,'history') <- history.stamp()
attr(X,'R-script') <- readLines('CORDEX.ESD.exp.1.pr.R')

save(file='CORDEX.ESD.exp1.pr.esd.rda',X)

print('finished')
