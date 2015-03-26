# Script for setting up and running CORDEX ESD experiment 1 for precip
# R.E. Benestad

print('CORDEX.ESD.exp.1.pr')

# Number of PCAs used to represent the predictand data - determines the degree
# of detail but also the robustness of the results
npca <- 20

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
                       lon=c(-90,-30),lat=c(-35,-15)),FUN='mean')
attr(olr,'unit') <- "W * m**-2"
eof.olr <- EOF(olr)

# Temperature
t2m <- retrieve('data/ERAINT/ERAINT_t2m_mon.nc',lon=c(-90,-30),lat=c(-35,-15))
es <- annual(C.C.eq(t2m),FUN='mean')
eof.es <- EOF(es)

# Mean sea level pressure
slp <- annual(retrieve('data/ERAINT/ERAINT_slp_mon.nc',
                       lon=c(-90,-30),lat=c(-35,-15)),FUN='mean')
eof.slp <- EOF(slp)

t2m.c1 <- subset(t2m,it=c(1979,1995));  t2m.v1 <- subset(t2m,it=c(1996,2006));
slp.c1 <- subset(t2m,it=c(1979,1995));  slp.v1 <- subset(t2m,it=c(1996,2006));
olr.c1 <- subset(t2m,it=c(1979,1995));  olr.v1 <- subset(t2m,it=c(1996,2006));


# Process the precipitation - predictand:
print('prepare the predictand')
print('estimate annual wet-day mean')
mu <- annual(Pr,FUN='wetmean',nmin=100)
pca.mu <- PCA(mu)
print('estimate annual wet-day frequency')
fw <- annual(Pr,FUN='wetfreq',nmin=100)
pca.fw <- PCA(fw)

pca.mu <- subset(pca.mu,pattern=1:npca)
pca.fw <- subset(pca.fw,pattern=1:npca)

print('the actual downsscaling')
# All post-processing is finished. Proceed with the actual downscaling:
print('wet-day mean')
z.mu <- DS(pca.mu,list(olr=eof.olr,es=eof.es),
           m='cordex-esd-exp1',eofs=1:10,detrend=FALSE,verbose=FALSE)
print('wet-day frequency')
z.fw <- DS(pca.fw,list(olr=eof.olr,es=eof.es,slp=eof.slp),
           m='cordex-esd-exp1',eofs=1:10,detrend=FALSE,verbose=FALSE)

# Post-processing - get the data into the right format:
print('post-processing of the results')
#mu.ds <- attr(z.mu,'evaluation')
#fw.ds <- attr(z.fw,'evaluation')
mu.ds <- pca2station(z.mu,what='xval')
fw.ds <- pca2station(z.fw,what='xval')

y <- pca2station(z.mu)
x <- matchdate(mu,y)

# Check the results:
print('check the results')
dev.new(width=5,height=9)
par(bty='n',las=1,oma=rep(0.25,4),mfcol=c(2,1),cex=0.5)
plot(coredata(anomaly(x)),coredata(anomaly(y)),
     pch=19,col=rgb(0.5,0.5,1,0.5),
     xlab=expression(paste('Observed ',mu,(mm/day))),
     ylab=expression(paste('Downscaled ',mu,(mm/day))),
     main='Wet-day mean precipitation',
     sub=paste('predictand: CLARIS; #PCA=',npca))
grid()
x <- c(coredata(anomaly(matchdate(mu,mu.ds))))
y <- c(coredata(anomaly(mu.ds)))
points(x,y,col=rgb(0,0,1,0.5),lwd=2)
abline(lm(y ~ x),col=rgb(0,0,1),lty=2)
ok <-  is.finite(x) & is.finite(y)
mtext(side=4,paste('r=',round(cor(x[ok],y[ok]),3)),las=3)

# Wet-day frequency:

y <- pca2station(z.fw)
x <- matchdate(fw,y)

# Check: Figure: scatter plot
plot(coredata(anomaly(x)),coredata(anomaly(y)),
     pch=19,col=rgb(0.4,0.4,0.8,0.5),
     xlab=expression(paste('Observed ',f[w])),
     ylab=expression(paste('Downscaled ',f[w])),
     main='Wet-day frequency',
     sub=paste('predictand: CLARIS; #PCA=',npca))
grid()
x <- c(coredata(anomaly(matchdate(fw,fw.ds))))
y <- c(coredata(anomaly(fw.ds)))
points(x,y,col=rgb(0,0,1,0.5),lwd=2)
abline(lm(y ~ x),col=rgb(0,0,1),lty=2)
ok <-  is.finite(x) & is.finite(y)
mtext(side=4,paste('r=',round(cor(x[ok],y[ok]),3)),las=3)



dev.copy2eps(file='CORDEX.ESD.exp1.pr.esd.eps')

X <- list(mu.exp1=mu.ds,fw.exp1=fw.ds,z.mu=z.mu,z.fw=z.fw,mu0=mu,fw0=fw)
attr(X,'description') <- 'cross-validation'
attr(X,'experiment') <- 'CORDEX ESD experiment 1'
attr(X,'predictand_file') <- 'claris.Pr.rda'
attr(X,'predictor_file') <- c('ERAINT_olr_mon.nc','ERAINT_t2m_mon.nc','ERAINT_slp_mon.nc')
attr(X,'predictor_domain') <- 'lon=c(-90,-30),lat=c(-35,-15)'
attr(X,'history') <- history.stamp()
attr(X,'R-script') <- readLines('CORDEX.ESD.exp.1.pr.R')

save(file='CORDEX.ESD.exp1.pr.esd.rda',X)

