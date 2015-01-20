# Script for setting up and running CORDEX ESD experiment 1
# 

# load the predictands: CLARIS precip
load('~/Dropbox/Public/CORDEX-ESDM/CORDEX-ESDM-data-clumps/claris.Pr.rda')
attr(Pr,'location')[77:81] <- c("Aerodromo de Pedro Juan Caballero","Aerodromo de Concepcion",
                                "Villarrica del Espedritu Santo","Aerodromo de Pilar",
                                "Encarnacion")

# retrieve the predictors
# Out-going long-wave radiation
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

# Process the precipitation - predictand:
mu <- annual(Pr,FUN='wetmean',nmin=100)
pca.mu <- PCA(mu)
fw <- annual(Pr,FUN='wetfreq',nmin=100)
pca.fw <- PCA(fw)

# The experiment results are provided by the cross-validation 
#z.mu <- DS(pca.mu,list(olr=eof.olr,es=eof.es,slp=eof.slp),
#           eofs=1:20,m='cordex-esd-exp1',detrend=FALSE,verbose=TRUE)
#mu.ds <- pca2station(z.mu)
#z.fw <- DS(pca.fw,list(olr=eof.olr,es=eof.es,slp=eof.slp),
#           eofs=1:20,m='cordex-esd-exp1',detrend=FALSE,verbose=TRUE)



z.mu <- DS(pca.mu,list(olr=eof.olr,es=eof.es),
           m='cordex-esd-exp1',eofs=1:10,detrend=FALSE)
z.fw <- DS(pca.fw,list(olr=eof.olr,es=eof.es,slp=eof.slp),
           m='cordex-esd-exp1',eofs=1:10,detrend=FALSE)
mu.ds <- attr(z.mu,'evaluation')
fw.ds <- attr(z.fw,'evaluation')

y <- pca2station(z.mu)
x <- matchdate(mu,y)

dev.new(width=5,height=9)
par(bty='n',las=1,oma=rep(0.25,4),mfcol=c(2,1),cex=0.5)
plot(coredata(x),coredata(y),
     pch=19,col=rgb(0,0,1,0.5),
     xlab=expression(paste('Observed ',mu,(mm/day))),
     ylab=expression(paste('Downscaled ',mu,(mm/day))),
     main='Wet-day mean precipitation',
     sub=paste('predictand: CLARIS; #PCA=',npca))
grid()

y <- pca2station(z.fw)
x <- matchdate(fw,y)

                                        # Check: Figure: scatter plot
x <- matchdate(subset(st4s,it=season),txm.ds)
plot(coredata(fw.ds[,1]),coredata(fw.ds[,2]),
     pch=19,col=rgb(0,0,1,0.5),
     xlab=expression(paste('Observed ',f[w])),
     ylab=expression(paste('Downscaled ',f[w])),
     main='Wet-day frequency',
     sub=paste('predictand: CLARIS; #PCA=',npca))
grid()

mu.ds <- attrcp(mu,mu.ds)
attr(mu.ds,'description') <- 'cross-validation'
attr(mu.ds,'experiment') <- 'CORDEX ESD experiment 1'
class(mu.ds) <- class(mu)
fw.ds <- attrcp(fw,fw.ds)
class(fw.ds) <- class(fw)
attr(fw.ds,'description') <- 'cross-validation'
attr(fw.ds,'experiment') <- 'CORDEX ESD experiment 1'
attr(X,'predictand_file') <- 'claris.Pr.rda'
attr(X,'predictor_file') <- c('ERAINT_olr_mon.nc','ERAINT_t2m_mon.nc','ERAINT_slp_mon.nc')
attr(X,'predictor_domain') <- 'lon=c(-90,-30),lat=c(-35,-15)'
attr(X,'history') <- history.stamp()
attr(X,'R-script') <- readLines('CORDEX.ESD.exp.1.pr.R')


