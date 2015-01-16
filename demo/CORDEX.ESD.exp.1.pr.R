# Script for setting up and running CORDEX ESD experiment 1
# 

# load the predictands: CLARIS precip
load('~/Dropbox/Public/CORDEX-ESDM-data-clumps/claris.Pr.rda')

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
mu <- annual(Pr,FUN='wetmean')
pca.mu <- PCA(mu)
fw <- annual(Pr,FUN='wetfreq')
pca.fw <- PCA(fw)

# The experiment results are provided by the cross-validation 
#z.mu <- DS(pca.mu,list(olr=eof.olr,es=eof.es,slp=eof.slp),
#           eofs=1:20,m='cordex-esd-exp1',detrend=FALSE,verbose=TRUE)
#mu.ds <- pca2station(z.mu)
#z.fw <- DS(pca.fw,list(olr=eof.olr,es=eof.es,slp=eof.slp),
#           eofs=1:20,m='cordex-esd-exp1',detrend=FALSE,verbose=TRUE)


for (i in 1:dim(mu)[2]) {
  z.mu <- DS(subset(mu,is=i),list(olr=eof.olr,es=eof.es),
             m='cordex-esd-exp1',eofs=1:10,detrend=FALSE)
  z.fw <- DS(subset(fw,is=1),list(olr=eof.olr,es=eof.es,slp=eof.slp),
             m='cordex-esd-exp1',eofs=1:10,detrend=FALSE)
  if (i == 1) {
    mu.ds <- attr(z.mu,'evaluation')[,2]
    fw.ds <- attr(z.fw,'evaluation')[,2]
  } else {
    mu.ds <- merge(mu.ds,attr(z.mu,'evaluation')[,2])
    fw.ds <- merge(fw.ds,attr(z.fw,'evaluation')[,2])
  }
}
mu.ds <- attrcp(mu,mu.ds)
attr(mu.ds,'description') <- 'cross-validation'
attr(mu.ds,'experiment') <- 'CORDEX ESD experiment 1'
class(mu.ds) <- class(mu)
fw.ds <- attrcp(fw,fw.ds)
class(fw.ds) <- class(fw)
attr(fw.ds,'description') <- 'cross-validation'
attr(fw.ds,'experiment') <- 'CORDEX ESD experiment 1'


