# Script for setting up and running CORDEX ESD experiment 1
# 
npca <- 4

# load the predictands: CLARIS precip
print('Get the CLARIS data')
load('~/Dropbox/Public/CORDEX-ESDM/CORDEX-ESDM-data-clumps/claris.Tx.rda')

# Process the precipitation - predictand as annual mean and annual standard deviation:
# Perhaps change to use seasonal rather than annual?
print('Estimate seasonal statistics')
mt4s <- as.4seasons(Tx,FUN='mean')
st4s <- as.4seasons(anomaly(Tx),FUN='sd')

# retrieve the predictors
print('Get the predictor data')
# Out-going long-wave radiation
olr <- as.4seasons(retrieve('data/ERAINT/ERAINT_olr_mon.nc',
                       lon=c(-90,-30),lat=c(-35,-15)),FUN='mean')
attr(olr,'unit') <- "W * m**-2"

# Temperature
t2m <- as.4seasons(retrieve('data/ERAINT/ERAINT_t2m_mon.nc',
                            lon=c(-90,-30),lat=c(-35,-15)))

# Mean sea level pressure
slp <- as.4seasons(retrieve('data/ERAINT/ERAINT_slp_mon.nc',
                       lon=c(-90,-30),lat=c(-35,-15)),FUN='mean')

# Loop through the seasons:
for (season in c('djf','mam','jja','son')) {
  print(paste('Season:',season))
  
  eof.slp <- EOF(subset(slp,it=season))
  eof.olr <- EOF(subset(olr,it=season))
  eof.t2m <- EOF(subset(t2m,it=season))
# Combine the local predictand information in the form of PCAs
  pca.mt <- PCA(subset(mt4s,it=season))
  pca.st <- PCA(subset(st4s,it=season))

  pca.mt <- subset(pca.mt,pattern=1:npca)
  pca.st <- subset(pca.st,pattern=1:npca)

# The experiment results are provided by the cross-validation deined by parameter 'm' (used in crossval).
# Impoertan to set detrend=FALSE to retain original data in the cross-validation.
# First downscale the mean values:
  z.mt <- DS(pca.mt,eof.t2m,
               m='cordex-esd-exp1',detrend=FALSE,verbose=TRUE)
# The results are in the form of a PCA-object - convert back to a group of stations            
  txm.ds <- pca2station(z.mt)

# Repeat the downscaling for the standard deviation:
  z.st <- DS(pca.st,list(t2m=eof.t2m,slp=eof.slp,olr=eof.olr),
               m='cordex-esd-exp1',detrend=FALSE,verbose=TRUE)
  txs.ds <- pca2station(z.st)

# Extract the predicted cross-validation results which follow the experiment:
  d <- dim(attr(tx.ds,'evaluation'))

# only grab the series of predicted values - not the original data used for calibration
  exp1.tx <- attr(tx.ds,'evaluation')[,seq(2,d[2],by=2)]
# copy the the original attributes
  exp1.tx <- attrcp(amt,exp1.tx)
  class(exp1.tx) <- class(mt)

  eval(parse(paste('X$',season,' <- txm.ds',sep='')))
  eval(parse(paste('Y$',season,' <- txs.ds',sep='')))
}

attributes(Y) <- NULL;
dim(Y) <- dim(X)
attr(X,'standard.deviation') <- Y
# add some new attributes describing the results:
attr(X,'description') <- 'cross-validation'
attr(X,'experiment') <- 'CORDEX ESD experiment 1'
attr(X,'method') <- 'esd'
attr(X,'url') <- 'https://github.com/metno/esd'
attr(X,'information') <- 'PCAs used to represent seasonal mean and standard deviation for all CLARIS stations, and DS applied to the PCs'
attr(X,'predictand_file') <- 'claris.Tx.rda'
attr(X,'predictor_file') <- c('ERAINT_olr_mon.nc','ERAINT_t2m_mon.nc','ERAINT_slp_mon.nc')
attr(X,'predictor_domain') <- 'lon=c(-90,-30),lat=c(-35,-15)'
attr(X,'history') <- history.stamp()
attr(X,'R-script') <- readLines('CORDEX.ESD.exp1.tx.R')

save(file='CORDEX.ESD.exp1.tx.esd.rda',X)
