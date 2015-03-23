# Script for setting up and running CORDEX ESD experiment 1 for Tmax
# R.E. Benestad

# Number of PCAs used to represent the predictand data- determines the degree
# of detail but also the robustness of the results
npca <- 4

# load the predictands: CLARIS precip
print('Get the CLARIS data')
load('claris.Tx.rda')

Tx0 <- Tx # Keep original copy

print('Take the anomalies')
Tx <- anomaly(Tx0)
clim <- Tx0 - Tx # climatology

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

X <- list(description='cordex-esd-exp1')

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
             m='cordex-esd-exp1',detrend=FALSE,verbose=FALSE)
# The results are in the form of a PCA-object - convert back to a group of stations            
  txm.ds <- pca2station(z.mt)
  
  # Check: Figure: scatter plot
  x <- matchdate(subset(mt4s,it=season),txm.ds)
  dev.new(width=5,height=9)
  par(bty='n',las=1,oma=rep(0.25,4),mfcol=c(2,1),cex=0.5)
  plot(coredata(anomaly(x)),coredata(anomaly(txm.ds)),
       pch=19,col=rgb(1,0,0,0.5),
       xlab=expression(paste('Observed ',T[2*m],(degree*C))),
       ylab=expression(paste('Downscaled ',T[2*m],(degree*C))),
       main=paste(toupper(season),' mean temperature'),
       sub=paste('predictand: CLARIS; #PCA=',npca))
  grid()

# Repeat the downscaling for the standard deviation: use a mix of predictors.
  z.st <- DS(pca.st,list(t2m=eof.t2m,slp=eof.slp,olr=eof.olr),
               m='cordex-esd-exp1',detrend=FALSE,verbose=TRUE)

  txs.ds <- pca2station(z.st)

  # Check: Figure: scatter plot
  x <- matchdate(subset(st4s,it=season),txm.ds)
  plot(coredata(anomaly(x)),coredata(anomaly(txs.ds)),
       pch=19,col=rgb(1,0,0,0.5),
       xlab=expression(paste('Observed ',T[2*m],(degree*C))),
       ylab=expression(paste('Downscaled ',T[2*m],(degree*C))),
       main=paste(toupper(season),' standard deviation'),
       sub=paste('predictand: CLARIS; #PCA=',npca))
  grid()
# Extract the predicted cross-validation results which follow the experiment:
  d <- dim(attr(txs.ds,'evaluation'))

# only grab the series of predicted values - not the original data used for calibration
  exp1.txm <- attr(txm.ds,'evaluation')[,seq(2,d[2],by=2)]
  print('add climatology to txm')
  exp1.txm <- exp1.txm + matchdate(clim,it=exp1.txm)
  exp1.txs <- attr(txs.ds,'evaluation')[,seq(2,d[2],by=2)]
# copy the the original attributes
  exp1.txm <- attrcp(mt4s,exp1.txm)
  exp1.txs <- attrcp(mt4s,exp1.txs)
  class(exp1.txm) <- class(mt4s)
  class(exp1.txs) <- class(mt4s)

  eval(parse(text=paste('X$txm.',season,' <- txm.ds',sep='')))
  eval(parse(text=paste('X$txs.',season,' <- txs.ds',sep='')))
  eval(parse(text=paste('X$txm0.',season,' <- mt4s',sep='')))
  eval(parse(text=paste('X$txs0.',season,' <- st4s',sep='')))
  
  # The independent validation is contained in exp1.txm (seasonal mean)
  # and exp1.txm (seasonal standard deviation)
  eval(parse(text=paste('X$exp1.txm.',season,' <- exp1.txm',sep='')))
  eval(parse(text=paste('X$exp1.txs.',season,' <- exp1.txs',sep='')))
}

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
attr(X,'R-script') <- readLines('CORDEX.ESD.exp.1.tx.R')

save(file='CORDEX.ESD.exp1.tx.esd.rda',X)
