# Script for setting up and running CORDEX ESD experiment 1 for Tmax
# R.E. Benestad

# Number of PCAs used to represent the predictand data- determines the degree
# of detail but also the robustness of the results
#npca <- 7; lon <- c(-90,-30); lat <- c(-35,-15); eofs <- 1:10 # m:r=[];s:r=[]
#npca <- 3; lon <- c(-90,-30); lat <- c(-35,-15); eofs <- 1:5 # m:r=[];s:r=[]
#npca <- 15; lon <- c(-90,-30); lat <- c(-35,-15); eofs <- 1:10 # m:r=[];s:r=[]
#npca <- 7; lon <- c(-90,-30); lat <- c(-35,-15); eofs <- 1:10 # m:r=[];s:r=[]
#npca <- 3; lon <- c(-70,-15); lat <- c(-37,-17); eofs <- 1:5  # m: r=[]; fw: r= []
npca <- 17; lon <- c(-70,-15); lat <- c(-37,-17); eofs <- 1:5  # m: r=[0.653,0.714]; s: r= [0.496,0.774]

# load the predictands: CLARIS precip
print('Get the CLARIS data')
load('claris.Tx.rda')

# The location names of the last stations contain unsupported characters
attr(Tx,'location')[77:81] <- c("Aerodromo de Pedro Juan Caballero","Aerodromo de Concepcion",
                                "Villarrica del Espedritu Santo","Aerodromo de Pilar",
                                "Encarnacion")

# The experiment protocol: http://wcrp-cordex.ipsl.jussieu.fr/images/pdf/guidelines/CORDEX_ESD_Experiment1.pdf
# Tier 1: calibration: 1979-1995; validation: 1996-2006 
# Tier 2: 5-fold cross-validation: [1979,1983], [1984,1988],[1989,1993],[1994,1998],[1999,2003] 

# Limit to the prescribed interval
Tx <- subset(Tx,it=c(1979,2007))
#Tx0 <- Tx # Keep original copy

## Daily anomalies:
#print('Take the anomalies')
#Tx <- anomaly(Tx0)
#clim <- Tx0 - Tx # climatology

# Process the precipitation - predictand as annual mean and annual standard deviation:
# Perhaps change to use seasonal rather than annual?
print('Estimate seasonal statistics')
Mt4s <- as.4seasons(Tx,FUN='mean',nmin=30)
St4s <- as.4seasons(anomaly(Tx),FUN='sd',nmin=30)

#Tx0 <- Tx # Keep original copy
print('Take the anomalies')
Mt4s0 <- Mt4s
Mt4s <- anomaly(Mt4s0)
clim <- Mt4s0 - Mt4s  # climatology

# retrieve the predictors
print('Get the predictor data')
# Out-going long-wave radiation
olr <- as.4seasons(retrieve('data/ERAINT/ERAINT_olr_mon.nc',
                       lon=lon,lat=lat),FUN='mean')
attr(olr,'unit') <- "W * m**-2"

# Temperature
t2m <- as.4seasons(retrieve('data/ERAINT/ERAINT_t2m_mon.nc',
                            lon=lon,lat=lat))

# Mean sea level pressure
slp <- as.4seasons(retrieve('data/ERAINT/ERAINT_slp_mon.nc',
                       lon=lon,lat=lat),FUN='mean')

olr <- subset(olr,it=c(1979,2007))
slp <- subset(slp,it=c(1979,2007))
t2m <- subset(t2m,it=c(1979,2007))

X <- list(description='cordex-esd-exp1')

# Loop through the seasons:
for (season in c('djf','mam','jja','son')) {

  ## Prepare the predictors: use all available data
  eof.slp <- EOF(subset(slp,it=season))
  eof.olr <- EOF(subset(olr,it=season))
  eof.t2m <- EOF(subset(t2m,it=season))

  ## Extract the predictands for the specified season
  mt4s <- subset(Mt4s,it=season)
  st4s <- subset(St4s,it=season)

  ## Extract the predictands for the calibration period:
  mt4sc <- subset(mt4s,it=c(1979,1995))
  st4sc <- subset(st4s,it=c(1979,1995))

  ## Missing values are a problem for PCA - replace missing anomalies with zero
  miss.mc <- !is.finite(coredata(mt4sc))
  coredata(mt4sc)[miss.mc] <- 0
  miss.sc <- !is.finite(coredata(st4sc))
  coredata(st4sc)[miss.sc] <- 0

  ## Combine the local predictand information in the form of PCAs
  pca.mt <- PCA(mt4sc,neofs=npca)
  pca.st <- PCA(st4sc,neofs=npca)

  print(paste('Tear 1: Season:',season,' calibrate on the period:',
              start(pca.mt),'-',end(pca.mt)))
  z.mt1 <- DS(pca.mt,eof.t2m,rmtrend=FALSE,m=NULL,verbose=FALSE,eofs=eofs)
  z.st1 <- DS(pca.st,list(t2m=eof.t2m,slp=eof.slp,olr=eof.olr),
             rmtrend=FALSE,m=NULL,verbose=FALSE,eofs=eofs)

  ## Predict the results.
  mt.tier1 <-as.station(predict(z.mt1,newdata=eof.t2m,verbose=FALSE))
  st.tier1 <-as.station(predict(z.st1,newdata=list(t2m=eof.t2m,slp=eof.slp,olr=eof.olr)))
  
  print(paste('Tear 2: Season:',season))

  ## Missing values are a problem for PCA - replace missing anomalies with zero
  miss.m <- !is.finite(coredata(mt4s))
  coredata(mt4s)[miss.m] <- 0
  miss.s <- !is.finite(coredata(st4s))
  coredata(st4s)[miss.s] <- 0
  
  ## Combine the local predictand information in the form of PCAs
  pca.mt <- PCA(mt4s,neofs=npca)
  pca.st <- PCA(st4s,neofs=npca)
  
  ## The experiment results are provided by the cross-validation defined by
  ## parameter 'm' (used in crossval).
  ## Importan to set rmtrend=FALSE to retain original data in the cross-validation.
  ## First downscale the mean values:
  z.mt2 <- DS(pca.mt,eof.t2m,m='cordex-esd-exp1',rmtrend=FALSE,verbose=FALSE,eofs=eofs)

  ## The results are in the form of a PCA-object - convert back to a group of stations     
  txm.ds <- pca2station(z.mt2)
  exp1.txm <- pca2station(z.mt2,what='xval')

  ## Extract the predicted cross-validation results which follow the experiment:
  ## only grab the series of predicted values - not the original data used for calibration
  mt.tier2 <- pca2station(z.mt2,what='xval',verbose=TRUE)

  # Repeat the downscaling for the standard deviation: use a mix of predictors.
  z.st2 <- DS(pca.st,list(t2m=eof.t2m,slp=eof.slp,olr=eof.olr),
               m='cordex-esd-exp1',rmtrend=FALSE,verbose=FALSE,eofs=eofs)
  txs.ds <- pca2station(z.st2)
  exp1.txs <- pca2station(z.st2,what='xval')

  ## Reset missing values to missing values
  coredata(mt4sc)[miss.mc] <- NA
  coredata(st4sc)[miss.sc] <- NA
  coredata(mt4s)[miss.m] <- NA
  coredata(st4s)[miss.s] <- NA
  
  # Check: Figure: scatter plot
  
  x <- matchdate(subset(mt4s,it=season),txm.ds)
  dev.new(width=5,height=9)
  par(bty='n',las=1,oma=rep(0.25,4),mfcol=c(2,1),cex=0.5)
  plot(coredata(x),coredata(txm.ds),
       pch=19,col=rgb(1,0.2,0.2,0.25),
       xlab=expression(paste('Observed ',T[2*m],(degree*C))),
       ylab=expression(paste('Downscaled ',T[2*m],(degree*C))),
       main=paste(toupper(season),' mean temperature'),
       sub=paste('predictand: CLARIS; #PCA=',npca))
  grid()

  if (season(exp1.txm)[1]=='djf') exp1.txm <- subset(exp1.txm,it=c(1980,2004)) else
                                  exp1.txm <- subset(exp1.txm,it=c(1979,2003))
  x <- c(coredata(matchdate(x,exp1.txm)))
  y <- c(coredata(exp1.txm))
  points(x,y,col=rgb(1,0,0,0.25),lwd=2)
  abline(lm(y ~ x),col=rgb(0.8,0,0),lty=2)
  ok <-  is.finite(x) & is.finite(y)
  mtext(side=4,paste('r=',round(cor(x[ok],y[ok]),3)),las=3)

  
  # Check: Figure: scatter plot
  x <- matchdate(subset(st4s,it=season),txm.ds)
  plot(coredata(x),coredata(txs.ds),
       pch=19,col=rgb(1,0.2,0.2,0.25),
       xlab=expression(paste('Observed ',T[2*m],(degree*C))),
       ylab=expression(paste('Downscaled ',T[2*m],(degree*C))),
       main=paste(toupper(season),' standard deviation'),
       sub=paste('predictand: CLARIS; #PCA=',npca))
  grid()

  if (season(exp1.txm)[1]=='djf') exp1.txs <- subset(exp1.txs,it=c(1980,2004)) else
                                  exp1.txs <- subset(exp1.txs,it=c(1979,2003))
  x <- c(coredata(matchdate(x,exp1.txs)))
  y <- c(coredata(exp1.txs))
  points(x,y,col=rgb(1,0,0,0.25),lwd=2)
  abline(lm(y ~ x),col=rgb(0.8,0,0),lty=2)
  ok <-  is.finite(x) & is.finite(y)
  mtext(side=4,paste('r=',round(cor(x[ok],y[ok]),3)),las=3)

  # The independent validation is contained in exp1.txm (seasonal mean)
  # and exp1.txm (seasonal standard deviation)
  eval(parse(text=paste('X$mean.tier1.',season,' <- mt.tier1',sep='')))
  eval(parse(text=paste('X$sd.tier1.',season,' <- st.tier1',sep='')))
  eval(parse(text=paste('X$mean.tier2.',season,' <- exp1.txm',sep='')))
  eval(parse(text=paste('X$sd.tier2.',season,' <- exp1.txs',sep='')))
  eval(parse(text=paste('X$mean.obs.',season,' <- mt4s',sep='')))
  eval(parse(text=paste('X$sd.obs.',season,' <- st4s',sep='')))

}

print("Finished looping and downscaling; organise the data using rbind:")
## Tier 1
txm1.4s <- c(X$mean.tier1.djf,X$mean.tier1.mam,
             X$mean.tier1.jja,X$mean.tier1.son)

X$tier1 <- subset(txm.tier1,it=c(1996,2006))

## Tier 2
txm2.4s <- c(X$mean.tier2.djf,X$mean.tier2.mam,
             X$mean.tier2.jja,X$mean.tier2.son)
t2 <- index(txm2.4s)
txm.tier2 <- zoo(txm2.4s,order.by=t2) + matchdate(clim,it=t2)
txm.tier2 <- attrcp(X$mean.tier2.djf,txm.tier2)
class(txm.tier2) <- class(X$mean.tier2.djf)
X$tier2 <- subset(txm.tier2,it=c(1979,2003))

## Observed temperature
obs <- c(X$mean.obs.djf,X$mean.obs.mam,
         X$mean.obs.jja,X$mean.obs.son)
t0 <- index(obs)
txm.obs <- zoo(obs,order.by=t0) + matchdate(clim,it=t0)
txm.obs <- attrcp(X$mean.obs.djf,txm.obs)
class(txm.obs) <- class(X$mean.obs.djf)
x$obs <- txm.obs

# add some new attributes describing the results:
print('add new attributes')
attr(X,'description') <- 'cross-validation'
attr(X,'experiment') <- 'CORDEX ESD experiment 1'
attr(X,'method') <- 'esd'
attr(X,'url') <- 'https://github.com/metxo/esd'
attr(X,'information') <- 'PCAs used to represent seasonal mean and standard deviation for all CLARIS stations, and DS applied to the PCs'
attr(X,'predictand_file') <- 'claris.Tx.rda'
attr(X,'predictor_file') <- c('ERAINT_olr_mon.nc','ERAINT_t2m_mon.nc','ERAINT_slp_mon.nc')
attr(X,'predictor_domain') <- 'lon=c(-90,-30),lat=c(-35,-15)'
attr(X,'history') <- history.stamp()
attr(X,'R-script') <- readLines('CORDEX.ESD.exp.1.tx.R')

print('Save the results')
save(file='CORDEX.ESD.exp1.tx.esd.rda',X)

print('finished')
