## Re-compute the global mean
library(esd)

globalmean <- function(path='CMIP5.monthly/rcp45',ref=1961:1990,usefnames=TRUE,
                       annual=TRUE,pattern='tas_',param='tas',relative=FALSE,
                       select=NULL,lon=NULL,lat=NULL,FUN='mean',anomaly=TRUE) {

  fnames <- list.files(path=path,pattern=pattern,full.name=TRUE)
  fnms <- list.files(path=path,pattern=pattern)
  if (!is.null(select)) {fnames <- fnames[select];fnms <- fnms[select]}
  if (!usefnames) fnms <- paste('gcm',1:length(fnames),sep='.')
  n <- length(fnames)
  if (annual) nt <- 240 else nt <- 240*12
  X <- matrix(rep(NA,n*nt),n,nt)
  if (annual) yr <- 1861:2100 else
              yr <- sort(rep(1861:2100,12)) + round((rep(1:12,240)-0.5)/12,2)
  meta <- list()
  for (i in 1:n) {
    gcm <- retrieve(fnames[i],param=param,lon=lon,lat=lat)
    lon.rng <- range(lon)
    lat.rng <- range(lat)
    gcm <- regrid(gcm,is=list(lon=seq(lon.rng[1],lon.rng[2],by=1),lat=seq(lat.rng[1],lat.rng[2],by=1))) # Added AM 21.02.2017 - All gcms must be at the same grid
    gcmnm <- attr(gcm,'model_id')
    run <- attr(gcm,'realization')
    rip <- attr(gcm,'parent_experiment_rip')
    d <- attr(gcm,'dimensions')
    cal <- attr(gcm,'calendar')
    print(paste(i,n,gcmnm,run,paste(d,collapse='-'),
                min(year(gcm)),max(year(gcm)),fnames[i]))
    if (annual) gcm <- annual(gcm)
    y <- aggregate.area(gcm,FUN=FUN)
    if (annual) {
      i1 <- is.element(yr,year(y))
      i2 <- is.element(year(y),yr)
    } else {
      i1 <- is.element(yr,year(y)+round((month(y)-0.5)/12,2))
      i2 <- is.element(year(y) + round((month(y)-0.5)/12,2),yr)
    }
    if (anomaly) ya <- anomaly(y,ref=ref) else ya <- y
    if(relative) {
      ya <- ya/attr(ya,"climatology")*100
      attr(ya,"unit") <- "%"
    }
    if (i==1) plot(ya) else lines(ya)
    X[i,i1] <- coredata(ya)[i2]
    gcmnm <- gsub('-','.',gcmnm)
    cline <- paste('meta$',fnms[i],
             ' <- list(GCM=gcmnm,run=run,rip=rip,d=d,calendar=cal,',
             'mean=attr(ya,"climatology"))',sep='')
    eval(parse(text=cline))
  }
  if (!annual) yr <- as.Date(paste(trunc(yr),round(12*(yr-trunc(yr))+0.5),'01',sep='-'))
  global.t2m.cmip5 <- zoo(t(X),order.by=yr)
  attr(global.t2m.cmip5,'metadata') <- meta

  attr(global.t2m.cmip5,'aspect') <- 'anomalies'
  attr(global.t2m.cmip5,'baseline') <- ref

  if(relative) {
    attr(global.t2m.cmip5,"unit") <- "%"
    attr(global.t2m.cmip5,"aspect") <- "relative anomaly"
  } else {
    attr(global.t2m.cmip5,'unit') <- attr(gcm,'unit')
    attr(global.t2m.cmip5,'aspect') <- 'anomalies'
  }
  attr(global.t2m.cmip5,'experiment_id') <-  attr(gcm,'experiment_id')
  attr(global.t2m.cmip5,'variable') <- attr(gcm,'variable')
  attr(global.t2m.cmip5,'history') <- match.call()
  invisible(global.t2m.cmip5)
}


AC <- function(path='CMIP5.monthly/rcp45',
               pattern='tas_',lon=c(-50,30),lat=c(40,70),
               reanalysis='air.mon.mean.nc',select=NULL) {

  fnames <- list.files(path=path,pattern=pattern,full.name=TRUE)
  if (!is.null(select)) fnames <- fnames[select]
  rea <- retrieve(reanalysis,lon=lon,lat=lat)
  Y <- aggregate(rea,month,FUN='mean')
  n <- length(fnames)
  meta <- list()
  for (i in 1:n) {
    gcm <- retrieve(fnames[i],lon=lon,lat=lat)
    d <- attr(gcm,'dimensions')
    y <- annual(aggregate.area(gcm,FUN='mean'))
    dT <- mean(window(zoo(y),start=2070,end=2099)) -
          mean(window(zoo(y),start=1961,end=1990))    
    gcm <- subset(gcm,it=range(year(rea)))
    gcmnm <- attr(gcm,'model_id')
    run <- attr(gcm,'realization')
    rip <- attr(gcm,'parent_experiment_rip')
    cal <- attr(gcm,'calendar')
    print(paste(i,n,gcmnm,run,paste(d,collapse='-')))
    y <- aggregate(gcm,month,FUN='mean')
    Y <- combine(Y,y)
    gcmnm <- gsub('-','.',gcmnm)
    cline <- paste('meta$',gcmnm,'.',run,
             ' <- list(GCM=attr(gcm,"model_id"),run=run,rip=rip,d=d,calendar=cal,dT=dT)',sep='')
    eval(parse(text=cline))
  }

  attr(Y,'metadata') <- meta
  attr(Y,'aspect') <- 'mean-seasonal-cycle'
  attr(Y,'baseline') <- range(year(rea))
  attr(Y,'experiment_id') <-  attr(gcm,'experiment_id')
  attr(Y,'variable') <- attr(gcm,'variable')
  attr(Y,'unit') <- attr(gcm,'unit')
  attr(Y,'history') <- match.call()
  invisible(Y)
}

## Examples of how to call these functions:

if (FALSE) {
  annual.cycle.cmip5 <- AC()
  ceof <- EOF(annual.cycle.cmip5)
  plot(ceof)
  save(file='annual.cycle.cmip5.rda',annual.cycle.cmip5)
}

if (FALSE) {
  global.t2m.cmip5.rcp45 <- globalmean(path='CMIP5.monthly/rcp45')
  global.t2m.cmip5.rcp85 <- globalmean(path='CMIP5.monthly/rcp85')
  global.t2m.cmip5.rcp26 <- globalmean(path='CMIP5.monthly/rcp26')
  ##global.t2m.cmip3.sresa1b <- globalmean(path='CMIP3.monthly/SRESA1b')

  reanalysis <- aggregate.area(annual(retrieve('air.mon.mean.nc')),FUN='mean')
  obs <- anomaly(reanalysis,ref=1961:1990)
  index(obs) <- year(obs)
  
#  data(global.t2m.cmip3)
#  global.t2m.cmip3 <- global.t2m.cmip3 - mean(window(global.t2m.cmip3,start=1961,end=1990))

  global.t2m.gcm <- list(global.t2m.cmip5.rcp45=global.t2m.cmip5.rcp45,
                         global.t2m.cmip5.rcp85=global.t2m.cmip5.rcp85,
                         global.t2m.cmip5.rcp26=global.t2m.cmip5.rcp26,
                         global.t2m.cmip3.sresa1b=global.t2m.cmip3)
  attr(global.t2m.gcm,'obs') <- obs

  save(file='global.t2m.gcm.rda',global.t2m.gcm)
}

if (FALSE) {
  ## Regional area mean temperature
  t2m.cmip3.sresa1b <- globalmean(path='CMIP3.monthly/SRESA1b',lon=c(-30,30),lat=c(50,70),ref=c(2000,2015),
                                  select=-c(7:10,42:50),usefnames=FALSE)
  ## Some of the CMIP3 runs do not follow standard format and make the code crash - exclude those.
  t2m.cmip5.rcp45 <- globalmean(path='CMIP5.monthly/rcp45',lon=c(-30,30),lat=c(50,70),ref=c(2000,2015))
  t2m.cmip5.rcp85 <- globalmean(path='CMIP5.monthly/rcp85',lon=c(-30,30),lat=c(50,70),ref=c(2000,2015))
  t2m.cmip5.rcp26 <- globalmean(path='CMIP5.monthly/rcp26',lon=c(-30,30),lat=c(50,70),ref=c(2000,2015))
  reanalysis <- aggregate.area(annual(retrieve('air.mon.mean.nc',lon=c(-30,30),lat=c(50,70))),FUN='mean')
  obs <- anomaly(reanalysis,ref=2000:2015)
  index(obs) <- year(obs)
    t2m.gcm <- list(t2m.cmip5.rcp45=t2m.cmip5.rcp45,
                  t2m.cmip5.rcp85=t2m.cmip5.rcp85,
                  t2m.cmip5.rcp26=t2m.cmip5.rcp26,
                  t2m.cmip3.sresa1b=t2m.cmip3.sresa1b)
  attr(t2m.gcm,'obs') <- obs
                  
  save(file='t2m.gcm.rda',t2m.gcm)
  graph(t2m.gcm)
  summary(t2m.gcm)
  lines(index(t2m.gcm[[1]]),apply(coredata(t2m.gcm[[1]]),1,'mean'),lwd=3,col='wheat')
  lines(index(t2m.gcm[[2]]),apply(coredata(t2m.gcm[[2]]),1,'mean'),lwd=3,col='red')
  lines(index(t2m.gcm[[3]]),apply(coredata(t2m.gcm[[3]]),1,'mean'),lwd=3,col='green')
  lines(index(t2m.gcm[[4]]),apply(coredata(t2m.gcm[[4]]),1,'mean',na.rm=TRUE),lwd=3,col='grey')
  legend(1850,4.5,c('RCP4.5','RCP8.5','RCP2.6','SRESA1b'),col=c('wheat','red','green','grey'),lty=1,lwd=3,bty='n')
}

arcticwarming <- function(presaved=TRUE,arcticrcp85data='t2m.70to90N.rcp85.rda') {
  if (presaved) {
    if (!file.exists(arcticrcp85data))
       download.file('https://ndownloader.figshare.com/files/5431400',arcticrcp85data)
    load(arcticrcp85data)
  } else t2m.70to90N.rcp85 <- globalmean(path='CMIP5.monthly/rcp85',lat=c(70,90),annual=FALSE)
  djf <- aggregate(subset(t2m.70to90N.rcp85,it='djf'),year,FUN='mean')
  plot(djf,plot.type='single',main='RCP8.5 mean winter temperature 70N-90N',
       xlab='Year',ylab=expression(degree*C))
  grid()
}
