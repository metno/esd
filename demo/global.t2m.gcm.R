## Re-compute the global mean
library(esd)

globalmean <- function(path='CMIP5.monthly/rcp45',
                     pattern='tas_',select=NULL,lon=NULL,lat=NULL) {

  fnames <- list.files(path=path,pattern=pattern,full.name=TRUE)
  fnms <- list.files(path=path,pattern=pattern)
  if (!is.null(select)) {fnames <- fnames[select]; fnms <- fnms[select]}
  n <- length(fnames)
  X <- matrix(rep(NA,n*240),n,240)
  yr <- 1861:2100
  meta <- list()
  for (i in 1:n) {
    gcm <- retrieve(fnames[i],lon=lon,lat=lat)
    gcmnm <- attr(gcm,'model_id')
    run <- attr(gcm,'realization')
    d <- attr(gcm,'dimensions')
    cal <- attr(gcm,'calendar')
    print(paste(i,n,gcmnm,run,paste(d,collapse='-')))
    y <- aggregate.area(annual(gcm),FUN='mean')
    i1 <- is.element(yr,year(y))
    i2 <- is.element(year(y),yr)
    ya <- anomaly(y,ref=1961:1990)
    if (i==1) plot(ya) else lines(ya)
    X[i,i1] <- coredata(ya)[i2]
    gcmnm <- gsub('-','.',gcmnm)
    cline <- paste('meta$',fnms[i],
             ' <- list(GCM=attr(gcm,"model_id"),run=run,d=d,calendar=cal,',
             'mean=attr(ya,"climatology"))',sep='')
    eval(parse(text=cline))
  }

  global.t2m.cmip5 <- zoo(t(X),order.by=yr)
  attr(global.t2m.cmip5,'metadata') <- meta
  attr(global.t2m.cmip5,'aspect') <- 'anomalies'
  attr(global.t2m.cmip5,'baseline') <- '1961-1990'
  attr(global.t2m.cmip5,'experiment_id') <-  attr(gcm,'experiment_id')
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
    cal <- attr(gcm,'calendar')
    print(paste(i,n,gcmnm,run,paste(d,collapse='-')))
    y <- aggregate(gcm,month,FUN='mean')
    Y <- combine(Y,y)
    gcmnm <- gsub('-','.',gcmnm)
    cline <- paste('meta$',gcmnm,'.',run,
             ' <- list(GCM=attr(gcm,"model_id"),run=run,d=d,calendar=cal,dT=dT)',sep='')
    eval(parse(text=cline))
  }

  attr(Y,'metadata') <- meta
  attr(Y,'aspect') <- 'mean-seasonal-cycle'
  attr(Y,'baseline') <- range(year(rea))
  attr(Y,'experiment_id') <-  attr(gcm,'experiment_id')
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

if (TRUE) {
  global.t2m.cmip5.rcp45 <- globalmean(path='CMIP5.monthly/rcp45')
  global.t2m.cmip5.rcp85 <- globalmean(path='CMIP5.monthly/rcp85')
  global.t2m.cmip5.rcp26 <- globalmean(path='CMIP5.monthly/rcp26')
  ##global.t2m.cmip3.sresa1b <- globmean(path='CMIP3.monthly/SRESA1b')

  reanalysis <- aggregate.area(annual(retrieve('air.mon.mean.nc')),FUN='mean')
  obs <- anomaly(reanalysis,ref=1961:1990)
  index(obs) <- year(obs)
  
  data(global.t2m.cmip3)
  global.t2m.cmip3 <- global.t2m.cmip3 - mean(window(global.t2m.cmip3,start=1961,end=1990))

  global.t2m.gcm <- list(global.t2m.cmip5.rcp45=global.t2m.cmip5.rcp45,
                         global.t2m.cmip5.rcp85=global.t2m.cmip5.rcp85,
                         global.t2m.cmip5.rcp26=global.t2m.cmip5.rcp26,
                         global.t2m.cmip3.sresa1b=global.t2m.cmip3)
  attr(global.t2m.gcm,'obs') <- obs

  save(file='global.t2m.gcm.rda',global.t2m.gcm)
}
