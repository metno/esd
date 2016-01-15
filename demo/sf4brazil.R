## Rasmus E. Benestad. Oslo & Dhaka, 2015, 2016.
## This R-script reads and prepares gridded precipitation data as well as
## seasonal forecasts from the ECMWF (monthly fields). It then performs model
## output statistics (MOS) to downscale the seasonal forecasts. The strategy is
## to use each ensemble member in the regression analysis, which implies that
## the size of the calibration sample is ensemble size * number of years. The
## downscaling is applied to 3-month aggregated forecasts. The MOS makes use of PCA
## for representing the predictands (http://dx.doi.org/10.3402/tellusa.v67.28326).

## Clear the work space
rm(list=ls()); gc(reset=TRUE)

library(esd)
library(ncdf4)

## Copy the latest files from operational folder into storage
system(paste('cp sesong/globses/fcmean.mmsa.',format(Sys.time(),'%Y-%m'),
             '*.nc sesong/.',sep=''))

## Set the region:
test <- FALSE
lon <- c(-180,-30)
lat <- c(-50,30)

### This section contains functions --------------------------------------

upper <- function(x) mean(x,na.rm=TRUE) + 2*sd(x,na.rm=TRUE)
lower <- function(x) mean(x,na.rm=TRUE) - 2*sd(x,na.rm=TRUE)



## Reads the netCDF files with forecasts from ECMWF. Assumes the files contain
## drift-corrected anomalies.
readfcst <- function(fname,lon=NULL,lat=NULL,verbose=FALSE) {
  if (verbose) print(paste('readfcst:',fname))
  ncid <- nc_open(fname)
  y <- ncvar_get(ncid)
  yatt <- ncatt_get( ncid, names(ncid$var))
  for (i in 1:length(yatt))
    attr(y,names(yatt)[i]) <- yatt[[i]]
  #str(y)
  dims <- names(ncid$dim)
  for (i in 1:length(dims)) {
    x <- ncvar_get(ncid,dims[i])
    attr(y,dims[i]) <- x
  }
  nc_close(ncid)
  attr(y,'time') <- as.Date(attr(y,'time')/24 + julian(as.Date("1900-01-01")))
  #print(names(attributes(y)))
  
  ## Select subset:
  if (!is.null(lon)) 
    ix <- (attr(y,'longitude') >= min(lon)) &
          (attr(y,'longitude') <= max(lon)) else
    ix <- is.finite(attr(y,'longitude'))
  if (!is.null(lat)) 
    iy <- (attr(y,'latitude') >= min(lat)) &
          (attr(y,'latitude') <= max(lat)) else
    iy <- is.finite(attr(y,'latitude'))
  z <- y[ix,iy,,]
  z <- attrcp(y,z)
  attr(z,'time') <- as.Date(attr(y,'time'))
  attr(z,'longitude') <- attr(y,'longitude')[ix]
  attr(z,'latitude') <- attr(y,'latitude')[iy]
  attr(z,'dimensions') <- dim(z)
  attr(z,'file_name') <- fname
  attr(z,'source') <- 'ECMWF seasonal forecast'
  attr(z,'variable') <- 't2m'
  attr(z,'greenwich') <- FALSE
  attr(z,'unit') <- 'deg C'
  attr(z,'aspect') <- 'anomaly'
  class(z) <- 'forecast'
  if (verbose) print(names(attributes(z)))
  return(z)
}

## Combine the data from individual hind/forecast files into one hypermatrix
## Keep the same year+season, but divide the members by days.
## 1 season = 3 months = 90 different days.
combine.forecast <- function(x,y,verbose=FALSE) {
  if (verbose) print('combine.forecast:')
  dx <- attr(x,'dimensions')
  if (length(dx)==4) dx <- c(dx,1)
  dim(x) <- c(dx[1]*dx[2]*dx[3]*dx[4],dx[5])
  dy <- attr(y,'dimensions')
  if (verbose) {
    print(paste(dx,dy,collapse=' - '))
    print(attr(x,'time'))
    print(attr(y,'time'))
  }
  if (dx[3] != dy[3]) {
    print(paste('Change in number of members:',dx[3],'->',dy[3]))
    if (verbose) print('Take a subset of y:')
    if (verbose) print(dim(y))
    ty <- attr(y,'time')
    y <- y[,,1:dx[3],]; dy[3] <- dx[3]
    attr(y,'time') <- ty
    if (verbose) print(dim(y))
  }
  if (length(dy)==4) dy <- c(dy,1)
  dim(y) <- c(dy[1]*dy[2]*dy[3]*dy[4],dy[5])
  if (verbose) print(paste('cbind',dim(x),dim(y),collapse=' - '))
  z <- cbind(x,y)
  dim(z) <- c(dy[1],dy[2],dy[3],dy[4],dx[5]+dy[5])
  for (i in 1:length(attributes(x))) {
    attnm <- names(attributes(x))[i]
    if (substr(attnm,1,3) != 'dim') attr(z,attnm) <- attributes(x)[[i]]
  }
  attr(z,'time') <- as.Date(cbind(attr(x,'time'),attr(y,'time')))
  if (verbose) print(attr(z,'time'))
  attr(z,'dimensions') <-c(dy[1],dy[2],dy[3],dy[4],dx[5]+dy[5])
  attr(z,'file_name') <- c(attr(x,'file_name'),attr(y,'file_name'))
  class(z) <- class(x)
  return(z)
}

## Convert the forecast hypermatrix to the 'esd' field object by selecting one
## forecast memeber from the ensemble and one lead time (lag)
forecast2field <- function(x,member=1,lag=1,verbose=FALSE) {
  if (verbose) print('forecast2field:')
  d <- attr(x,'dimensions')[-c(3,4)]
  attnm <- names(attributes(x))
  if (verbose) print(attnm)
  z <- x[,,member,lag,]
  dim(z) <- c(d[1]*d[2],d[3])
  if (verbose) print(attr(x,'time'))
  y <- zoo(t(z),order.by=as.Date(attr(x,'time')[lag,]))
  for (i in 1:length(attnm)) {
    if (sum(is.element(substr(attnm[i],1,3),c('dim','time')))==0 )
      attr(y,attnm[i]) <- attributes(x)[[i]] else
      if (verbose) print(paste('skip',attnm[i]))
  } 

  attr(y,'time') <- NULL
  attr(y,'dimensions') <- d
  attr(y,'lag') <- lag
  attr(y,'member') <- member
  class(y) <- c('field','monthly','zoo')
  return(y)
}



threemonths <- function(Z1,Z2,Z3,m1,m2,m3) {
    ## Combines the forecasts for monthly series and makes sure they
    ## are assigned same year.
    index(Z1) <- year(Z1)
    ## If the months cover two years, make sure they all refer to the
    ## year of the first month:
    if (m2 > m1) index(Z2) <- year(Z2)  else
                 index(Z2) <- year(Z2) - 1
    if (m3 > m1) index(Z3) <- year(Z3)  else
                 index(Z3) <- year(Z3) - 1
    #print(index(Z1)); print(index(Z2)); print(index(Z2)) 
    Z <- 0.333*(zoo(Z1)+zoo(Z2)+zoo(Z3))
    class(Z) <- class(Z1)
    Z <- attrcp(Z1,Z)
    invisible(Z)
}

## Extract the same three months as from seasonal forecasts from
## observations.
obsthreemonths <- function(x,mons,FUN='mean') {
  y1m <- as.monthly(subset(x,it=mons[1]),FUN=FUN)
                                        # pick the three months following the
  y2m <- as.monthly(subset(x,it=mons[2]),FUN=FUN)
                                        # forecast. But the months may spill
  y3m <- as.monthly(subset(x,it=mons[3]),FUN=FUN)
                                        # over to the next year (eg Dec-Jan-Feb)
  print('Aggregate over next three forecast months')
  m1 <- (1:12)[is.element(month.abb,mons[1])]
  m2 <- (1:12)[is.element(month.abb,mons[2])]
  m3 <- (1:12)[is.element(month.abb,mons[3])]
  Ys <- threemonths(y1m,y2m,y3m,m1,m2,m3)
  invisible(Ys)
}

downscaleSFC <- function(sa24,FUN='mean',param='T2M',
                         plot=FALSE,verbose=FALSE,n=3,
                         path='sesong',pattern='fcmean.mmsa.') {

  ## Process the predictors: --------------------------------------
  ## Retrieve the predictors from sf system 4 (MOS): from /opdata
  print('Read the predictors from seasonal forecasts')
  pred.list <- list.files(path=path,pattern=pattern,
                          full.names=TRUE)
  if (length(pred.list)==0)
    stop('Cannot find the files with seasonal forecasts')
    
  ## Extract T2M (or SLP) fields:
  pred.list <- pred.list[grep(param,pred.list)]

  ## Extract the last month:
  lf <- pred.list[length(pred.list)]
  pred.list <- pred.list[grep(substr(lf,nchar(lf)-9,nchar(lf)),pred.list)]

  ## Read all the hind/forecast files and combine into one array
  for (i in 1:length(pred.list)) {
    print(pred.list[i])
    if (i==1) X <- readfcst(pred.list[i],lon=lon,lat=lat) else {
      x <- readfcst(pred.list[i],lon=lon,lat=lat)
      X <- combine.forecast(X,x,verbose=verbose)
    }
  }

  print(as.Date(attr(X,'time')))

  ## Loop through the different ensemble members and make predictor

  dfcst <- dim(X)
  Nens <- dim(X)[3]
  nlag <- dim(X)[4]

  ## Lag 1 month:
  il <- 1
  print('Convert to field')
  for (im in 1:41) {
    Z1 <- forecast2field(X,member=im,lag=1)
    Z2 <- forecast2field(X,member=im,lag=2)
    Z3 <- forecast2field(X,member=im,lag=3)

    ## Which months are foreccasted?
    m1 <- month(Z1)[1]; m2 <- month(Z2)[1]; m3 <- month(Z3)[1]
    mons <- month.abb[c(m1,m2,m3)]

    Z <- threemonths(Z1,Z2,Z3,m1,m2,m3)
    print(paste(paste(mons,collapse='-'),'forecasts'))
    ## Get the seasonal aggregate
    print(index(Z))
    ## Each member will have one data point per year, but use the day
    ## to represent the different model members. The same days need to
    ## be found in the observations - by repeating them for each
    ## ensemble member.
    index(Z) <- as.Date(paste(year(Z),'-01-01',sep='')) + im -1

    print(paste('Combine ensemble members: different days to member',im))
    
    if (im==1) sf4specs <- Z else
    sf4specs <- c(zoo(sf4specs),zoo(Z))
  }

  
  ## Fix the attributes and turn the zoo object into a field object
  print('Set the attributes')
  class(sf4specs) <- class(Z)
  sf4specs <- attrcp(Z,sf4specs)
  attr(sf4specs,'dimensions') <-
    c(length(lon(sf4specs)),length(lat(sf4specs)),length(index(sf4specs)))
  print('Estimate the SA24s')
  eof.fcst <- EOF(sf4specs)
  class(eof.fcst)[3] <- "day" # a fix to allow for MOS for ensemble forecasting

  ## Check the results for the seasonal forecasts:
  if (plot) plot(eof.fcst)

  ## Process the predictands:---------------------------------------------

  ## Aggregate 3-month means from the stations:
  print('Extract predictands for the three forecast months')
  y1m <- as.monthly(subset(sa24,it=mons[1]),FUN=FUN)
                                        # pick the three months following the
  y2m <- as.monthly(subset(sa24,it=mons[2]),FUN=FUN)
                                        # forecast. But the months may spill
  y3m <- as.monthly(subset(sa24,it=mons[3]),FUN=FUN)
                                        # over to the next year (eg Dec-Jan-Feb)
  print('Aggregate over next three forecast months')

  Ys <- threemonths(y1m,y2m,y3m,m1,m2,m3)
  index(Ys) <- as.Date(paste(year(Ys),'-01-01',sep=''))
  attr(Ys,'dimensions') <- attr(y1m,'dimensions')
  if (verbose) print(index(Ys))

  ## Add copies of the observations, but with different days to match
  ## the dates of the different ensemble members:
  for (i in 1:attr(sf4specs,'member')) {
    if (i==1) Y <- Ys else
    Y <- c(zoo(Y),zoo(Ys))
    index(Ys) <- index(Ys) + 1
  }
  class(Y) <- class(sa24)
  Y <- attrcp(sa24,Y)
  attr(Y,'dimensions') <- c(attr(y1m,'dimensions')[1:2],length(index(Y)))

  ## Synchronise predictands to the predictors:
  ## Use data synchronised with the seasonal forecasts
  Y <- matchdate(Y,sf4specs)
  
## Express as EOFs - these are used for the downscaling
  print('EOF of predictand')
  eof.mu <- EOF(Y)
  class(eof.mu)[3] <- "day" # a fix to allow for MOS for ensemble forecasting

  ## Downscale the data using PCA data:
  print("Downscale the station data using seasonal fcst as predictor")
  Z <- DS(eof.mu,eof.fcst,n=n)
  predictand <- attr(Z,'pattern')
  predictand <- attrcp(eof.mu,predictand)
  predictand -> attr(Z,'pattern')
  if (plot) {dev.new(); plot(Z)}
  attr(Z,'eof.fcst') <- eof.fcst
  attr(Z,'mons') <- mons
  attr(Z,'sf4specs') <- sf4specs
  attr(Z,'Y') <- Y
  print("Calibration of downscaling model finished")
  invisible(Z)
}

display.sf <- function(x,FUN) {
  y1 <- subset(x,it=length(index(x)))
  y3 <- subset(x,it=length(index(x)))
  mm1 <- (1:12)[is.element(month.abb,attr(x,'mons')[1])] 
  mm3 <- (1:12)[is.element(month.abb,attr(x,'mons')[3])] 
  index(y1) <- as.Date(paste(year(y1),mm1,'01',sep='-'))
  if (mm3 < mm1) yr3 <- year(y1) + 1 else  yr3 <- year(y1)
  index(y3) <- as.Date(paste(yr3,mm3,'01',sep='-'))
  y <- c(y1,y3)
  y <- attrcp(y1,y)
  if (FUN=='sum') y <- y*3
  if (FUN=='wetfreq') y <- y*100
  attr(y,'variable') <- switch(FUN,
                               'sum'=expression(P[tot]),
                               'wetmean'=expression(mu),
                               'wetfreq'=expression(f[w]))
  attr(y,'unit') <- switch(FUN,
                               'sum'='mm',
                               'wetmean'='mm/day',
                               'wetfreq'="'%'")
  attr(y,'dimensions') <- attr(y1,'dimensions')
  class(y) <- class(y1)
  map(y,colbar=list(pal='precip'))
  figlab(paste(attr(x,'mons'),collapse='-'),ypos=0.99)
}


sf <- function(eof,anomaly=NULL,param='T2M',save=TRUE,FUN='mean',
               path='sesong',outpath='specs.SF',pattern='fcmean.mmsa.',n=3,
               verbose=FALSE,plot=FALSE) {
  
  ## Default: if temperature, plot anomalies, but for precip show
  ## the full values (anomaly is applied to daily data and the
  ## seasonal forecasts for temperature are provided as anomalies)
  if (verbose) print('sf')
  
## Downscaled seasonal forecasts based on MOS and PCA representation
## of predictands:
  Z <- downscaleSFC(eof,param=param,FUN=FUN,path=path,pattern=pattern,
                    n=n,plot=plot,verbose=verbose)

  sfc.mos <- as.field(predict(Z,newdata=attr(Z,'eof.fcst')))
  attr(sfc.mos,'mons') <- attr(Z,'mons')
  print(paste('Seasonal forecast for',end(sfc.mos)))

## Save the results for further processing
        
  rname <- paste(outpath,'/specs.results.',varid(eof)[1],sep='')
  print(paste('Saving results in ',path,'.',varid(eof)[1],'.rda',sep=''))
  if (save) save(file=paste(rname,'.rda',sep=''),Z,sfc.mos)

## Use aggregate to extract statistics such as the mean, quantile, etc.
  fc.mean <- aggregate(sfc.mos,year,'mean')
#fc.upper <- aggregate(sfc.mos,year,'upper')
#fc.lower <- aggregate(sfc.mos,year,'lower')

## Plot the time series:
  dev.new()
  display.sf(sfc.mos,FUN)
  dev.copy2pdf(file=paste(rname,'.hcst.pdf',sep=''))
  
  return(sfc.mos)
}

## End of the section with functions -------------------------------

##-------------------------------------------------------------------
## Carry out the analysis
## Retrieve the predictands: station data -> PCA

## Retrieve the predictand data from PSD:
url <- 'ftp://ftp.cdc.noaa.gov/Datasets.other/south_america/sa24.daily.1.1940-2012.nc'

## Only download if the data is not stored locally
if (!file.exists('~/Downloads/sa24daily1940-2012.nc')) {
  print(paste('Download the predictand from',url))
  download.file(url,'~/Downloads/sa24daily1940-2012.nc')
}

## retrieve the predictand
print('Gridded precipitation')
sa24 <- retrieve('~/Downloads/sa24daily1940-2012.nc',param='precip',
                 lon=c(-50,-32),lat=c(-20,0))

print('Select season and time interval')
sa24 <- subset(sa24,it=c(1940,2011))

print('MOS for precipitation sum')
pt <- sf(sa24,FUN='sum',plot=TRUE)

print('MOS for wet-day mean precipitation')
mu <- sf(sa24,FUN='wetmean')

print('MOS for precipitation wet-day frequency')
fw <- sf(sa24,param='MSL',FUN='wetfreq')
