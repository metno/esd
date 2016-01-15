## Rasmus E. Benestad. Dhaka, November 16-19, 2015.
## This R-script reads and prepares station data from the Bangladesh
## Meteorological Department (BMD) as well as seasonal forecasts from the
## EXMWF (monthly fields). It then performs model output statistics (MOS)
## to downscale the seasonal forecasts. The strategy is to use each ensemble
## member in the regression analysis, which implies that the size of the
## calibration sample is ensemble size * number of years. The downscaling
## is applied to 3-month aggregated forecasts. The MOS makes use of PCA for
## representing the predictands (http://dx.doi.org/10.3402/tellusa.v67.28326).
## (Based on sf4specs.R).

library(esd)
library(ncdf4)

## Set the region:

lat <- c(-30,30)
lon <- c(50,180)

### This section contains functions --------------------------------------

upper <- function(x) mean(x,na.rm=TRUE) + 2*sd(x,na.rm=TRUE)
lower <- function(x) mean(x,na.rm=TRUE) - 2*sd(x,na.rm=TRUE)

readBMD <- function(fname='~/Dropbox/BMD/rain',verbose=TRUE,exclude=NULL) {
  ## Read the data from the BMD and prepare as esd station object.
  param <- switch(substr(fname,nchar(fname)-3,nchar(fname)),
                  'rain'='precip','tmax'='tmax','tmin'='tmin')
  unit <- switch(param,
                  'precip'='mm/day','tmax'='degC','tmin'='degC')
  lname <- switch(param,
                  'precip'='precipitation',
                  'tmax'='maximum temperature',
                  'tmin'='minimum temperature')
  if (verbose) print(param)
  Y <- read.table(fname,na.strings="***",
                  col.names=c('stid','year','month','day',param))
  meta <- read.table('~/Dropbox/BMD/StList',header=TRUE,sep=',',as.is=TRUE)
  stations <- as.character(rownames(table(Y$stid)))
  il <- 0
  if (!is.null(exclude)) stations <- stations[!is.element(stations,exclude)]
  if (verbose) print(stations)
  for (is in stations) {
    it <- is.element(Y$stid,is)
    im <- is.element(as.character(meta$ID),is)
    y <- Y[it,]
    if (verbose) print(paste(is,meta$Name[im]))
    z <- zoo(y[[5]],
             order.by=as.Date(paste(y[[2]],y[[3]],y[[4]],sep='-')))
    cz <- coredata(z); cz[z < 0] <- NA; cz -> coredata(z)
    z <- as.station(z,param=param,unit=unit,
                    loc=meta$Name[im],stid=as.numeric(is),
                    lon=meta$LonDeg[im] + meta$LonMin[im]/60,
                    lat=meta$LatDeg[im] + meta$LatMin[im]/60,
                    alt=meta$Alt[im],longname=lname,cntr='Bangladesh')
    if (il==0) pr.bmd <- z else
               pr.bmd <- combine(pr.bmd,z)
    il <- il+1
  }
  return(pr.bmd)
}


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

gridmap <- function(Y,breaks=NULL,verbose=FALSE,colbar=NULL) {
  require(LatticeKrig)

  y <- apply(annual(Y,FUN='sum'),2,'mean',na.rm=TRUE)

  ## Get data on the topography on the 5-minute resolution
  data(etopo5)
  etopo5 <- subset(etopo5,
                   is=list(lon=range(lon(Y))+c(-1,1),
                           lat=range(lat(Y))+c(-1,1)))
  ## Mask the sea: elevations below 1m below sea level is masked.
  etopo5[etopo5<=-1] <- NA

  ## Set the grid to be the same as that of etopo5:
  grid <- structure(list(x=lon(etopo5),y=lat(etopo5)),class='gridList')

  ## Flag dubplicated stations:
  ok <- !(duplicated(lon(Y)) & duplicated(lat(Y)))

  ## Spread in the  90-percente interval changing
  obj <- LatticeKrig( x=cbind(lon(Y)[ok],lat(Y)[ok]),
                      y=y[ok],Z=alt(Y)[ok])

  ##  obj <- LatticeKrig( x=cbind(lon[ok],lat[ok]), y=z[2,ok],Z=alt[ok])
  w <- predictSurface(obj, grid.list = grid,Z=etopo5)
  w$z[is.na(etopo5)] <- NA

## Get rid of packages that have functions of same name:
  detach("package:LatticeKrig")
  detach("package:fields")
  detach("package:spam")
  detach("package:grid")
  detach("package:maps")

  ## Convert the results from LatticeKrig to esd:
  W <- w$z
  attr(W,'variable') <- varid(Y)[1]
  attr(W,'unit') <- unit(Y)[1]
  attr(W,'longitude') <- w$x
  attr(W,'latitude') <- w$y
  class(W) <- class(etopo5)

  ## Make a projection that zooms in on the Barents region

  rev <- switch(varid(Y)[1],'t2m'=FALSE,'precip'=TRUE)
  Wx <- max(abs(W),na.rm=TRUE)
  if (is.null(breaks)) breaks <- round(seq(-Wx,Wx,length=31),2) 
  map(W,xlim=range(lon(W)),ylim=range(lat(W)),
      colbar=colbar,verbose=verbose)
  invisible(W)
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

downscaleSFC <- function(tx.bmd,FUN='mean',param='T2M',
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
  if (length(pred.list)==0) {print('No matching seasonal forecasts')  ; browser()}

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
    
    if (im==1) sf4bmd <- Z else
    sf4bmd <- c(zoo(sf4bmd),zoo(Z))
  }

  
  ## Fix the attributes and turn the zoo object into a field object
  print('Set the attributes')
  class(sf4bmd) <- class(Z)
  sf4bmd <- attrcp(Z,sf4bmd)
  attr(sf4bmd,'dimensions') <-
    c(length(lon(sf4bmd)),length(lat(sf4bmd)),length(index(sf4bmd)))
  print('Estimate the EOFs')
  eof.fcst <- EOF(sf4bmd)
  class(eof.fcst)[3] <- "day" # a fix to allow for MOS for ensemble forecasting

  ## Check the results for the seasonal forecasts:
  if (plot) plot(eof.fcst)


  ## Process the predictands:---------------------------------------------

  ## Aggregate 3-month means from the stations:
  print('Extract predictands for the three forecast months')
  y1m <- as.monthly(subset(tx.bmd,it=mons[1]),FUN=FUN)
                                        # pick the three months following the
  y2m <- as.monthly(subset(tx.bmd,it=mons[2]),FUN=FUN)
                                        # forecast. But the months may spill
  y3m <- as.monthly(subset(tx.bmd,it=mons[3]),FUN=FUN)
                                        # over to the next year (eg Dec-Jan-Feb)
  print('Aggregate over next three forecast months')

  Ys <- threemonths(y1m,y2m,y3m,m1,m2,m3)
  index(Ys) <- as.Date(paste(year(Ys),'-01-01',sep=''))
  if (verbose) print(index(Ys))

  ## Add copies of the observations, but with different days to match
  ## the dates of the different ensemble members:
  for (i in 1:attr(sf4bmd,'member')) {
    if (i==1) Y <- Ys else
    Y <- c(zoo(Y),zoo(Ys))
    index(Ys) <- index(Ys) + 1
  }
  class(Y) <- class(tx.bmd)
  Y <- attrcp(tx.bmd,Y)

  ## Synchronise predictands to the predictors:
  ## Use data synchronised with the seasonal forecasts
  Y <- matchdate(Y,sf4bmd)

  ## Fill in missing data and estimate PCA of the station data:
  Y <- pcafill(Y)
  pca <- PCA(Y)
  class(pca)[3] <- "day" # a fix to allow for MOS for ensemble forecasting

  ## Downscale the data using PCA data:
  print("Downscale the station data using seasonal fcst as predictor")
  Z <- DS(pca,eof.fcst,n=n)
  if (plot) {dev.new(); plot(Z)}
  attr(Z,'eof.fcst') <- eof.fcst
  attr(Z,'mons') <- mons
  attr(Z,'sf4bmd') <- sf4bmd
  attr(Z,'Y') <- Y
  invisible(Z)
}

display.sf <- function(x,x0=NULL,is=1) {
  y <- subset(x,is=is)
  z <- subset(attr(x,'Y'),is=is)
  Y <- coredata(y)
  t <- year(y)
  par(bty='n')
  varid <- switch(varid(x)[1],'tmax'='T[max]','tmin'='T[min]',
                  't2m'='T[2*m]','pr'='precip','precip'='precip')
  if (unit(x)[1]=='degC') unit <- 'degree*C' else unit <- unit(x)[1]
  ylab <- eval(parse(text=paste('expression(',varid,
                       ' * ~(',unit[1],'))')))
  plot(t,y,cex=3,pch=19,col=rgb(0.5,0.5,0.5,0.15),
       main=loc(y),ylab=ylab,xlab='year')
  grid()
  if (!is.null(x0)) {
    y0 <- subset(x0,is=is)
    points(year(y0)+0.33,coredata(y0),cex=1.5,lwd=2,
           pch=4,col=rgb(0.25,0.25,0.75,0.15))
  }
  points(year(z),coredata(z),cex=3,pch=19)
}

## Estimate the daily mean temperature from the daily min and max
dailymeanT <- function(tx.bmd,tn.bmd,verbose=FALSE) {
  ## make sure that the mean is estimated using the same stations
  if (verbose) print('dailymeanT')
  isn <- is.element(stid(tn.bmd),stid(tx.bmd))
  isx <- is.element(stid(tx.bmd),stid(tn.bmd))
  tn.bmd <- subset(tn.bmd,is=isn)
  tx.bmd <- subset(tx.bmd,is=isx)
  if (verbose) {print(stid(tn.bmd)); print(stid(tx.bmd))}
  if (sum(stid(tn.bmd) - stid(tx.bmd)) != 0)
    stop('dailymeanT: stations are not corresponding')
  t2m <- 0.5*(zoo(tn.bmd)+zoo(tx.bmd))
  class(t2m) <- class(tn.bmd)
  t2m <- attrcp(tn.bmd,t2m)
  attr(t2m,'variable') <- 't2m'
  attr(t2m,'longname') <- 'daily mean temperature'
  attr(t2m,'info') <- '0.5*(Tn+Tx)'
  attr(t2m,'history') <- history.stamp(tn.bmd)
  invisible(t2m)
}
  
## Distribution plots histograms for terciles by default. 
distribution.sf <- function(x,x0=NULL,Y3m=NULL,col=c('blue','green','red'),
                            breaks=NULL,is=1,it=NULL) {
  y <- subset(x,is=is)
  if (is.null(it)) it <- year(y)[length(index(y))]
  y <- subset(y,it=it)
  par(bty='n')
    varid <- switch(varid(x)[1],'tmax'='T[max]','tmin'='T[min]',
                  't2m'='T[2*m]','pr'='precip','precip'='precip')
  if (unit(x)[1]=='degC') unit <- 'degree*C' else unit <- unit(x)[1]
  xlab <- eval(parse(text=paste('expression(',varid,
                       ' * ~(',unit[1],'))')))
  if (is.null(breaks)) {
    yx <- max(abs(coredata(y)),na.rm=TRUE)
    if (is.null(Y3m)) breaks <- round(seq(-yx,yx,length=7)) else {
      ## Use the observations to define the categories low, neutral, high
      obs <- subset(Y3m,is=is)
      breaks <- qnorm(p=c(0.1,0.333,0.666,0.9),
                      mean=mean(coredata(obs),na.rm=TRUE))
    }
  }
  y <- coredata(y)
  y[y < min(breaks)] <- min(breaks)
  y[y > max(breaks)] <- max(breaks)
  print(breaks); print(summary(y))
  hist(y,col=col,xlab=xlab,freq=FALSE,breaks=breaks,
       xlim=range(breaks),
       main=paste(paste(attr(x,'mons'),collapse='-'),it))
  if (!is.null(x0)) {
    y0 <- subset(x0,is=is)
    y0 <- coredata(subset(y0,it=it))
    y0[y0 < min(breaks)] <- min(breaks)
    y0[y0 > max(breaks)] <- max(breaks)
    h0 <- hist(y0,breaks=breaks,plot=FALSE)
    lines(h0$mids,h0$density,lwd=5,
          col=rgb(0.25,0.25,0.75,0.25))
  }
}

sf <- function(tx.bmd,is=1,anomaly=NULL,param='T2M',save=TRUE,FUN='mean',
               path='sesong',pattern='fcmean.mmsa.',n=3,
               verbose=FALSE,plot=FALSE) {
  
  ## Default: if temperature, plot anomalies, but for precip show
  ## the full values (anomaly is applied to daily data and the
  ## seasonal forecasts for temperature are provided as anomalies)
  if (verbose) print('sf')
  if (is.null(anomaly)) anomaly <- is.T(tx.bmd)
  if (anomaly) tx.bmd <- anomaly(tx.bmd)
  
## Downscaled seasonal forecasts based on MOS and PCA representation
## of predictands:
Z <- downscaleSFC(tx.bmd,param=param,FUN=FUN,path=path,pattern=pattern,
                  n=n,plot=plot,verbose=verbose)

sfc.mos <- as.station(predict(Z,newdata=attr(Z,'eof.fcst')))
attr(sfc.mos,'variable') <- varid(tx.bmd)
  
## Get the seasonal forecasts interpolated directly from the models if the
  ## predictand is the daily mean temperature:
if (varid(tx.bmd)[1]=='t2m')
  sfc.int <- regrid(attr(Z,'sf4bmd'),is=attr(Z,'Y')) else
  sfc.int <- NULL

## Save the results for furthr processing
  rname <- paste('sesong4bmd/bmd.results.',varid(tx.bmd)[1],sep='')
print(paste('Saving results in bmd.results.',varid(tx.bmd)[1],'.rda',sep=''))
if (save) save(file=paste(rname,'.rda',sep=''),Z,sfc.mos,sfc.int)

## Use aggregate to extract statistics such as the mean, quantile, etc.
fc.mean <- aggregate(sfc.mos,year,'mean')
#fc.upper <- aggregate(sfc.mos,year,'upper')
#fc.lower <- aggregate(sfc.mos,year,'lower')
  
colbar <- list(breaks=pretty(coredata(fc.mean),n=11),pal=varid(tx.bmd)[1])
  
dev.new()
gridmap(subset(fc.mean,it=length(index(fc.mean))),colbar=colbar)
points(lon(sfc.mos),lat(sfc.mos))
text(lon(sfc.mos),lat(sfc.mos),loc(sfc.mos),pos=1,cex=0.7)
lab <- paste('ECMWF MOS',paste(attr(Z,'mons'),collapse='-'),'starting in',
             year(fc.mean)[length(index(fc.mean))])
figlab(lab,xpos=0.35,ypos=0.99)
dev.copy2pdf(file=paste(rname,'.map.pdf',sep=''))

if (!is.null(sfc.int)) {
  print('Maps for the interpolated seasonal forecasts')
  fc0.mean <- aggregate(sfc.int,year,'mean')
  y0 <- subset(fc0.mean,it=length(index(fc.mean)))
  attr(y0,'longitude') <- lon(fc.mean)
  attr(y0,'latitude') <- lat(fc.mean)
  gridmap(y0,colbar=colbar)
  lab0 <- paste('ECMWF interpolated',paste(paste(attr(Z,'mons'),collapse='-')),
                'starting in',year(fc0.mean)[length(index(fc0.mean))])
  figlab(lab0,xpos=0.35,ypos=0.99)
  dev.copy2pdf(file=paste(rname,'.int.map.pdf',sep=''))

  FC0 <- aggregate(attr(Z,'sf4bmd'),year,'mean')
  FC0 <- subset(FC0,it=length(index(FC0)))
  map(subset(FC0,is=list(lon=range(lon(fc.mean)),lat=range(lat(fc.mean)))),colbar=colbar)
  lab0 <- paste('ECMWF original',paste(paste(attr(Z,'mons'),collapse='-')),
                'starting in',year(FC0)[length(index(FC0))])
  figlab(lab0,xpos=0.35,ypos=0.99)
  dev.copy2pdf(file=paste(rname,'.original.map.pdf',sep='')) 
}
## Plot the time series:

dev.new()
display.sf(sfc.mos,sfc.int,is=is)
dev.copy2pdf(file=paste(rname,'.hcst.pdf',sep=''))
  
## Distribution
dev.new()

## Use the observations to define upper, neutral and lower categories
Y3m <- obsthreemonths(subset(tx.bmd,it=c(1981,2010)),
                             attr(sfc.mos,'mons'),FUN='mean')
distribution.sf(sfc.mos,sfc.int,Y3m,is=is,col=c('blue','green','red'))
dev.copy2pdf(file=paste(rname,'.tercile.pdf',sep=''))
}

## End of the section with functions -------------------------------


##-------------------------------------------------------------------
## Carry out the analysis
## Retrieve the predictands: station data -> PCA

if (!file.exists('bmd.rda')) {
  ## Read the stations data from BMD and convert to esd station
  ## and save in an R-binary file to speed up ext analyses
  pr.bmd <- readBMD()
  ## Some of the temperature records have problems - wrong dates!
  tx.bmd <- readBMD('~/Dropbox/BMD/tmax',exclude='11927')
  tn.bmd <- readBMD('~/Dropbox/BMD/tmin',exclude=c('11927','41977'))
  save(file='bmd.rda',pr.bmd,tx.bmd,tn.bmd)
} else load('bmd.rda')

## Make maps to examine the station data
#gridmap(annual(subset(pr.bmd,it=c(1990,2015)),FUN='sum'))
#gridmap(annual(subset(tx.bmd,it=c(1990,2015))))
#gridmap(annual(subset(tn.bmd,it=c(1990,2015))))

## Temperature
print('MOS for temperature')
t2m <- dailymeanT(tx.bmd,tn.bmd)
sf(t2m,is=1)

print('MOS for precipitation')
sf(pr.bmd,param='MSL',FUN='sum',is=1)
