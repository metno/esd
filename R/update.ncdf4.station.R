#' update.ncdf4.station
#' 
#' The function adds new days of station data to an existing netCDF with daily station data
#' and then updates the summary statistics. Using this function to update netCDF files with 
#' station data can be an effective solution for when reading large volumes of data from a 
#' database is time-consuming. 
#' 
#' @param x station object
#' @param fname file name of the netCDF file with stations to be updated
#' @param vebose For diagnostics
#' 
#' @aliases update.ncdf4.station
#' @author R.E. Benestad
#' @export update.ncdf4.station
update.ncdf4.station <- function(x, file, verbose=TRUE,torg='1899-12-31') {
  if (verbose) print(paste('update.ncdf4.station: ',match.call()))
  if (!file.exists(file)) return(paste('Cannot find',file))
  if (verbose) print(dim(x))
  
  ## Retrieve the metadata, e.g. the summary statistics of 
  Y <- retrieve.stationsummary(file)
  if ( (esd::is.precip(x) & attr(Y,'variable') != 'precip') |
       (esd::is.T(x) & attr(Y,'variable') != 't2m') ) warning('Check variable of reanalysis and input')
  
  ## Identify common stations from the station numbers.
  is <- is.element(Y$station.id,stid(x))
  js <- is.element(stid(x),Y$station.id)
  x0 <- x
  if (!inherits(x,'station')) class(x) <- c('station','day','zoo')
  #browser()
  #x <- subset(x,is=js)
  x <- x[,js]; x <- attrcp(x0,x); class(x) <- class(x0)
  if (verbose) {
    print(paste('Update',sum(is),'stations of',length(Y$station.id)))
    print(attr(Y,'period')); print(range(index(x)))
  }
  
  ## Update the summary statistics
  exclude <- c("location","longitude","latitude","altitude","country","station.id","first.year")
  sumsta <- names(Y)
  for (i in 1:length(exclude)) sumsta <- sumsta[-grep(exclude[i],sumsta)]
  if (verbose) print(sumsta)
  nold <- Y$number.valid[is]
  nval <- apply(coredata(x),2,'nv')
  
  ## Open the netCDF files for appending along the time axis.
  # https://rdrr.io/cran/ncdf4/man/ncvar_add.html
  ncid <- nc_open(file)
  ## Get the times
  tim <- ncvar_get(ncid,'time'); nt <- length(tim)
  tunit <- ncatt_get(ncid,'time','units')
  if (length(grep('days since',tunit$value))) 
    t <- as.Date(substr(tunit$value,12,21)) + tim else
      if (length(grep('months since',tunit$value))) 
        t <- seq(as.Date(substr(tunit$value,14,23)),max(tim),'1 month') else
          if (length(grep('years since',tunit$value))) 
            t <- seq(as.Date(substr(tunit$value,14,23)),max(tim),'1 year')
  if (verbose) print(range(t))
  
  Y$number.valid[is] <- nold + nval
  Y$last.year[is] <- lastyear(x)
  if (verbose) print('read the old data to complete the statistics')
  missval <- ncid$var[[1]]$missval  
  offset <- ncid$var[[1]]$addOffset
  scale <- ncid$var[[1]]$scaleFact
  torg <- sub('days since ','',tunit$value)
  nc_close(ncid)
  
  ## Retrieve the actual data in the old file: Only worry about those stations that are to be updated
  if (verbose) print('read the data from the oroginal netCDF file.')
  y <- retrieve.station(file)
  ## Need to use similar spatial dimension
  if (verbose) {print(dim(y)); print('Need to pad missing data with NAs')}
  nd <- dim(x)[1]; ns <- length(Y$station.id)
  x2 <- matrix(rep(NA,nd*ns),nd,ns)
  x2[,is] <- coredata(x)
  x2 <- zoo(x2,order.by=index(x)); x2 <- attrcp(y,x2); class(x2) <- class(y)
  x <- x2; rm('x2')
  ## Update the statistics from the old data
  if (verbose) print(paste('dimensions for original data',paste(dim(y),collapse=' x ')))
  ## combine the old and new
  if (verbose) print('Combine non-overlapping periods')
  ## Make sure that the two don't overlap in time: only add part of x that is new.
  it <- !is.element(index(x),index(y))
  if (sum(it)==0) return('No new data')
  #x <- subset(x,it=it)
  x <- x[it,]; x <- attrcp(y,x); class(x) <- class(y)
  ## Combine old and new data 
  if (verbose) {print(paste(sum(it),'days added to old data')); print(dim(y)); print(dim(y))}
  y <- c(zoo(y),zoo(x))
  y <- attrcp(x,y); class(y) <- class(x)
  if (verbose) print(paste('dimensions for combined data',paste(dim(y),collapse=' x ')))
  if (verbose) print('annual stats')
  seasons <- c('DJF','MAM','JJA','SON')
  Y$mean[is] <- (nold*Y$mean[is] + nval*apply(coredata(x),2,'mean',na.rm=TRUE)[is])/(nold + nval)
  Y$max[!is.finite(Y$max)] <- missval
  Y$max[is] <- max(cbind(Y$max,coredata(x))[,is],na.rm=TRUE) #?!

  Y$records[is] <- apply(anomaly(subset(y,is=is)),2,'arec') #?!
  if (!is.precip(x)) {
    if (verbose) print('not precip')
    ## The new mean estimated from weighted old and mean of new obs.
    Y$min[is] <- min(cbind(Y$min,coredata(x))[,is],na.rm=TRUE) #?!
    Y$min[!is.finite(Y$min)] <- missval
    Y$sd[is] <- apply(coredata(anomaly(subset(y,is=is))),2,'sd',na.rm=TRUE)
    Y$trend[is] <- trend.coef(annual(subset(y,is=is))) #?!
    lelr <- lastelementrecord(-subset(y,is=is))
    if (verbose) print('seasonal stats')
    for (sea in seasons) {
      if (verbose) print(sea)
      Y[[paste0('mean_',sea)]][is] <- apply(coredata(subset(y,is=is,it=sea)),2,'mean',na.rm=TRUE)
      Y[[paste0('sd_',sea)]][is] <- apply(coredata(subset(anomaly(y),is=is,it=sea)),2,'sd',na.rm=TRUE) #?!
      Y[[paste0('trend_',sea)]][is] <- trend.coef(annual(subset(y,is=is,it=sea),nmin=75))
    }
    Y$lows[is] <- apply(-anomaly(subset(y,is=is)),2,'arec') #?!
  } else {
    if (verbose) print('precip')
    Y$mean[is] <-  apply(annual(subset(y,is=is),'sum'),2,'mean',na.rm=TRUE)
    Y$wetmean[is] <- apply(annual(subset(y,is=is),'wetmean'),2,'mean',na.rm=TRUE)
    Y$wetfreq[is] <- apply(100*annual(subset(y,is=is),'wetfreq'),2,'mean',na.rm=TRUE)
    Y$trend[is] <- trend.coef(annual(subset(y,is=is),FUN='sum'))
    Y$trend_wetmean[is] <- trend.coef(annual(subset(y,is=is),FUN='wetmean'))
    Y$trend_wetfreq[is] <- trend.coef(annual(subset(y,is=is),FUN='wetfreq'))
    Y$sigma2[is] <- apply(subset(y,is=is),2,'rainvar')
    Y$trend_sigma2[is] <- rainvartrend(subset(y,is=is))
    if (verbose) print('seasonal stats')
    for (sea in seasons) {
      if (verbose) print(sea)
      Y[[paste0('mean_',sea)]][is] <- apply(annual(subset(y,is=is,it=sea),'sum',nmin=90),2,'mean',na.rm=TRUE)
      Y[[paste0('wetmean_',sea)]][is] <- apply(annual(subset(y,is=is,it=sea),'wetmean',nmin=90),2,'mean',na.rm=TRUE) #?!
      Y[[paste0('wetfreq_',sea)]][is] <- apply(100*annual(subset(y,is=is,it=sea),'wetfreq',nmin=90),2,'mean',na.rm=TRUE) #?!
      Y[[paste0('trend_',sea)]][is] <- trend.coef(annual(subset(y,is=is,it=sea),FUN='sum',nmin=90)) #?!
      Y[[paste0('trend_wetmean_',sea)]][is] <- trend.coef(annual(subset(y,is=is,it=sea),FUN='wetmean',nmin=90)) #?!
      Y[[paste0('trend_wetfreq_',sea)]][is] <- trend.coef(annual(subset(y,is=is,it=sea),FUN='wetfreq',nmin=90)) #?!
      if (verbose) print('Sigma2')
      Y[[paste0('sigma2_',sea)]][is] <- apply(subset(y,is=is,it=sea),2,'rainvar')
      Y[[paste0('trend_sigma2_',sea)]][is] <- rainvartrend(subset(y,is=is,it=sea),nmin=90) #?!
    }
    lr <- sapply(y,'lastrains')
    ld <- sapply(y,'lastdry')
    ## Mean wet/dry-spell length
    if (verbose) print('Spell')
    t <- index(y)
    ss <- try(spell(y,1))
    if (inherits(ss,'spell')) { 
      mwsl <- colMeans(subset.station(ss,is=list(param='wet')),na.rm=TRUE)
      mdsl <- colMeans(subset.station(ss,is=list(param='dry')),na.rm=TRUE)
    } else {mwsl <- rep(NA,sum(is)); mdsl <- rep(NA,sum(is))}
    if (verbose) print('finished precip stuff')
  }
  
  ## Check if last element is a record
  lehr <- lastelementrecord(y)
  time <- julian(index(x)) - julian(as.Date(torg)) 
  
  ## Write the updated summary statistics over the old values.
  if (verbose) print('Get netCDF-handles')
  ncid <- nc_open(file, write=TRUE)
  nctid <- ncid$var[["time"]]
  ncvar <- ncid$var[[1]]
  lyrid <- ncid$var[["last"]]
  nvid <- ncid$var[["number"]]
  dimS <- ncid$dim[["stid"]]; ns <- dimS$len
  dimT <- ncid$dim[["time"]]; nt <- dimT$len
  meanid <- ncid$var[["summary_mean"]]
  meanid.djf <- ncid$var[["summary_mean_DJF"]]
  meanid.mam <- ncid$var[["summary_mean_MAM"]]
  meanid.jja <- ncid$var[["summary_mean_JJA"]]
  meanid.son <- ncid$var[["summary_mean_SON"]]
  tdid <- ncid$var[["summary_trend"]]
  tdid.djf <- ncid$var[["summary_trend_DJF"]]
  tdid.mam <- ncid$var[["summary_trend_MAM"]]
  tdid.jja <- ncid$var[["summary_trend_JJA"]]
  tdid.son <- ncid$var[["summary_trend_SON"]]
  maxid <- ncid$var[["summary_max"]]
  minid <- ncid$var[["summary_min"]]
  nhrid <- ncid$var[["summary_records"]]
  lehrid <- ncid$var[["last_element_highest"]]
  
  if (!is.precip(x)) {
    if (verbose) print('Not precipitation')
    sdid <- ncid$var[["summary_sd"]]
    sdid.djf <- ncid$var[["summary_sd_DJF"]]
    sdid.mam <- ncid$var[["summary_sd_MAM"]]
    sdid.jja <- ncid$var[["summary_sd_JJA"]]
    sdid.son <- ncid$var[["summary_sd_SON"]]
    nlrid <- ncid$var[["summary_lows"]]
    lelrid <- ncid$var[["last_element_lowest"]]
  } else if (is.precip(x)) {
    if (verbose) print('Precipitation')
    muid <- ncid$var[["summary_wetmean"]]
    muid.djf <- ncid$var[["summary_wetmean_DJF"]]
    muid.mam <- ncid$var[["summary_wetmean_MAM"]]
    muid.jja <- ncid$var[["summary_wetmean_JJA"]]
    muid.son <- ncid$var[["summary_wetmean_SON"]]
    fwid <- ncid$var[["summary_wetfreq"]]
    fwid.djf <- ncid$var[["summary_wetfreq_DJF"]]
    fwid.mam <- ncid$var[["summary_wetfreq_MAM"]]
    fwid.jja <- ncid$var[["summary_wetfreq_JJA"]]
    fwid.son <- ncid$var[["summary_wetfreq_SON"]]
    tdmuid <- ncid$var[["summary_trend_wetmean"]]
    tdmuid.djf <- ncid$var[["summary_trend_wetmean_DJF"]]
    tdmuid.mam <- ncid$var[["summary_trend_wetmean_MAM"]]
    tdmuid.jja <- ncid$var[["summary_trend_wetmean_JJA"]]
    tdmuid.son <- ncid$var[["summary_trend_wetmean_SON"]]
    tdfwid <- ncid$var[["summary_trend_wetfreq"]]
    tdfwid.djf <- ncid$var[["summary_trend_wetfreq_DJF"]]
    tdfwid.mam <- ncid$var[["summary_trend_wetfreq_MAM"]]
    tdfwid.jja <- ncid$var[["summary_trend_wetfreq_JJA"]]
    tdfwid.son <- ncid$var[["summary_trend_wetfreq_SON"]]
    lrid <- ncid$var[["summary_lastrains"]]
    ldid <- ncid$var[["summary_lastdry"]]
    sigma2id <- ncid$var[["summary_sigma2"]]
    sigma2id.djf <- ncid$var[["summary_sigma2_DJF"]]
    sigma2id.mam <- ncid$var[["summary_sigma2_MAM"]]
    sigma2id.jja <- ncid$var[["summary_sigma2_JJA"]]
    sigma2id.son <- ncid$var[["summary_sigma2_SON"]]
    tsigma2id <- ncid$var[["summary_trend_sigma2"]]
    tsigma2id.djf <- ncid$var[["summary_trend_sigma2_DJF"]]
    tsigma2id.mam <- ncid$var[["summary_trend_sigma2_MAM"]]
    tsigma2id.jja <- ncid$var[["summary_trend_sigma2_JJA"]]
    tsigma2id.son <- ncid$var[["summary_trend_sigma2_SON"]]
    mwslid <- ncid$var[["summary_mean_wetdur"]]
    mdslid <- ncid$var[["summary_mean_drydur"]]
  }
  ## Only overwrite the stations
  # same <- is.element(Y$station.id,stid(x))
  # if (verbose) print(paste(sum(same),'stations with similar IDs'))
  # s1 <- min((1:length(same))[same])           # start
  # c1 <- max((1:length(same))[same]) - s1 + 1  # count
  # if (verbose) print(c(s1,c1))
  s1 <- 1; c1 <- length(Y$station.id)
  
  if (verbose) print('over-write summary statistics')
  ncvar_put( ncid, lyrid, Y$last.year,start=s1,count=c1)
  ncvar_put( ncid, meanid, Y$mean,start=s1,count=c1)
  ncvar_put( ncid, meanid.djf, Y$mean_DJF,start=s1,count=c1)
  ncvar_put( ncid, meanid.mam, Y$mean_MAM,start=s1,count=c1)
  ncvar_put( ncid, meanid.jja, Y$mean_JJA,start=s1,count=c1)
  ncvar_put( ncid, meanid.son, Y$mean_SON,start=s1,count=c1)
  ncvar_put( ncid, tdid, Y$trend,start=s1,count=c1)
  ncvar_put( ncid, tdid.djf, Y$trend_DJF,start=s1,count=c1)
  ncvar_put( ncid, tdid.mam, Y$trend_MAM,start=s1,count=c1)
  ncvar_put( ncid, tdid.jja, Y$trend_JJA,start=s1,count=c1)
  ncvar_put( ncid, tdid.son, Y$trend_SON,start=s1,count=c1)
    ncvar_put( ncid, maxid, Y$max,start=s1,count=c1)
  ncvar_put( ncid, nhrid, Y$records,start=s1,count=c1)
  ncvar_put( ncid, nvid, Y$number.valid,start=s1,count=c1)
  #print('.?.')
    #ncvar_put( ncid, lehrid, Y$lehr,start=s1,count=c1) # last element highest
  if (!is.precip(x)) {
    if (verbose) print('extras for non-precipitation')
    ncvar_put( ncid, minid, Y$min,start=s1,count=c1)
    #print('.'); print(names(Y));
    Y$sd[!is.finite(Y$sd)] <- missval
    #print(Y$sd); str(sdid)
    if (!is.null(sdid)) ncvar_put( ncid, sdid, Y$sd,start=s1,count=c1)
    if (!is.null(sdid)) ncvar_put( ncid, sdid.djf, Y$sd_DJF,start=s1,count=c1)
    if (!is.null(sdid)) ncvar_put( ncid, sdid.mam, Y$sd_MAM,start=s1,count=c1)
    if (!is.null(sdid)) ncvar_put( ncid, sdid.jja, Y$sd_JJA,start=s1,count=c1)
    if (!is.null(sdid)) ncvar_put( ncid, sdid.son, Y$sd_SON,start=s1,count=c1)
    if (!is.null(nlrid)) ncvar_put( ncid, nlrid, Y$lows,start=s1,count=c1)
    #ncvar_put( ncid, lelrid, lelr,start=s1,count=c1)
  }
  if (is.precip(x)) {
    if (verbose) print('extras for precipitation')
    ncvar_put( ncid, muid, Y$wetmean,start=s1,count=c1)
    ncvar_put( ncid, muid.djf, Y$wetmean_DJF,start=s1,count=c1)
    ncvar_put( ncid, muid.mam, Y$wetmean_MAM,start=s1,count=c1)
    ncvar_put( ncid, muid.jja, Y$wetmean_JJA,start=s1,count=c1)
    ncvar_put( ncid, muid.son, Y$wetmean_SON,start=s1,count=c1)
    ncvar_put( ncid, fwid, Y$wetfreq,start=s1,count=c1)
    ncvar_put( ncid, fwid.djf, Y$wetfreq_DJF,start=s1,count=c1)
    ncvar_put( ncid, fwid.mam, Y$wetfreq_MAM,start=s1,count=c1)
    ncvar_put( ncid, fwid.jja, Y$wetfreq_JJA,start=s1,count=c1)
    ncvar_put( ncid, fwid.son, Y$wetfreq_SON,start=s1,count=c1)
    ncvar_put( ncid, tdfwid, Y$trend_wetfreq,start=s1,count=c1)
    ncvar_put( ncid, tdfwid.djf,  Y$trend_wetfreq_DJF,start=s1,count=c1)
    ncvar_put( ncid, tdfwid.mam,  Y$trend_wetfreq_MAM,start=s1,count=c1)
    ncvar_put( ncid, tdfwid.jja,  Y$trend_wetfreq_JJA,start=s1,count=c1)
    ncvar_put( ncid, tdfwid.son,  Y$trend_wetfreq_SON,start=s1,count=c1)
    ncvar_put( ncid, tdmuid, Y$trend_wetmean,start=s1,count=c1)
    ncvar_put( ncid, tdmuid.djf, Y$trend_wetmean_DJF,start=s1,count=c1)
    ncvar_put( ncid, tdmuid.mam, Y$trend_wetmean_MAM,start=s1,count=c1)
    ncvar_put( ncid, tdmuid.jja, Y$trend_wetmean_JJA,start=s1,count=c1)
    ncvar_put( ncid, tdmuid.son, Y$trend_wetmean_SON,start=s1,count=c1)
    ncvar_put( ncid, lrid, lr,start=s1,count=c1)
    ncvar_put( ncid, ldid, ld,start=s1,count=c1)
    ncvar_put( ncid, sigma2id, Y$sigma2,start=s1,count=c1)
    ncvar_put( ncid, sigma2id.djf, Y$sigma2_DJF,start=s1,count=c1)
    ncvar_put( ncid, sigma2id.mam, Y$sigma2_MAM,start=s1,count=c1)
    ncvar_put( ncid, sigma2id.jja, Y$sigma2_JJA,start=s1,count=c1)
    ncvar_put( ncid, sigma2id.son, Y$sigma2_SON,start=s1,count=c1)
    ncvar_put( ncid, tsigma2id, Y$trend_sigma2,start=s1,count=c1)
    ncvar_put( ncid, tsigma2id.djf, Y$trend_sigma2_DJF,start=s1,count=c1)
    ncvar_put( ncid, tsigma2id.mam, Y$trend_sigma2_MAM,start=s1,count=c1)
    ncvar_put( ncid, tsigma2id.jja, Y$trend_sigma2_JJA,start=s1,count=c1)
    ncvar_put( ncid, tsigma2id.son, Y$trend_sigma2_SON,start=s1,count=c1)
    if (verbose) print('Mean spell length')
    ncvar_put( ncid, mwslid, mwsl,start=s1,count=c1)
    ncvar_put( ncid, mdslid, mdsl,start=s1,count=c1)
  } 
  
  ## Add the new data to the old
  if (verbose) print('write new data...')
  n <- length(index(x))
  start <- c(1,nt+1)
  count <- c(ns,n)
  if (verbose) print(rbind(start,count))
  # X <- matrix(rep(NA,ns*n),n,ns)
  # X[,is] <- coredata(subset(x,is=js))
  if (verbose) print(paste('add new times',index(x)[1],'-',index(x)[length(index(x))],
                           'to get',index(y)[1],'-',index(y)[length(index(y))]))
  ##https://stackoverflow.com/questions/30084261/extend-dimensions-in-netcdf-file-using-r
  #ncvar_put( ncid, varid='time', index(x),start=start[2],count=count[2])
  ncvar_put( ncid, varid='time', index(y) - julian(as.Date(torg)),start=1,count=length(index(y)))
  if (verbose) {print('add new data'); print(dim(y))}
  #ncvar_put( ncid, ncvar, t(X),start=start,count=count)
  y[!is.finite(y)] <- missval
  y <- round((y - offset)/scale)
  ncvar_put( ncid, ncvar, t(coredata(y)))
  nc_close(ncid)
  if (verbose) print('success :-)')
}

## Test function
#' @export test.update.ncdf4.station
test.update.ncdf4.station <- function(param='t2m',l=1,is=1:20,verbose=FALSE,plot=TRUE) {
  print('test: update.ncdf4.station')
  if (file.exists('test0.nc')) file.remove('test0.nc')
  if (file.exists('test1.nc')) file.remove('test1.nc')
  x0 <- station.thredds(param=param,is=is)
  x1 <- subset(x0,it=range(year(x0))-c(0,l))
  print('test file: test0.nc'); print(range(index(x0))); print(dim(x0))
  T1 <- Sys.time()
  write2ncdf4.station(x0,file='test0.nc')
  T2 <- Sys.time()
  print(paste('Time it took to save original netCDF-file with all data was',round(T2-T1),'s'))
  print('test file: test1.nc'); print(range(index(x1)))
  write2ncdf4.station(x1,file='test1.nc')
  x2 <- subset(x0,it=max(year(x0))-c(l,0))
  print('update.ncdf4.station')
  print(dim(x2))
  t1 <- Sys.time()
  update.ncdf4.station(x2,file='test1.nc',verbose=verbose)
  t2 <- Sys.time()
  print(paste('The time taken to update file was',round(t2-t1),'s'))
  rm('x2')
  ## Read the updated file
  print('--- Check the files ---')
  x0 <- retrieve.station('test0.nc')
  print('...')
  x1 <- retrieve.station('test1.nc')
  print('Diagnostics...')
  print(range(index(x0))); print(range(index(x1)))
  print(dim(x0)); print(dim(x1))
  print(paste(mean(zoo(x0) - zoo(x1),na.rm=TRUE),'mean difference'))
  if (plot) {
    plot(zoo(x0) - zoo(x1))
    par(mfcol=c(3,1))
    image(index(x0),is,coredata(x0),main='Original file')
    image(index(x1),is,coredata(x1),main='Updated file')
    image(index(x0),is,coredata(x1)- coredata(x0),main='Difference')
  }
  print('Compare data dimensions:') 
  print(dim(x0)); print(dim(x1))
  Y0 <- retrieve.stationsummary('test0.nc')
  Y1 <- retrieve.stationsummary('test1.nc')
  print('Compare metadata+summary - differences:')
  print(setdiff(names(Y0),names(Y1)))
  for (ii in names(Y0)) {if (!is.null(Y0[[ii]]) & !is.null(Y1[[ii]]) & is.numeric(Y0[[ii]])) 
    print(paste(ii,sum(Y0[[ii]]-Y1[[ii]] > 0, na.rm=TRUE))) else print(paste(ii,setdiff(Y0[[ii]],Y1[[ii]])))}
  file.remove('test0.nc'); file.remove('test1.nc');
  print(attr(Y0,'period'))
  print(attr(Y1,'period'))
  print('Finished test')
}

