## A script that saves the gridded data from downscaling of ensembles as netCDF files
## This script takes results from the R-script cdsehistssp.R which have already been
## converted to station objects.
library(esd)
verbose <- FALSE
path <- './'
output.path <- paste(path,'netCDF_files',sep='/')
region <- 'Europe'             ## Region marker for the final netCDf files
## Move to the catalogue for the project 
setwd(path)
dsepath <- './'
pattern <- 'dse.ecad.'  ## Text pattern of data files to work with
textmarker <- 'ecad.'              ## The text pattern before the part indicating SSP/CRP scenario
dse.results <- list.files(path=dsepath,pattern=pattern,full.names = TRUE)

## This function extracts the interval where all members of the ensemble are present. It differs from 'allN'
## as it works on the PCA-based ensemble results as opposed to single stations. 
allMs <- function(x,verbose=FALSE) {
  if (verbose) cat('allMs \n')
  nt <- length(index(x))
  cls <- class(x)
  ## Remove members with short time span
  nv <- unlist(lapply(x,function(x) if (length(dim(x))==2) y <-nv(x[,1]) else  
    y <- NA))[-c(1,2,length(x))]
  im <- nv >= quantile(nv,probs = 0.75)
  if (verbose) cat('keep',sum(im),' members of',length(im),' \n')
  x <- subset(x,im=im)
  ## Select the time interval where the ensemble is complete.
  if (verbose) cat('select time interval \n')
  it1 <- unlist(lapply(x,function(x) min(index(x[is.finite(x)]))))[-c(1,2,length(x))]
  it2 <- unlist(lapply(x,function(x) max(index(x[is.finite(x)]))))[-c(1,2,length(x))]
  it <- c(max(it1,na.rm=TRUE),min(it2,na.rm=TRUE))
  if (verbose) cat(it,'\n')
  x <- subset(x,it=it)
  class(x) <- cls
  if (verbose) cat('finished allMs')
  return(x)
}

## Function that returns Shapiro-Wilk normality test for each year 
ensembletest <- function(x) {
  ## Assume class(dse[[1]])
  ## [1] "dsensemble" "zoo"
  d <- dim(x)
  y <- rep(NA,d[1]); p <- y
  for (it in 1:d[1]) {
    z <- coredata(x)[it,]
    z <- z[is.finite(z)]
    #print(summary(z))
    y[it] <- shapiro.test(z)$statistic
    p[it] <- shapiro.test(z)$p.value
  }
  #print(summary(y))
  Y <- mean(y,na.rm=TRUE)
  attr(Y,'p.value') <- mean(p,na.rm=TRUE)
  return(Y)
}

ranktest <- function(x,plot=FALSE) {
  #print(dim(x))
  y <- subset(zoo(attr(x,'station')),it=range(index(x)))
  index(y) <- year(y)
  x <- subset(x,it=range(year(y)))
  z <- merge(y,zoo(x),all=FALSE)
  r <- apply(z,1,function(x) order(x)[1])/dim(z)[2]
  if (plot) {
    plot(z,plot.type='single',col='grey',main='Rank test')
    lines(y,type='b')
    text(index(y),y,r,sub=1)
    hist(r)
  }
  #print(c(length(y),length(r),range(r)))
  rank.test <- ks.test(r,'punif')$statistic
  #attr(rank.test,'p.value') <- ks.test(r,'punif')$p.value
  return(rank.test)
}

trendtest <- function(x,plot=FALSE) {
  #print(dim(x))
  y <- subset(zoo(attr(x,'station')),it=range(index(x)))
  index(y) <- year(y)
  x <- subset(x,it=range(year(y)))
  z <- merge(y,zoo(x),all=FALSE)
  t <- coredata(apply(z,1,trend.coef))
  z_score <- (t[1] - mean(t[-1])) / sd(t[-1])
  if (plot) {hist(t[-1]); lines(rep(t[1],2),c(0,100),lty=2,lwd=2)}
  #print(z_score)
  return(z_score)
}

if (!file.exists(output.path)) dir.create(output.path)

for (one in dse.results) {
  load(one)
  cat(one,' \n')
  dse.pca <- Z; rm('Z')
  dse <- as.station(dse.pca)
  locs <- unlist(lapply(dse,function(z) loc(z)))
  lons <- unlist(lapply(dse,function(z) lon(z)))
  lats <- unlist(lapply(dse,function(z) lat(z)))
  alts <- unlist(lapply(dse,function(z) alt(z)))
  cntrs <- cntr(dse.pca$pca)
  ## Test whether the ensemble spread is normally distributed for each year and site
  cat('apply quality tests \n')
  ensdist <- unlist(lapply(dse,ensembletest))
  rankscores <- unlist(lapply(dse,ranktest))
  trendscores <- unlist(lapply(dse,trendtest))
  
  ## Extract the ensemble mean for each site (list element)
  z_mean <- lapply(dse, function(z) zoo(rowMeans(coredata(z), na.rm = TRUE), index(z)))
  z_sd <- lapply(dse, function(z) zoo(apply(coredata(z),1,sd,na.rm = TRUE),index(z)))
  ## Extract required metadata
  em.ssp <- do.call(merge, c(z_mean, all = TRUE))
  es.ssp <- do.call(merge, c(z_sd, all = TRUE))
  ## Make a station object of the combined results
  cat('Make station objects \n')
  ## Ensemble mean
  longname <- attr(dse.pca$pca,'longname')
  param <- varid(attr(dse[[1]],'station'))
  #if (length(grep(param,one))==0) browser('Check param')
  unit <- unit(attr(dse[[1]],'station'))
  ssp <- sub("\\..*","",sub(paste0(".*",textmarker), "", one))
  cat('ssp= ',ssp,'\n')
  
  em.ssp <- as.station(em.ssp,param=param,unit=unit,
                       location=loc(dse),longname =longname,
                       lon= lons, lat=lats, alt=alts, cntr = cntrs)
  ## Ensemble spread
  es.ssp <- as.station(es.ssp,param=param,unit=unit,
                       location=loc(dse),longname =longname,
                       lon= lons, lat=lats, alt=alts, cntr = cntrs)
  ## Grid ensemble spread information
  attr(ensdist,'longitude') <- c(lons)
  attr(ensdist,'latitude') <- c(lats)
  attr(ensdist,'altitude') <- c(alts)
  attr(ensdist,'location') <- locs
  attr(ensdist,'variable') <- 'shapiro.test'
  attr(ensdist,'unit') <- 'statistic'
  Ensdist <- gridmap(ensdist,verbose=verbose)
  attr(Ensdist,'history') <- history.stamp()
  attr(Ensdist,'call') <- 'dse2ncdf.R' 
  map(Ensdist,main='Ensemble distribution normal?')
  points(lon(ensdist),lat(ensdist),cex=0.5)
  attr(Ensdist,'longname') <- 'Test for normal distribution in ensemlbe spread'
  yrs <- paste(range(year(attr(dse[[1]],'station'))),collapse='-')
  normfile <- paste(param,'DSEns_Normal-test',region,yrs,ssp,'.nc',sep='_')
  write2ncdf4(Ensdist,file=normfile)
  
  ## Grid ensemble trend test
  attr(trendscores,'longitude') <- c(lons)
  attr(trendscores,'latitude') <- c(lats)
  attr(trendscores,'altitude') <- c(alts)
  attr(trendscores,'location') <- locs
  attr(trendscores,'variable') <- 'z_score'
  attr(trendscores,'unit') <- 'statistic'
  Trendscores <- gridmap(trendscores,verbose=verbose)
  attr(Trendscores,'history') <- history.stamp()
  attr(Trendscores,'call') <- 'dse2ncdf.R' 
  map(Trendscores,main='Observed consistent with ensemble trends?')
  points(lon(trendscores),lat(trendscores),cex=0.5)
  attr(Trendscores,'longname') <- 'z_score <- (obs_trend - mean(ensemble_trend)) / sd(ensemble_trend)'
  trendfile <- paste(param,'DSEns_trend-test',region,yrs,ssp,'.nc',sep='_')
  write2ncdf4(Trendscores,file=trendfile)
  
  ## Grid ensemble rank test
  attr(rankscores,'longitude') <- c(lons)
  attr(rankscores,'latitude') <- c(lats)
  attr(rankscores,'altitude') <- c(alts)
  attr(rankscores,'location') <- locs
  attr(rankscores,'variable') <- 'ks.test_uniform'
  attr(rankscores,'unit') <- 'statistics'
  Rankscores <- gridmap(rankscores,verbose=verbose)
  attr(Rankscores,'history') <- history.stamp()
  attr(Rankscores,'call') <- 'dse2ncdf.R' 
  map(Rankscores,main='Observed consistent with ensemble spread')
  points(lon(rankscores),lat(rankscores),cex=0.5)
  attr(Rankscores,'longname') <- 'Rank-test between ensemble members and observations'
  rankfile <- paste(param,'DSEns_rank-test',region,yrs,ssp,'.nc',sep='_')
  write2ncdf4(Rankscores,file=rankfile)
  
  ## Estimate PCA 
  cat('compute PCAs \n')
  pca.em.ssp <- PCA(pcafill(em.ssp),n=5)
  pca.es.ssp <- PCA(pcafill(es.ssp),n=5)
  plot(pca.em.ssp,new=FALSE,main=paste('Ensemble mean',param,ssp))
  plot(pca.es.ssp,new=FALSE,main=paste('Ensemble spread',param,ssp))
  ## Grid the PCA and transform to EOF
  tmp.dse.eof <- paste('eof.dse',param,region,ssp,'rda',sep='.')
  cat('compute EOFs \n')
  if (!file.exists(tmp.dse.eof)) {
    eof.em.ssp <- gridmap(pca.em.ssp)
    attr(eof.em.ssp,'variable') <- attr(eof.em.ssp,'variable')[1]
    attr(eof.em.ssp,'unit') <- unit[1]
    eof.es.ssp <- gridmap(pca.es.ssp)
    attr(eof.es.ssp,'variable') <- attr(eof.es.ssp,'variable')[1]
    attr(eof.es.ssp,'unit') <- attr(eof.es.ssp,'unit')[1]
    save(eof.em.ssp,eof.es.ssp,file=tmp.dse.eof)
  } else load(tmp.dse.eof)
  plot(eof.em.ssp,new=FALSE)
  plot(eof.es.ssp,new=FALSE)
  ## Translate the EOFs into gridded field data (one lon-lat grid for each year)
  cat ('EOF -> field \n')
  zm <- as.field(eof.em.ssp,anomaly=FALSE,verbose=verbose)
  zs <- as.field(eof.es.ssp,anomaly=FALSE,verbose=verbose)
  yrs <- paste(range(year(zm)),collapse='-')
  mz <- map(zm,main=paste('Simulated annual ',param),
            colbar=list(pal='precip.ipcc'),new=FALSE)
  points(lons,lats,cex=0.5)
  
  tz <- map(zm,FUN='trend',
            main=paste('Simulated future trends in annual',param),
            colbar=list(breaks=seq(-50,50,by=10),pal='precip.ipcc'),new=FALSE)
  points(lons,lats,cex=0.5)
  
  s1 <- map(zs,main='SSP3-70 Ensemble spread 1991-2020',new=FALSE)
  points(lons,lats,cex=0.5)
  
  ncfile.em <- paste(param,'Ayear_DSEnsMean',region,yrs,ssp,'.nc',sep='_')
  ncfile.es <- paste(param,'Ayear_DSEnsSpread',region,yrs,ssp,'.nc',sep='_')
  
  cat(ncfile.em,'\n')
  ## Save to netCDF file
  attr(zm,'variable') <- paste0('ens_mean_',param)
  attr(zm,'unit') <- attr(dse.pca$pca,'variable')[1]
  attr(zm,'longname') <- attr(dse.pca$pca,'longname')
  attr(zm,'source') <- paste(paste(paste(locs,cntrs,lons,lats,sep=', '),collapse='; '),
                             'GCM runs:',paste(attr(dse,'model_id'),collapse='; '))
  attr(zm,'type') <- 'climate statistics; threshold of 1 mm/day'
  attr(zm,'aspect') <- 'aggregated over hydrological Jan-Dec'
  attr(zm,'description') <- paste('Ensemble of',dim(dse[[1]])[2],
                                  'downscaled multi-model CMIP6 runs',ssp,
                                  'calibrated on in-situ rain gauge data and ERA5',
                                  'and subject to kriging')
  attr(zm,'info') <- 'DOI:10.5194/hess-29-45-2025'
  attr(zm,'greenwich') <- TRUE
  attr(zm,'frequency') <- 'annual'
  attr(zm,'institution') <- 'Met Norway'
  attr(zm,'author') <- 'R.E. Benestad'
  attr(zm,'history') <- history.stamp()
  attr(zm,'call') <- 'dse2ncdf.R' 
  write2ncdf4(zm,file=ncfile.em)
  cat(ncfile.es,'\n')
  zs <- attrcp(zm,zs)
  attr(zs,'variable') <- paste0('ens_sd_',param)
  attr(zs,'longname') <- 'Ensemble spread total annual rainfall (stdv)'
  write2ncdf4(zs,file=ncfile.es)
  
  ncfile <- file.path(output.path,paste0(param,'_Ens_',region,'_',
                      sub('\\.','_',sub('.rda','.nc',sub(paste0(".*",textmarker), "", one)))))
  cat(paste('cdo merge',ncfile.em,ncfile.es,normfile,trendfile,rankfile,ncfile),'\n')
  system(paste('cdo -O merge',ncfile.em,ncfile.es,normfile,trendfile,rankfile,ncfile))
  cat('-------------------- \n')
}

## Also: 
## (1) how close ids the ensemble spread to being normal 
## Function: shapiro.test(x)
## Hypothesis: H0 is that the data is normally distributed. A p<0.05 indicates 
## you should reject the null hypothesis (data is likely not normal).
## (2) an evaluation of the downscaled results.
## Use the rank of observations with the ensemble members
## ks.test(x, "punif", min = min_val, max = max_val)
## (3) an evaluation of trends
## z_score <- (obs_trend - mean(ensemble_trend)) / sd(ensemble_trend)

## Combine the different results and scores into one single netCDF file with CDO
# cdo merge file1.nc file2.nc file3.nc output_merged.nc
#cluster sesongdata for å dele opp i klimaregioner
#R Implementation: Use pam() from the cluster package.
#så, finn PCA for regionale prediktander og nedskalere med EOFer som dekker omtrent samme område.

