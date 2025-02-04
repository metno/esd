library(esd)

setwd('~/R')
scen <- 'ssp245'

path <- paste0('~/data/CMIP/CMIP6.monthly/',scen,'/')
ssps <- list.files(path=path,pattern='tas_Amon')
ssps <- ssps[grep('1850-2100',ssps)]

results <- list()
fname <- paste0(scen,'.mon.gmst.rda')
if (file.exists(fname)) {
  load(fname)
  print(names(results))
  ssps <- ssps[!is.element(ssps,names(results))]
}
ignore <- c('HR_ssp245_')
if (length(ignore)>0) ssps <- ssps[-grep(ignore,ssps)]

bad <- c()
for (gcm in ssps) {
  print(gcm)
  # x <- try(retrieve(file.path(path,gcm)))
  # if (!inherits(x,'try-error')) { 
  # y <- aggregate.area(x,FUN='mean')
  system(paste('cdo fldmean',file.path(path,gcm),'gmst.nc'))
  ncid <- nc_open('gmst.nc')
  x <- ncvar_get(ncid,varid='tas') - 273.15
  t <- ncvar_get(ncid,varid='time')
  tunit <- ncatt_get(ncid,varid='time','units')$value
  origin <- sub(" 00:00:00.0","",sub(".*since ", "", tunit))
  if (length(grep('hours',tunit))) t <- t/24
  nc_close(ncid)
  y <- zoo(x,order.by=as.Date(t,origin=origin))
  nc_close(ncid)
  results[[gcm]] <- y
  plot(zoo(y),main=gcm)
  save(results,file=fname)
  #} else {print(paste(gcm,'is bad')); bad <- c(bad,file.path(path,gcm))}
  rm('x'); gc(reset=TRUE)
}

