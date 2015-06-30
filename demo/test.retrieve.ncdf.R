## R-script made to test the reading of GCMs in the ensemble.

library(esd)

retrieve <- function(x,...) {
  require(ncdf)
  y <- retrieve.ncdf(x,...)
  return(y)
}

path <- 'CMIP5.monthly/'; rcp <- 'rcp45'
verbose <- FALSE
lon <- c(0,30)
lat <- c(50,70)
predictor="ERA40_t2m_mon.nc"
FUNX="mean"
pattern="tas_Amon_ens_"

t2m <- retrieve(ncfile=predictor,lon=lon,lat=lat,verbose=verbose) 

  #rm("predictor"); gc(reset=TRUE)
  #t2m <- t2m.ERA40(lon=lon,lat=lat)
t2m <- subset(as.4seasons(t2m,FUN=FUNX),it='djf')
T2m <- as.anomaly(as.annual(aggregate.area(t2m,FUN=FUNX),FUN=FUNX))

##T2M <- as.4seasons(t2m,FUN=FUNX)
  # Fix - there is a bug with 'T2M <- as.4seasons(t2m,FUN=FUNX)' - the date is not correct
##DJF <- subset(as.4seasons(t2m,FUN=FUNX),it='djf')
##MAM <- subset(as.4seasons(t2m,FUN=FUNX),it='mam')
#JJA <- subset(as.4seasons(t2m,FUN=FUNX),it='jja')
#SON <- subset(as.4seasons(t2m,FUN=FUNX),it='son')

##plot(MAM,xlim=as.Date(c('1900-01-01','2100-12-31')),ylim=c(-5,10))
plot(T2m,xlim=c(1900,2100),ylim=c(-4,6))

## browser()
## Ensemble GCMs
path <- file.path(path,rcp,fsep = .Platform$file.sep)
ncfiles <- list.files(path=path,pattern=pattern,full.name=TRUE)
N <- length(ncfiles)

if (verbose) print("loop...") 
  for (i in 1:N) {
    print(ncfiles[i])
    gcm <- retrieve(ncfile = ncfiles[i],
                          lon=lon,
                          lat=lat,verbose=verbose)
    print(attr(gcm,'model_id'))
    ## REB: 30.04.2014 - new lines...
    browser()
    print('...')
    gcm <- subset(as.4seasons(gcm,FUN=FUNX),it='djf')
    ##MAMGCM <- subset(as.4seasons(gcm,FUN=FUNX),it='mam')
    ##JJAGCM <- subset(as.4seasons(gcm,FUN=FUNX),it='jja')
    ##SONGCM <- subset(as.4seasons(gcm,FUN=FUNX),it='son')
    ##y  <- as.anomaly(as.annual(aggregate.area(gcm,FUN=FUNX),FUN=FUNX))
    y <- aggregate.area(as.anomaly(gcm),FUN='mean')
    if (diff(range(y,na.rm=TRUE)) > 10) {
        print(paste(attr(gcm,'model_id'),'. Max T(2m) diff.',
                    diff(range(y,na.rm=TRUE)),
                    'min time step',min(diff(index(gcm))),
                    'max time step',max(diff(index(gcm)))))
      print(attr(gcm,'calendar'))
  }
    
    lines(y,col=rgb(0.5,0.5,0.5,0.3))
}
lines(T2m,xlim=c(1900,2100),ylim=c(-4,6))
