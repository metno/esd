## Read gridded empirical-statistical downscaling (ESD) data from thredds
## This demo shows how ESD results derived from different regions and gridded to a 5-arc-minute resolution
## with ETOPO5 as a covariate in a kriging can be accessed from the OpenDap/thredds at Met Norway.
## The results consist of the ensemble mean, ensemble standard deviation, as well as three different quality measures:
## (1) if the ensemble distribution is approximately normal, (2) if the spread of trends within the ensemble members is 
## consistent with corresponding observed trend, and (3) if the observed interannual variability is consistent with the 
## ensemble spread. These results are described in more details in https://doi.org/10.5194/hess-29-45-2025 and 
## https://doi.org/10.1073/pnas.2503806122, in addition to other papers in progress. Information about  these results
## are  provided in the netCDF headers, and the attributes 'predictor' and 'predictand' provide details about the 
## calibration and the CMIP models involved. 
## Rasmus Benestad, April 22, 2026

library(esd)

ESDresults <- function(param='tas',region='Nordics',ssp='ssp370',verbose=TRUE,...,
                    urlpath='https://thredds.met.no/thredds/dodsC/metusers/rasmusb/'){
  fname <- switch(region,
                  'Nordics'=paste(param,'DSEns_Nordics_4seasons_1950-2100',ssp,sep='_'),
                  'Barents'=paste(param,'DSEns_BarentsSea',ssp,'4seasons_1950-2100',sep='_'),
                  'southeastAfrica'=paste(param,'Ayear_DSEns_CORDEXFPS_southeast_Africa_1940-2099',ssp,sep='_'),
                  'Ghana'=paste(param,'DSEns_Ghana',ssp,'4seasons_1950-2100',sep='_'))
  fname <- file.path(urlpath,paste(fname,'nc',sep='.'))
  if (verbose) cat('thredds: ',fname,'\n')
  params <- ncvars(fname)
  X <- list()
  for (param in params) X[[param]] <- retrieve(param=param,fname,...)
  return(X)
}

x <- ESDresults()
## Show the fist element/parameter in the list
map(x[[1]])
## Show quality on trends
map(x[['z_score_trend']])

## The Barents Sea region
y <- ESDresults(region='Barents')
map(y[[1]])

## The CORDEX FPS southeast Africa: total annual rainfall
z <- ESDresults(param='pr',region='southeastAfrica')
map(z[[1]])

## Results from springs: the wet-day frequency
w <- ESDresults(param='fw',region='Ghana')
map(w[[1]])