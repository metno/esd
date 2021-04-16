## CMIP6 data
## @RasmusBenestad, 2021-04-16
## This script reads CMIP6 data from Earth System Grid Federation (ESGF)
## this is used to download the ilfes that are too big for R and cause a problem
## The solution is to read in smaller chunks and then use CDO to combine them.
library(esd)
## Get the latest metadata from ESGF

overwrite <- FALSE   
expids <- c("ssp370","ssp126","ssp585")
params <- c("tas","slp","pr")
lon <- NULL; lat <- NULL
localpath <- '/lustre/storeB/project/CMIP/CMIP6.monthly/'
it1 <- c(2015,2040)
it2 <- c(2041,2070)
it3 <- c(2071,2100)

## Skip simulations that cause problems:
#skip <- c("tas_Amon_BCC-ESM1_ssp370_r2i1p1f1.nc","tas_Amon_BCC-ESM1_ssp370_r3i1p1f1.nc")
#skip <- c("tas_Amon_EC-Earth3_ssp370_r131i1p1f1.nc","tas_Amon_EC-Earth3_ssp370_r132i1p1f1.nc","tas_Amon_EC-Earth3_ssp370_r133i1p1f1.nc",
#          "tas_Amon_EC-Earth3_ssp370_r134i1p1f1.nc")
skip <- c()

## Loop thought the SSPs
for (expid in expids) { 
  ## Set the path for a file structure similar to the KNMI ClimateExplorer that works with DSensemble 
  path <- paste0(localpath,expid)
  if (!dir.exists(path)) dir.create(path)
  print(path)
  
  ## Loop through the parameters
  for (param in params) {
    metafname <- paste('meta',expid,param,as.Date(Sys.time()),'rda',sep='.')
    ## Get the ESGF metadata on the CMIP6 runs
    if (!file.exists(metafname)) {
      meta <- meta.ESGF(param=param,expid = expid)
      save(meta,file=metafname)
    } else {
      print(paste('Found',metafname))
      load(metafname)
    }
    yr1 <- as.numeric(substr(meta$period,1,4))
    yr2 <- as.numeric(substr(meta$period,8,11))
    ii1 <- (yr1 >= it1[1]) & (yr2 <= it1[2])
    ii2 <- (yr1 >= it2[1]) & (yr2 <= it2[2])
    ii3 <- (yr1 >= it3[1]) & (yr2 <= it3[2])
    meta1 <- meta[ii1,]
    meta2 <- meta[ii2,]
    meta3 <- meta[ii3,]
    print(table(meta$period))
    
    if (!is.null(meta)) { 
      print(meta$title)
      
      ## Extract the different simulations - there may be repeated ones in several chuncks
      simulations <- rownames(table(paste0(paste(param,'Amon',meta$model,expid,meta$member.id,sep='_'),'.nc')))
      n <- length(simulations)
      print(paste('Download',n,'model simulations for',expid,'and',param))
      downloaded <- list.files(path=path,pattern=paste0(param,'_Amon'))
      if (length(downloaded)==0) overwrite <- TRUE ## If the catalog is new, then download all
      for (i in 1:n) {
        ## Check if the file has already been downloaded - it it is, skip the download unless overwrite option 
        if ( ((length(grep(simulations[i],downloaded))==0) & (length(grep(simulations[i],skip))==0)) | 
             (overwrite) ) {
          print(c(i,n,simulations[i]))
          ## Read the data from the remote ESGF repository
          X <- retrieve.ESGF(im=i,lon=lon,lat=lat,meta=meta1)
          ## Save the file locally
          write2ncdf4(X,file='tmp1.nc')
          ## Read the data from the remote ESGF repository
          X <- retrieve.ESGF(im=i,lon=lon,lat=lat,meta=meta2)
          ## Save the file locally
          write2ncdf4(X,file='tmp2.nc')
          ## Read the data from the remote ESGF repository
          X <- retrieve.ESGF(im=i,lon=lon,lat=lat,meta=meta3)
          ## Save the file locally
          write2ncdf4(X,file='tmp3.nc')
          ## Combine the chunks: 
          system(paste('cdo mergetime tmp1.nc tmp2.nc tmp3.nc ',file.path(path,simulations[i])))
          ## Clean up the temporary files:
          system('rm tmp?.nc')
        } else (print(paste('Skipping',i,' - already downloaded',simulations[i])))
      }
    } else {
      print(paste('-------- Empty metadata',metafname,'--------'))
    }
  }
}
