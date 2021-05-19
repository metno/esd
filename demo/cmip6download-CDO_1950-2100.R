## CMIP6 data
## @RasmusBenestad, 2021-04-16
## This script reads CMIP6 data from Earth System Grid Federation (ESGF)
## this is used to download the ilfes that are too big for R and cause a problem
## The solution is to read in smaller chunks and then use CDO to combine them.
library(esd)
## Get the latest metadata from ESGF

## Function to check the time period covered
getinterval <- function(file,param='tas',verbose=FALSE) {
  if (!file.exists(file)) return(NULL)
  ncid <- nc_open(file)
  time <- check.ncdf4(ncid,param=param,verbose=verbose)$time 
  nc_close(ncid)
  print(range(time$vdate))
  return(time$vdate)
}

## Function to split the CMIP file names into model, ssp and RIPF
decipher <- function(x) {
  i <- gregexpr('_',x)[[1]]
  model <- substr(x,i[2]+1,i[3]-1)
  ssp <- substr(x,i[3]+1,i[4]-1)
  ripf <- substr(x,i[4]+1,i[4]+8)
  return(c(model,ssp,ripf))
}


overwrite <- FALSE   
expids <- c("ssp370","ssp126","ssp585","ssp245","ssp119")
params <- c("tas","psl","pr")
lon <- NULL; lat <- NULL
localpath <- '/lustre/storeB/project/CMIP/CMIP6.monthly/'

## Skip simulations that cause problems:
#skip <- c("tas_Amon_BCC-ESM1_ssp370_r2i1p1f1.nc","tas_Amon_BCC-ESM1_ssp370_r3i1p1f1.nc")
#skip <- c("tas_Amon_EC-Earth3_ssp370_r131i1p1f1.nc","tas_Amon_EC-Earth3_ssp370_r132i1p1f1.nc","tas_Amon_EC-Earth3_ssp370_r133i1p1f1.nc",
#          "tas_Amon_EC-Earth3_ssp370_r134i1p1f1.nc")
skip <- c()
system('rm tmp.cmip6.*.nc') 
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
    M <- dim(meta)[1]
    print(table(meta$period))
    
    if (!is.null(meta)) { 
      print(meta$title)
      
      ## Extract the different simulations - there may be repeated ones in several chunks
      simulations <- rownames(table(paste(param,'Amon',meta$model,expid,meta$member.id,sep='_')))
      MODS <- rownames(table(paste(meta$model,expid,meta$member.id,sep='_')))
      models <- paste(meta$model,expid,meta$member.id,sep='_')
      n <- length(MODS)
      print(paste('Download',n,'model simulations for',expid,'and',param))
      downloaded <- list.files(path=path,pattern=paste0(param,'_Amon'))
      if (length(downloaded)==0) overwrite <- TRUE ## If the catalog is new, then download all
      for (i in 1:n) {
        print(MODS[i])
        imatch <- regexpr(simulations[i],downloaded)
        ## Check if the file has already been downloaded - it it is, skip the download unless overwrite option 
        if ( (sum(imatch > 0)==0) |  (overwrite) ) {
          i1 <- regexpr(MODS[i],models)
          ims <- (1:length(models))[i1 > 0]
          print(meta$title[ims])
          print(c(i,n,simulations[i])); print(decipher(simulations[i])); print(ims)
          trueinterval <- c()
          for (ii in 1:length(ims)) { 
            print(meta$model[ims[ii]])
            ## Download the netCDF files for the same model and ripf from the remote ESGF repository
            tmpfile <- paste0('tmp.cmip6.',ii,'.nc')
            test <- try(download.file(meta$http[ims[ii]],tmpfile))
            if (tolower(substr(test,1,5))=="error") system('rm tmp.cmip6.*.nc') else
              yrs <- try(year(getinterval(tmpfile,param)))
              if (!inherits(yrs,'try-error'))
                trueinterval <- try(range(year(c(trueinterval,getinterval(tmpfile,param)))))
          }
          if (is.null(trueinterval)) trueinterval <- c(NA,NA)
          if (tolower(substr(test,1,5))=="error") print('Skip this file') else { 
            ncname <- paste0(paste(param,'Amon',MODS[i],paste(trueinterval,collapse='-'),sep='_'),'.nc')
            print(ncname)
            ## Combine the chunks: 
            if (length(ims)>1) system(paste('cdo -O mergetime tmp.cmip6.*.nc ',file.path(path,ncname))) else
              if (length(ims)==1) file.rename(tmpfile,file.path(path,ncname))
            ## Clean up the temporary files:
            system('rm tmp.cmip6.*.nc') 
          }
        } else (print(paste('Skipping',i,' - already downloaded',simulations[i])))
      }
    } else {
      print(paste('-------- Empty metadata',metafname,'--------'))
    }
  }
}
