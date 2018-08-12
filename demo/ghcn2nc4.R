## Save the GHCN station data as a netCDF4 file.
## Go through the data element-wise to generate a netCDF files for each. 
## Read 50 stations per increment and append to the netCDF file. 

#require(esd)
#source('~/R/esd/R/write2ncdf.R')
library(esd)
nmin <- 75
eles <- c(111,121,601)[2:3]
it <- seq(as.Date('1893-01-01'),as.Date('2017-12-1'),by=1)

for (ele in eles) {
  SS <- select.station(src='ghcnd',nmin=nmin,ele=ele)
  stids <- SS$station_id
  ns <- length(stids)
  print(paste('Put',ns,'stations world-wide in netCDF file')) 
  is <- seq(1,ns,by=50)
  if (max(is) < ns) is <- c(is,ns)
  
  param <- tolower(as.character(ele2param(ele,src='ghcnd')[5]))
  print(param)
  fname <- paste(param,'ghcnd','nc',sep='.')
  print(fname)
  append <- file.exists(fname) 
  for (id in is) {for (id in is) {
    iii <- seq(id,id+99,length=100)
    iii <- iii[iii <= ns]
    print('Read data');print(id)
    
    x <- try(station(SS[iii,]))
    
    if (!inherits(x,'try-error')) {
      print(loc(x))
      
      ## Quality check
      if ( (min(x,na.rm=TRUE) < -999) | (max(x,na.rm=TRUE)>2000) ) {
        print("Detected suspect data")
        print(range(x,na.rm=TRUE))
        xc <- coredata(x); xc[xc < -999] <- NA
        xc[xc > 2000] <- NA; coredata(x) <- as.matrix(xc)
      }
      write2ncdf4(x,fname,it=it,append=append,verbose=FALSE,stid_unlim=TRUE)
      print('added to netCDF file')      
    } else {
      print('Failed to get data from GHCN:')
      print(param)
      print(stids[iii])
    }
  }
  }
}

