## Save the ECA&D station data as a netCDF4 file - ecad2nc4
## Go through the data conuntry and element wise to generate several netCDF files which 
## then can be combined into one. 

from.scratch <- TRUE
require(esd)
#source('~/R/esd/R/write2ncdf.R')
SS <- select.station(src='ecad')
variables <- ls()
it <- c(1800,as.numeric(format(Sys.Date(),'%Y')))

cntrs <- rownames(table(SS$country))
cntrs <- gsub(" ",".",cntrs)
cntrs <- gsub("[","",cntrs,fixed=TRUE)
cntrs <- gsub("]","",cntrs,fixed=TRUE)
cntrs <- gsub(",",".",cntrs,fixed=TRUE)
if (sum(is.element(variables,'eles'))==0)  eles <- rownames(table(SS$element))
if (sum(is.element(variables,'nmin'))==0) nmin <- 30

for (ele in eles) {
  ii <- 1 ## counter to keep track of number of stations saved
  param <- tolower(as.character(ele2param(ele,src='ecad')[5]))
  print(param)
  fname <- paste(param,'ecad','ncx',sep='.')
  if (file.exists(fname) & from.scratch) file.remove(fname)
  ## Read stations country by country and put the data in netCDF files according to that order
  for (cntr in cntrs) {
    print(cntr)
    ss <- select.station(src='ecad',cntr=cntr,param=param,nmin=nmin)   ## single country
    Ss <- select.station(src='ecad',param=param,nmin=nmin)             ## All countries
    append <- file.exists(fname)
    
    if (!is.null(ss)) {
      x <- station(cntr=cntr,param=param,src='ecad',save2file=FALSE)
      if (sum(!is.na(unit(x)))==0) {
        units <- switch(toupper(param),'SD'='cm','CC'='octas','RR'='mm/day','FX'='m/s',
                        'DD'='degree','FG'='m/s','PP'='hPa','SS'='hours','HU'='percent')
        attr(x,'unit') <- units
      }
      if (!is.null(dim(x))) {
        if (!append) stano <- 1:dim(Ss)[1] else stano <- ii:(ii+dim(x)[2]-1)
      } else if (!is.null(x)) {
        if (!append) stano <- 1:dim(Ss)[1] else stano <- ii
      }
      print(rbind(loc(x),firstyear(x),lastyear(x)))
      ## Quality check for the data in general
      if ( (min(x,na.rm=TRUE) < -999) | (max(x,na.rm=TRUE)>2000) ) {
        print("Detected suspect data")
        print(range(x,na.rm=TRUE))
        xc <- coredata(x); xc[xc < -999] <- NA; xc[xc > 2000] <- NA; coredata(x) <- as.matrix(xc)
        rm("xc"); gc(reset=TRUE)
      }
      
      if (dim(ss)[1]==1) {
        ## If there is only one station, then R does not treat the data as matrices ...
        save(x,file='temp.ecad2ncd4.rda')
      } else 
        if (dim(x)[2] > 0) {
          print(paste('Add',dim(x)[2],'stations to',fname))
          if (file.exists('temp.ecad2ncd4.rda')) {
            ## If a single station was previously aved, then add it to the next group and tidy up.  
            x2 <- x; load('temp.ecad2ncd4.rda')
            file.remove('temp.ecad2ncd4.rda')
            xx <- combine.stations(x,x2)
            x <- xx; rm("xx"); rm("x2")
          }
          ## additional quality check for temperature and precipitation:
          if (is.T(x)) {
            print('Quality check for temperature')
            xa <- coredata(anomaly(x))
            suspect <- abs(xa) > 40
            suspect[is.na(suspect)] <- FALSE
            print(paste('There were',sum(suspect),'suspect data points'))
            coredata(x)[suspect] <- NA
            rm('suspect'); gc(reset=TRUE)
          }
          if (is.precip(x)) {
            print('Quality check for precipitation')
            suspect <- (coredata(x) > 1000) | (coredata(x) < 0)
            suspect[is.na(suspect)] <- FALSE
            print(paste('There were',sum(suspect),'suspect data points'))
            coredata(x)[suspect] <- NA
            rm('suspect'); gc(reset=TRUE)
          }
          ## If there are many stations, there may be a memory problem - save the data in
          ## chuncks of 300 stations
          ns <- dim(x)[2]
          is  <- seq(1,ns,by=300)
          print('Saving stations in chuncks up to 300 stations at a time'); print(is)
          for (i in is) {
            xx <- subset(x,is=seq(i,min(ns,i+299),by=1))
            print(loc(xx))
            print(range(index(xx)))
            #map(xx,FUN='nv',new=FALSE)
            write2ncdf4(xx,file=fname,it=it,append=append,verbose=TRUE,stid_unlim=TRUE)
          }
        }
      if (!is.null(dim(x))) ii <- ii + dim(x)[2] else if (!is.null(x)) ii <- ii + 1
    } else x <- NULL
    rm("x"); gc(reset=TRUE)
    print(ii)
  }
  system(paste0('mv ',fname,' /lustre/storeB/project/ECAD/',substr(fname,1,nchar(fname)-1)))
}
