## Save the ECA&D station data as a netCDF4 file - ecad2nc4
## Go through the data conuntry and element wise to generate several netCDF files which 
## then can be combined into one. 

#require(esd)
#source('~/R/esd/R/write2ncdf.R')

variables <- ls()
it <- seq(as.Date('1900-01-01'),as.Date('2018-05-31'),by='day')

# cntrs <- rownames(table(SS$country))
# cntrs <- gsub(" ",".",cntrs)
# cntrs <- gsub("[","",cntrs,fixed=TRUE)
# cntrs <- gsub("]","",cntrs,fixed=TRUE)
# cntrs <- gsub(",",".",cntrs,fixed=TRUE)
if (sum(is.element(variables,'eles'))==0)  eles <- rev(rownames(table(SS$element)))
if (sum(is.element(variables,'nmin'))==0) nmin <- NULL

for (ele in eles) {
  SS <- select.station(src='ecad',ele=ele)
  stids <- SS$station_id
  ns <- length(stids)
  is <- seq(1,ns,by=50)
  if (max(is) < ns) is <- c(is,ns)
  
  param <- tolower(as.character(ele2param(ele,src='ecad')[5]))
  print(param)
  fname <- paste(param,'ecad','nc',sep='.')
  if (file.exists(fname)) file.remove(fname)
  for (id in is) {
    iii <- seq(id,id+49,length=50)
    iii <- iii[iii <= ns]
    
    append <- file.exists(fname)
    
    x <- try(station(SS[iii,],save2file=FALSE))
    
    #print(range(it)); print(range(index(subset(x,is=apply(x,1,'nv')>0))))
    if (!inherits(x,'try-error')) {
      units <- switch(toupper(param),'SD'='cm','CC'='octas','RR'='mm/day','FX'='m/s',
                      'DD'='degree','FG'='m/s','PP'='hPa','SS'='hours','HU'='percent')
      attr(x,'unit') <- units
      
      print(rbind(loc(x),firstyear(x),lastyear(x)))
      ## Quality check
      if ( (min(x,na.rm=TRUE) < -999) | (max(x,na.rm=TRUE)>2000) ) {
        print("Detected suspect data")
        print(range(x,na.rm=TRUE))
        xc <- coredata(x); xc[xc < -999] <- NA; xc[xc > 2000] <- NA; coredata(x) <- as.matrix(xc)
        rm("xc"); gc(reset=TRUE)
      }
      # if (length(x) > 0) write2ncdf4(x,fname,it=it,
      #                                stid=stano,append=append,verbose=FALSE,stid_unlim=TRUE)
      if (length(x) > 0) write2ncdf4(x,fname,it=it,append=append,verbose=FALSE,stid_unlim=TRUE)
      if (!is.null(dim(x))) ii <- ii + dim(x)[2] else if (!is.null(x)) ii <- ii + 1
    } else x <- NULL
    rm("x"); gc(reset=TRUE)
  }
}