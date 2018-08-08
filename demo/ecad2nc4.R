## Save the ECA&D station data as a netCDF4 file - ecad2nc4
## Go through the data conuntry and element wise to generate several netCDF files which 
## then can be combined into one. 

#require(esd)
#source('~/R/esd/R/write2ncdf.R')
SS <- select.station(src='ecad')
variables <- ls()

cntrs <- rownames(table(SS$country))
cntrs <- gsub(" ",".",cntrs)
cntrs <- gsub("[","",cntrs,fixed=TRUE)
cntrs <- gsub("]","",cntrs,fixed=TRUE)
cntrs <- gsub(",",".",cntrs,fixed=TRUE)
if (sum(is.element(variables,'eles'))==0)  eles <- rev(rownames(table(SS$element)))
if (sum(is.element(variables,'nmin'))==0) nmin <- NULL

for (ele in eles) {
  ii <- 1 ## counter to keep track of number of stations saved
  param <- tolower(as.character(ele2param(ele,src='ecad')[5]))
  print(param)
  fname <- paste(param,'ecad','nc',sep='.')
  if (file.exists(fname)) file.remove(fname)
  for (cntr in cntrs) {
    print(cntr)
#    meta <- read.table(file.path(paste('data.ECAD/ECA_nonblend',param,sep='_'),'sources.txt'),
#                       skip=22,header=TRUE,sep=',')
    ss <- select.station(src='ecad',cntr=cntr,param=param,nmin=nmin)   ## single country
    Ss <- select.station(src='ecad',param=param,nmin=nmin)             ## All countries
    append <- file.exists(fname)
    
    if (!is.null(ss)) {
      x <- station(cntr=cntr,param=param,src='ecad',save2file=FALSE)
      if (!append) it <- seq(min(c(as.Date('1900-01-01'),index(x))),max(c(as.Date('2018-05-31'),index(x))),by='day')
      print(range(it)); print(range(index(subset(x,is=apply(x,1,'nv')>0))))
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
      print(rbind(loc(x),round(apply(x,2,'nv')/365.25),firstyear(x),lastyear(x)))
      ## Quality check
      if ( (min(x,na.rm=TRUE) < -999) | (max(x,na.rm=TRUE)>2000) ) {
        print("Detected suspect data")
        print(range(subset(x,is=apply(x,1,'nv')>0),na.rm=TRUE))
        xc <- coredata(x); xc[xc < -999] <- NA; xc[xc > 2000] <- NA; coredata(x) <- as.matrix(xc)
      }
      # if (length(x) > 0) write2ncdf4(x,fname,it=it,
      #                                stid=stano,append=append,verbose=FALSE,stid_unlim=TRUE)
      if (length(x) > 0) write2ncdf4(x,fname,it=it,append=append,verbose=FALSE,stid_unlim=TRUE)
      if (!is.null(dim(x))) ii <- ii + dim(x)[2] else if (!is.null(x)) ii <- ii + 1
    } else x <- NULL
    print(ii)
  }
}