## Save the ECA&D station data as a netCDF4 file - ecad2nc4
## Go through the data conuntry and element wise to generate several netCDF files which 
## then can be combined into one. 

#require(esd)
#source('~/R/esd/R/write2ncdf.R')
SS <- select.station(src='ecad')

cntrs <- rownames(table(SS$country))
eles <- rownames(table(SS$element))
ii <- 1

for (ele in eles) {
  param <- tolower(as.character(ele2param(ele,src='ecad')[5]))
  print(param)
  for (cntr in cntrs) {
    ss <- select.station(src='ecad',cntr=cntr,param=param)   ## single country
    Ss <- select.station(src='ecad',param=param)             ## All countries
    fname <- paste(param,'ecad','nc',sep='.')
    append <- file.exists(fname)
    
    if (!is.null(ss)) {
      x <- station(cntr=cntr,param=param,src='ecad')
      if (!append) stano <- 1:dim(Ss)[1] else stano <- ii:(ii+dim(x)[2]-1)
      if (length(x) > 0) write2ncdf4(x,fname,tim=seq(as.Date('1900-01-01'),as.Date('2018-02-28'),by=1),
                                     stano=stano,append=append,verbose=TRUE)
      ii <- ii + dim(x)[2]
    }
  }
}