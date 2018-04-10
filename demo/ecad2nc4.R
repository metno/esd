## Save the ECA&D station data as a netCDF4 file - ecad2nc4
## Go through the data conuntry and element wise to generate several netCDF files which 
## then can be combined into one. 

require(esd)
ss <- select.station(src='ecad')

cntrs <- rownames(table(ss$country))
eles <- rownames(table(ss$element))

for (ele in eles) {
  param <- tolower(as.character(ele2param(ele,src='ecad')[5]))
  print(param)
  for (cntr in cntrs) {
    ss <- select.station(src='ecad',cntr=cntr,param=param)
    if (!is.null(ss)) {
      fname <- paste(param,'ecad',cntr,'nc',sep='.')
      x <- station(cntr=cntr,param=param,src='ecad')
      if (length(x) > 0) write2ncdf4(x,fname,tim=seq(as.Date('1900-01-01'),as.Date('2018-02-28'),by=1))
    }
  }
}