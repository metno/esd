## A script that uses climate data operators (CDO) to produce monthly mean from
## daily mean data
## Rasmus Benestad
##  https://www.unidata.ucar.edu/software/netcdf/workshops/2012/third_party/CDO.html
dailyfield2monthlyfield <- function(fname,cline='cdo monmean',
                                    output='dailyfield2monthlyfield_out.nc') {
  print('This function only works if you have CDO installed')
  system(paste(cline,fname,output))
  x <- retireve(output)
  invisible(x)
}
