## A script that uses climate data operators (CDO) to produce monthly mean from
## daily mean data
## Rasmus Benestad
## (thanks to Oskar Landgren for advice)
## CDO source:
##  https://www.unidata.ucar.edu/software/netcdf/workshops/2012/third_party/CDO.html
dailyfield2monthlyfield <- function(fname,cline='cdo monmean',
                                    output='dailyfield2monthlyfield_out.nc') {
  print('This function only works if you have CDO installed')
  system(paste(cline,fname,output))
  x <- retireve(output)
  invisible(x)
}

## load a test-file:
url <- 'ftp://ftp.cdc.noaa.gov/Datasets.other/south_america/sa24.daily.1.1940-2012.nc'

## Only download if the data is not stored locally
if (!file.exists('sa24daily1940-2012.nc')) {
  print(paste('Download the predictand from',url))
  download.file(url,'sa24daily1940-2012.nc')
}

library(esd)
X <- dailyfield2monthlyfield('sa24daily1940-2012.nc')
map(X)
