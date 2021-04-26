## Estimate the wet-day mean rainfall for ERA5
## Rasmus Benestad, 2021-04-06

library(esd)

aggregatedailyfield <- function(fname,cline='cdo yearmean',
                                output='mu_year.nc') {
  print('This function only works if you have CDO installed')
  print(paste(cline,fname,output))
  system(paste(cline,fname,output))
  x <- retrieve(output)
  invisible(x)
}


wetdaymean <- function(ifile,cline='cdo yearmean',fname='wetdaymean.nc') {
  print('calculate mask')
  system(paste('cdo gtc,0.001',ifile,'mumask.nc')) ## ERA5 tp has units of m/day.
  ## Use mask on ifile so all pixles with non-zero values are set to 1 i all time steps
  ## and a value 0 i all others -> mask
  print('apply mask')
  system(paste('cdo ifthen mumask.nc',ifile, fname))
  ## Apply mask on ifile and store in ofile (set masked values to NA).
  X <- aggregatedailyfield(fname=fname,cline=cline)
  Y <- aggregatedailyfield(fname='mumask.nc',cline=cline,output='fw_year.nc')
  file.remove(fname); file.remove('mumask.nc')
}

fname <- readLine('Name (and full path) of netCDF file with daily precipitation >')
X <- wetdaymean(fname)
print(range(index(X)))
map(X)



