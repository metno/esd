#' R-script that downloads daily data from the Copernicus Climate Data Store (CDS) using
#' the CDS set-up and python scripts through the API. The files will be stored
#' as netCDF files. This script assumes that CDO is installed: 
#' https://www.unidata.ucar.edu/software/netcdf/workshops/most-recent/third_party/CDO.html.
#' It only works on Linux platforms...
#' See https://cds.climate.copernicus.eu/api-how-to
#'
#' @aliases ERA5.CDS
#'
#' @param X variable name in CDS call, e.g. 'total_precipitation', '2m_temperature', or 'mean_sea_level_pressure'
#' @param it the years to extract.
#' @param varnm variable name for local data file.
#' @param AREA the area/region to extract [south,west,north,east]
#' @param FNAME the name of the local files for storing the data
#' @param FUN the function for CDO to aggregate the data, eg 'monsum', 'daymean',monmean', 'yearsum', 'yearmax', etc. If NULL, then leave the data as they are (e.g. daily data).  
#' @param verbose a boolean; if TRUE print information about progress
#' @examples
#' ERA5.CDS(X='2m_temperature',varnm='t2m',it=2015:2018,AREA="['50','0','60','10']",
#'          FUN='daymean')
#' ERA5.CDS(X='total_precipitation',varnm='tp',it=2018,AREA="['0','50','10','60']",
#'          FUN='yearsum')
#' ERA5.CDS(X='mean_sea_level_pressure',varnm='slp',it=2018,AREA="['40','-50','60','30']",
#'          FUN='monmean')
#'
#' @export
ERA5.CDS <- function(X='total_precipitation',it=1979:2018,
                     varnm='tp', AREA="['-90','-180','90','180']",
                     FNAME="'ERA5_XXX_YYYY.nc'",FUN='monsum',
                     path='~/Downloads/',verbose=TRUE) { 
  if (!file.exists('~/.cdsapirc')) {
    print('You need to install the CDS API key according to the web site and then re-run the call...')
    browser('https://cds.climate.copernicus.eu/api-how-to#install-the-cds-api-key')
    return()
  }
  dir <- getwd()
  setwd(path)
  FNAME <- sub('XXX',varnm,FNAME)
  if (verbose) print(FNAME)
  
  for (yr in it) {
    data("py.script")
    filename <- paste0('get-era5-',varnm,'_cds_',yr,'.py')
    py.script <- gsub('FNAME',FNAME,py.script)
    py.script <- gsub('YYYY',as.character(yr),py.script)
    py.script <- gsub('AREA',AREA,py.script)
    py.script <- gsub('XXX',X,py.script)
    writeLines(py.script,con=filename)
    #   print(py.script[13])
    rm('py.script')
    system(paste('python',filename))
    if (!is.null(FUN)) {
      ## If FUN is provided for aggregation:
      system(paste('cdo -b 64 ',FUN,gsub('YYYY',as.character(yr),FNAME),'aggregated.nc'))
      file.rename('aggregated.nc',gsub("'","",gsub('YYYY',as.character(yr),FNAME)))
    }
    print(gsub('YYYY',as.character(yr),FNAME))
    #file.remove(filename)
  }
  if (verbose) print('merge the years to single file')
  if (length(it)>1) { 
    system(paste('cdo -b 64 mergetime ',gsub('YYYY','????',FNAME),gsub('YYYY','',FNAME)))
    file.remove(gsub('YYYY','????',FNAME))
  }
  if (verbose) print(paste0('Download finished: ',path,'/',gsub('YYYY','',FNAME)))
  setwd(dir)
}
