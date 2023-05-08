#' R-script that downloads daily data from the Copernicus Climate Data Store (CDS) using
#' the CDS set-up and python scripts through the API. The files will be stored
#' as netCDF files. This script assumes that CDO and python are installed: 
#' https://www.unidata.ucar.edu/software/netcdf/workshops/most-recent/third_party/CDO.html.
#' It only works on Linux platforms...
#' See https://cds.climate.copernicus.eu/api-how-to
#'
#' @aliases ERA5.CDS
#'
#' @param param variable name in CDS call, e.g. 'total_precipitation', '2m_temperature', 'mean_sea_level_pressure',
#' '10m_u_component_of_wind', '10m_v_component_of_wind', 'relative_humidity', 'dewpoint_depression', 'snow_depth'
#' @param it the years to extract.
#' @param varnm variable name for local data file.
#' @param lon longitude of the area/region to extract.
#' @param lat latitude of the area/region to extract.
#' @param FNAME the name of the local files for storing the data
#' @param FUN the function for CDO to aggregate the data, eg 'monsum', 'daymean',monmean', 'yearsum',
#' 'yearmax', etc. If NULL, then leave the data as they are (e.g. daily data). If a vector (e.g. FUN=c('daymean','daymin'm',daymas')) it will 
#' use CDO repeated times to estimate each statistic.
#' @param cleanup If true, remove the original netCDF-file with hourly data to avoid clogging up the disc.
#' @param path The path where the data are stored. Can be a symbolic link.
#' @param python The version of python to use
#' @param verbose a boolean; if TRUE print information about progress
#' @examples
#' \dontrun{
#' ERA5.CDS(param='2m_temperature',varnm='t2m',it=2015:2018,lon=c(0,10), lat=c(50,60),
#'          FUN='daymean')
#' ERA5.CDS(param='total_precipitation',varnm='tp',it=2018,lon=c(50,60),lat=c(0,10),
#'          FUN='yearsum')
#' ERA5.CDS(param='mean_sea_level_pressure',varnm='slp',it=2018,lon=c(-50,30),lat=c(40,60),
#'          FUN='monmean')
#'}
#' @export
ERA5.CDS <- function(param='total_precipitation',it=1979:2018,
                     varnm=NULL, lon=c(-180,180),lat=c(-90,90),
                     FNAME="'ERA5_XXX_YYYY.nc'",FUN='monsum',cleanup=TRUE,
                     path='~/Downloads/',python='python3',verbose=TRUE) { 
  system('pip install cdsapi')
  AREA <- paste0("['",min(lat),"','",min(lon),"','",max(lat),"','",max(lon),"']")
  if (verbose) print(AREA)
  if (!file.exists('~/.cdsapirc')) {
    print('You need to install the CDS API key according to the web site and then re-run the call...')
    system('firefox https://cds.climate.copernicus.eu/api-how-to#install-the-cds-api-key')
    return()
  }
  dir <- getwd()
  setwd(path)
  if (is.null(varnm)) {
    if (sum(is.element(c("total_precipitation", "2m_temperature", "mean_sea_level_pressure",
                         "10m_u_component_of_wind", "10m_v_component_of_wind", "relative_humidity",
                         "dewpoint_depression", "snow_depth"),param)>0)) {
      varnm <- switch(param,"total_precipitation"='tp', "2m_temperature"='t2m', 
                      "mean_sea_level_pressure"='slp',
                      "10m_u_component_of_wind"='u10', "10m_v_component_of_wind"='v10',
                      "relative_humidity"='rh', "dewpoint_depression"='dpt', "snow_depth"='sd')
    } else varnm <-'x'
  }
  FNAME <- sub('XXX',varnm,FNAME)
  if (verbose) print(FNAME)
  
  for (yr in it) {
    data("py.script")
    filename <- paste0('get-era5-',varnm,'_cds_',yr,'.py')
    py.script <- gsub('FNAME',FNAME,py.script)
    py.script <- gsub('YYYY',as.character(yr),py.script)
    py.script <- gsub('AREA',AREA,py.script)
    py.script <- gsub('XXX',param,py.script)
    writeLines(py.script,con=filename)
    #     print(py.script[13])
    rm('py.script')
    if (verbose) print(paste0(python,' ./',filename))
    system(paste0(python,' ./',filename))
    if (!is.null(FUN)) {
      ## If FUN is provided for aggregation:
      for (ifs in 1:length(FUN)) {
        system(paste('cdo -b 64 ',FUN[ifs],gsub('YYYY',as.character(yr),FNAME),'aggregated.nc'))
        file.rename('aggregated.nc',sub('.nc',paste0('_',FUN[ifs],'.nc'),gsub("'","",gsub('YYYY',as.character(yr),FNAME)),fixed=TRUE))
      }
    }
    print(gsub('YYYY',as.character(yr),FNAME))
    if (cleanup) {
      file.remove(filename)
      file.remove(gsub('YYYY',as.character(yr),FNAME))
    }
  }
  if (verbose) print('merge the years to single file')
  if (length(it)>1) { 
    system(paste('cdo -b 64 mergetime ',gsub('YYYY','????',FNAME),gsub('YYYY','',FNAME)))
    file.remove(gsub('YYYY','????',FNAME))
    if (verbose) print(paste0('Download finished: ',path,'/',gsub('YYYY','',FNAME)))
  } else if (verbose) print(paste('Download finished:',path))
  
  setwd(dir)
}
