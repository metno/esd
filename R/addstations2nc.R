#' Combine or append netCDF files for stations into one. Used with write2ncdf4.station to bypass memory demand/limitations.
#' The \code{write2ncdf4.station} with the argument append=TRUE also saves data chunk-wise. Another function is 
#' \code{ncrcat} which reads all data and attributes in listed files and sace them into a new file for all the 
#' inputs combined. 
#' Rasmus.Benestad, met.no, 2023-07-31
#' 
#' @seealso write2ncdf4.station as.station retrieve.station 
#' 
#' @param x a station object
#' @param fname list with file names of netCDF files with station data
#' @param out name of the output file
#' @param namelength - see write2ncdf4.station
#' @param prec precision
#' @param offset offset value for the data
#' @param scale scaling value for the data
#' @param missval missing value
#' @param torg time origin
#' @param it index time - selection of time interval unless NULL
#' @param stid selection of station IDs unless set to NULL
#' @param verbose used for debugging and diagnostics of executing the code
#'
#' @keywords utilities
#'
#' @examples
#' \dontrun{
#'   ncfiles <- list.files(path='~/Downloads/ecad2nc4',pattern='.ecad.nc',full.names=TRUE)
#'   nc4combine.stations(ncfiles,out=paste0('~/Downloads/',variable,'.ecad.nc'))
#' }
#'
#' @export
nc4add.stations <- function(x,fname,verbose=FALSE) {
  ## This routine adds stations to a netCDF file (fname) by first saving a temporary nc-file and then
  ## combining it with the old file
  if (verbose) print('nc4add.stations')
  write2ncdf4.station(x,file='addstations2nc.temp.nc')
  combinencstations(fname=list(fname,'addstations2nc.temp.nc'),'addstations2nc.new.nc',verbose)
  file.remove('addstations2nc.temp.nc')
  file.rename('addstations2nc.new.nc',fname)
}

#' @export
nc4combine.stations <- function(fname,out='nc4combine.stations.nc',namelength=24,prec='short',offset=0, 
                                scale=0.1,missval=-99,torg='1899-12-31',it=NULL,stid=NULL,nv=NULL,verbose=TRUE,
                                stid_unlim=FALSE,doi='NA',namingauthority='NA',processinglevel='NA',
                                creatortype='NA',creatoremail='NA',institution='NA',publishername='NA',
                                publisheremail='NA',publisherurl='NA',project='NA') {
  if (verbose) print('nc4combine.stations')
  stopifnot(length(fname) > 1)
  if (is.list(fname)) filenames <- names(fname) else if (is.character(fname)) filenames <- fname
  
  ## Read the station data, metadata and summary statistics
  i <- 0
  stats <- list()
  for (filename in filenames) {
    i <- i + 1
    if (verbose) print(filename)
    if (i == 1) {
      xval <- retrieve.station(filename)
      meta <- retrieve.stationsummary(filename)
      summarystats <- setdiff(names(meta),names(attributes(xval)))
      for (stat in summarystats) stats[[stat]] <- meta[[stat]]
    } else {
      xval <- combine.stations(xval,retrieve.station(filename))
      meta <- retrieve.stationsummary(filename)
      for (stat in summarystats) stats[[stat]] <- c(stats[[stat]],meta[[stat]])
    }
    if (!is.null(nv)) {
      nvx <- apply(xval,2,FUN='nv')
      xval <- subset(xval,is=nvx > nv)
      if (verbose) print(paste('Keep only stations with more than',nv,'valid data points: ns=',dim(xval)[2]))
    }
  }
  if (verbose) print('Combined the stations with their summary statistics- now create a new single combined netCDF file')
  if (is.null(attr(xval,'calendar'))) {
    attr(xval,'calendar') <- 'standard'
  }
  
  if (!is.na(attr(xval,'calendar'))) {
    calendar <- attr(xval,'calendar')
  } else {
    calendar <- 'standard'
  }
  insufficient <- rep(FALSE,length(stid(xval)))
  
  ## rename some of the variables: 'mean' -> 'ave'; 'trend' -> 'td'; '_' -> '.'
  statnames <- names(stats)
  if (verbose) {print('original statistics');print(statnames)}
  for (rnm in statnames[grep('_',statnames)]) {
    if (verbose) print(rnm)
    stats[[gsub('_','.',rnm)]] <- stats[[rnm]]
    stats[[rnm]] <- NULL
  }
  statnames <- names(stats)
  for (rnm in statnames[grep('trend',statnames)]) {
    if (verbose) (print(rnm))
    if (length(grep('mean|freq',rnm))>0) stats[[sub('trend.','td',rnm)]] <- stats[[rnm]] else
      if (length(grep('sigma',rnm))>0) stats[[sub('trend.','t',rnm)]] <- stats[[rnm]] else
        stats[[sub('trend','td',rnm)]] <- stats[[rnm]]
      stats[[rnm]] <- NULL
  }
  statnames <- names(stats)
  for (rnm in setdiff(statnames[grep('mean',statnames)],statnames[grep('wet|dry',statnames)])) {
    if (verbose) print(rnm)
    stats[[sub('mean','ave',rnm)]] <- stats[[rnm]]
    stats[[rnm]] <- NULL
  }
  statnames <- names(stats)
  for (rnm in statnames[grep('DJF|MAM|JJA|SON',statnames)]) {
    if (verbose) print(rnm)
    stats[[tolower(rnm)]] <- stats[[rnm]]
    stats[[rnm]] <- NULL
  }
  statnames <- names(stats)
  for (rnm in statnames[grep('wetfreq',statnames)]) {
    if (verbose) print(rnm)
    stats[[sub('wetfreq','fw',rnm)]] <- stats[[rnm]]
    stats[[rnm]] <- NULL
  }
  statnames <- names(stats)
  for (rnm in statnames[grep('wetmean',statnames)]) {
    if (verbose) print(rnm)
    stats[[sub('wetmean','mu',rnm)]] <- stats[[rnm]]
    stats[[rnm]] <- NULL
  }
  statnames <- names(stats)
  for (rnm in statnames[grep('sd',statnames)]) {
    if (verbose) print(rnm)
    stats[[sub('sd','std',rnm)]] <- stats[[rnm]]
    stats[[rnm]] <- NULL
  }
  stats[['lr']] <- stats[['lastrains']]; stats[['lastrains']] <- NULL
  stats[['ld']] <- stats[['lastdry']]; stats[['lastdry']] <- NULL
  stats[['mx']] <- stats[['max']]; stats[['max']] <- NULL
  stats[['mn']] <- stats[['mn']]; stats[['min']] <- NULL
  stats[['nhr']] <- stats[['records']]; stats[['records']] <- NULL
  stats[['nlr']] <- stats[['lows']]; stats[['lows']] <- NULL
  stats[['mwsl']] <- stats[['mean.wetdur']]; stats[['mean.wetdur']] <- NULL
  stats[['mdsl']] <- stats[['mean.drydur']]; stats[['mean.drydur']] <- NULL
  
  statnames <- names(stats)
  if (verbose) {print('renamed statistics');print(statnames); print(table(esd::src(xval)))}
  
  generate.station.ncfile(xval,out,stats,missval,offset,scale,torg,prec,it,verbose,stid_unlim,
                          namelength,calendar,doi,namingauthority,processinglevel,creatortype,
                          creatoremail,institution,publishername,publisheremail,
                          publisherurl,project,insufficient)
  
  if (verbose) print('nc4combine.stations successfuly finished')
}

## This function copies information in netCDF files into one assuming the same time
## coordinates - used for combining station data
#' @export
ncrcat.station <- function(fnames,ncout,namelength=24,verbose=TRUE) {
  ## Start with the first file in the list of files and use this as a template for the final file
  ## We assume that the netCDF files have the same structure (nc-files for stations) keeping the same 
  ## attributes.
  if (verbose) {print('ncrcat'); print(fnames)}
  
  if (verbose) print('Analyse the time spans in the netCDF files')
  ## Collect the times from all the files
  t.span <- c(); ns <- 0; lats <- c(); lons <- c(); alts <- c()
  for (fname in fnames) {
    #if (verbose) print(fname)
    ncadd <- nc_open(fname)
    ## Read the time information
    t.1 <- ncvar_get(ncid,'time')
    tunit.1 <- ncatt_get(ncid,'time','units')$value
    #if (verbose) print(sub('days since ','',tunit.1))
    if (length(t.span)==0) t.span <- range(as.Date(t.1,orig=sub('days since ','',tunit.1))) else
      t.span <- range(c(t.span,as.Date(t.1,orig=sub('days since ','',tunit.1))))
    stid <- ncvar_get(ncadd,'stationID')
    lats <- c(lats,ncvar_get(ncadd,'lat'))
    lons <- c(lons,ncvar_get(ncadd,'lon'))
    alts <- c(lons,ncvar_get(ncadd,'alt'))
    ns <- ns + length(stid)
    class.x <- ncatt_get( ncid, 0, 'class')$value
    nc_close(ncadd)
  }
  if (verbose) {print(t.span); print(ns)}
  
  time <- seq(t.span[1],t.span[2],by=1)
  nt <- length(t.new)
  
  ## Create a new netCDF file that spans all stations and times in the netCDF files
  
  dimS <- ncdim_def( name="stid", units="number",vals=1:ns,unlim=stid_unlim)
  dimT <- ncdim_def( name="time", units=paste("days since",torg), vals=time, calendar=calendar,unlim=TRUE)
  dimnchar   <- ncdim_def("nchar",   "", 1:namelength, create_dimvar=FALSE )
  latid <- ncvar_def(name="lat",dim=list(dimS), units="degrees_north", missval=missval,longname="latitude", 
                     prec="float",verbose=verbose)
  
  lonid <- ncvar_def(name="lon",dim=list(dimS), units="degrees_east", missval=missval,longname="longitude", 
                     prec="float",verbose=verbose)
  altid <- ncvar_def(name="alt",dim=list(dimS), units="meters", missval=missval,longname="altitude", 
                     prec=prec,verbose=verbose)
  
  locid <- ncvar_def(name="loc",dim=list(dimnchar,dimS),units="name",prec="char",longname="location",
                     verbose=verbose)
  stid <- ncvar_def(name="stationID",dim=list(dimnchar,dimS),units="number",prec="char",longname="station_id",
                    verbose=verbose)
  cntrid <- ncvar_def(name="cntr",dim=list(dimnchar,dimS),units="name",prec="char",longname="country",
                      verbose=verbose)
  ncvar <- ncvar_def(name=varid(x)[1],dim=list(dimS,dimT), units=ifelse(unitx[1]=="\u00B0C", "degC",unitx[1]),
                     longname=longname, prec=prec,compression=9,verbose=verbose)
  fyrid <- ncvar_def(name="first",dim=list(dimS), units="year", missval=missval,longname="first_year", 
                     prec="short",verbose=verbose)
  lyrid <- ncvar_def(name="last",dim=list(dimS), units="year", missval=missval,longname="last_year", 
                     prec="short",verbose=verbose)
  doid <- ncvar_def(name="days_old",dim=list(dimS), units="day", missval=missval,longname="days since last record and saving date", 
                    prec="integer",verbose=verbose)
  nvid <- ncvar_def(name="number",dim=list(dimS), units="count", missval=missval,longname="number_valid_data", 
                    prec="float",verbose=verbose)
  
  ## KMP 2018-11-02: devtools (run_examples) can only handle ASCII characters so I had to replace the 
  ## degree symbol with "\u00B0", but I'm not sure if it is going to work here.
  maxid <- ncvar_def(name="summary_max",dim=list(dimS), units= ifelse(attr(x,"unit")[1]=="\u00B0C", "degC",attr(x,"unit")[1]),
                     missval=missval,longname=varid(x)[1],prec="float",verbose=verbose)
  minid <- ncvar_def(name="summary_min",dim=list(dimS), ifelse(attr(x,"unit")[1]=="\u00B0C", "degC",attr(x,"unit")[1]), 
                     missval=missval,longname=varid(x)[1],prec="float",verbose=verbose)
  nhrid <- ncvar_def(name="summary_records",dim=list(dimS),
                     units=ifelse(attr(x,"unit")[1]=="\u00B0C", "degC",attr(x,"unit")[1]), 
                     missval=missval,longname="fraction_of_high_records",prec="float",verbose=verbose)
  lehrid <- ncvar_def(name="last_element_highest",dim=list(dimS),
                      units=ifelse(attr(x,"unit")[1]=="\u00B0C", "degC",attr(x,"unit")[1]), 
                      missval=missval,longname="If_last_element_is_a_record",prec="short",verbose=verbose)
  if (param=='t2m') {
    meanid <- ncvar_def(name="summary_mean",dim=list(dimS), units="degC", 
                        missval=missval,longname="annual_mean_temperature",prec="float",verbose=verbose)
    meanid.djf <- ncvar_def(name="summary_mean_DJF",dim=list(dimS), units="degC", 
                            missval=missval,longname="seasonal_mean_temperature_Dec-Feb",prec="float",verbose=verbose)
    meanid.mam <- ncvar_def(name="summary_mean_MAM",dim=list(dimS), units="degC", 
                            missval=missval,longname="seasonal_mean_temperature_Mar-May",prec="float",verbose=verbose)
    meanid.jja <- ncvar_def(name="summary_mean_JJA",dim=list(dimS), units="degC", 
                            missval=missval,longname="seasonal_mean_temperature_Jun-Aug",prec="float",verbose=verbose)
    meanid.son <- ncvar_def(name="summary_mean_SON",dim=list(dimS), units="degC", 
                            missval=missval,longname="seasonal_mean_temperature_Sep-Nov",prec="float",verbose=verbose)
    sdid <- ncvar_def(name="summary_sd",dim=list(dimS), units="degC", 
                      missval=missval,longname="annual_temperature_anomaly_standard_deviation",prec="float",verbose=verbose)
    sdid.djf <- ncvar_def(name="summary_sd_DJF",dim=list(dimS), units="degC", 
                          missval=missval,longname="temperature_anomaly_standard_deviation_Dec-Feb",prec="float",verbose=verbose)
    sdid.mam <- ncvar_def(name="summary_sd_MAM",dim=list(dimS), units="degC", 
                          missval=missval,longname="temperature_anomaly_standard_deviation_Mar-May",prec="float",verbose=verbose)
    sdid.jja <- ncvar_def(name="summary_sd_JJA",dim=list(dimS), units="degC", 
                          missval=missval,longname="temperature_anomaly_standard_deviation_Jun-Aug",prec="float",verbose=verbose)
    sdid.son <- ncvar_def(name="summary_sd_SON",dim=list(dimS), units="degC", 
                          missval=missval,longname="temperature_anomaly_standard_deviation_Sep-Nov",prec="float",verbose=verbose)
    nlrid <- ncvar_def(name="summary_lows",dim=list(dimS), units="degC", 
                       missval=missval,longname="fraction_of_low_records",prec="float",verbose=verbose)
    tdid <- ncvar_def(name="summary_trend",dim=list(dimS), units="degC/decade", 
                      missval=missval,longname="annual_mean_temperature",prec="float",verbose=verbose)
    tdid.djf <- ncvar_def(name="summary_trend_DJF",dim=list(dimS), units="degC/decade", 
                          missval=missval,longname="seasonal_mean_temperature_Dec-Feb",prec="float",verbose=verbose)
    tdid.mam <- ncvar_def(name="summary_trend_MAM",dim=list(dimS), units="degC/decade", 
                          missval=missval,longname="seasonal_mean_temperature_Mar-May",prec="float",verbose=verbose)
    tdid.jja <- ncvar_def(name="summary_trend_JJA",dim=list(dimS), units="degC/decade", 
                          missval=missval,longname="seasonal_mean_temperature_Jun-Aug",prec="float",verbose=verbose)
    tdid.son <- ncvar_def(name="summary_trend_SON",dim=list(dimS), units="degC/decade", 
                          missval=missval,longname="seasonal_mean_temperature_Sep-Nov",prec="float",verbose=verbose)
    ## KMP 2018-11-02: devtools (run_examples) can only handle ASCII characters so I had to replace the 
    ## degree symbol with "\u00B0", but I'm not sure if it is going to work here.
    lelrid <- ncvar_def(name="last_element_lowest",dim=list(dimS), 
                        units=ifelse(attr(x,"unit")[1]=="\u00B0C", "degC",attr(x,"unit")[1]), 
                        missval=missval,longname="If_last_element_is_a_record",prec="short",verbose=verbose)
    
  } else if (param=='precip') {
    meanid <- ncvar_def(name="summary_mean",dim=list(dimS), units="mm/year", 
                        missval=missval,longname="mean_annual_precipitation_sum",prec="float",verbose=verbose)
    meanid.djf <- ncvar_def(name="summary_mean_DJF",dim=list(dimS), units="mm/season", 
                            missval=missval,longname="mean_seasonal_precip_sum_Dec-Feb",prec="float",verbose=verbose)
    meanid.mam <- ncvar_def(name="summary_mean_MAM",dim=list(dimS), units="mm/season", 
                            missval=missval,longname="mean_seasonal_precip_sum_Mar-May",prec="float",verbose=verbose)
    meanid.jja <- ncvar_def(name="summary_mean_JJA",dim=list(dimS), units="mm/season", 
                            missval=missval,longname="mean_seasonal_precip_sum_Jun-Aug",prec="float",verbose=verbose)
    meanid.son <- ncvar_def(name="summary_mean_SON",dim=list(dimS), units="mm/season", 
                            missval=missval,longname="mean_seasonal_precip_sum_Sep-Nov",prec="float",verbose=verbose)
    muid <- ncvar_def(name="summary_wetmean",dim=list(dimS), units="mm/day", 
                      missval=missval,longname="mean_annual_precipitation_wet_day_mean",prec="float",verbose=verbose)
    muid.djf <- ncvar_def(name="summary_wetmean_DJF",dim=list(dimS), units="mm/day", 
                          missval=missval,longname="mean_seasonal_precip_wetmean_Dec-Feb",prec="float",verbose=verbose)
    muid.mam <- ncvar_def(name="summary_wetmean_MAM",dim=list(dimS), units="mm/day", 
                          missval=missval,longname="mean_seasonal_precip_wetmean_Mar-May",prec="float",verbose=verbose)
    muid.jja <- ncvar_def(name="summary_wetmean_JJA",dim=list(dimS), units="mm/day", 
                          missval=missval,longname="mean_seasonal_precip_wetmean_Jun-Aug",prec="float",verbose=verbose)
    muid.son <- ncvar_def(name="summary_wetmean_SON",dim=list(dimS), units="mm/day", 
                          missval=missval,longname="mean_seasonal_precip_wetmean_Sep-Nov",prec="float",verbose=verbose)
    fwid <- ncvar_def(name="summary_wetfreq",dim=list(dimS), units="mm/day", 
                      missval=missval,longname="mean_annual_precipitation_wet_day_frequency",prec="float",verbose=verbose)
    fwid.djf <- ncvar_def(name="summary_wetfreq_DJF",dim=list(dimS), units="%", 
                          missval=missval,longname="mean_seasonal_precip_wetfreq_Dec-Feb",prec="float",verbose=verbose)
    fwid.mam <- ncvar_def(name="summary_wetfreq_MAM",dim=list(dimS), units="%", 
                          missval=missval,longname="mean_seasonal_precip_wetfreq_Mar-May",prec="float",verbose=verbose)
    fwid.jja <- ncvar_def(name="summary_wetfreq_JJA",dim=list(dimS), units="%", 
                          missval=missval,longname="mean_seasonal_precip_wetfreq_Jun-Aug",prec="float",verbose=verbose)
    fwid.son <- ncvar_def(name="summary_wetfreq_SON",dim=list(dimS), units="%", 
                          missval=missval,longname="mean_seasonal_precip_wetfreq_Sep-Nov",prec="float",verbose=verbose)
    tdid <- ncvar_def(name="summary_trend",dim=list(dimS), units="annual mm/decade", 
                      missval=missval,longname="trend_in_annual_precipitation_sum",prec="float",verbose=verbose)
    tdid.djf <- ncvar_def(name="summary_trend_DJF",dim=list(dimS), units="seasonal mm/decade", 
                          missval=missval,longname="trend_in_seasonal_precip_sum_Dec-Feb",prec="float",verbose=verbose)
    tdid.mam <- ncvar_def(name="summary_trend_MAM",dim=list(dimS), units="seasonal mm/decade", 
                          missval=missval,longname="trend_in_seasonal_precip_sum_Mar-May",prec="float",verbose=verbose)
    tdid.jja <- ncvar_def(name="summary_trend_JJA",dim=list(dimS), units="seasonal mm/decade", 
                          missval=missval,longname="trend_in_seasonal_precip_sum_Jun-Aug",prec="float",verbose=verbose)
    tdid.son <- ncvar_def(name="summary_trend_SON",dim=list(dimS), units="seasonal mm/decade", 
                          missval=missval,longname="trend_in_seasonal_precip_sum_Sep-Nov",prec="float",verbose=verbose)
    tdmuid <- ncvar_def(name="summary_trend_wetmean",dim=list(dimS), units="mm/day/decade", 
                        missval=missval,longname="annual_precipitation_wetday_mean",prec="float",verbose=verbose)
    tdmuid.djf <- ncvar_def(name="summary_trend_wetmean_DJF",dim=list(dimS), units="mm/day/decade", 
                            missval=missval,longname="trend_in_seasonal_precip_wetmean_Dec-Feb",prec="float",verbose=verbose)
    tdmuid.mam <- ncvar_def(name="summary_trend_wetmean_MAM",dim=list(dimS), units="mm/day/decade", 
                            missval=missval,longname="trend_in_seasonal_precip_wetmean_Mar-May",prec="float",verbose=verbose)
    tdmuid.jja <- ncvar_def(name="summary_trend_wetmean_JJA",dim=list(dimS), units="mm/day/decade", 
                            missval=missval,longname="trend_in_seasonal_precip_wetmean_Jun-Aug",prec="float",verbose=verbose)
    tdmuid.son <- ncvar_def(name="summary_trend_wetmean_SON",dim=list(dimS), units="mm/day/decade", 
                            missval=missval,longname="trend_in_seasonal_precip_wetmean_Sep-Nov",prec="float",verbose=verbose)
    tdfwid <- ncvar_def(name="summary_trend_wetfreq",dim=list(dimS), units="%/decade", 
                        missval=missval,longname="annual_precipitation_wet_day_frequency",prec="float",verbose=verbose)
    tdfwid.djf <- ncvar_def(name="summary_trend_wetfreq_DJF",dim=list(dimS), units="%/decade", 
                            missval=missval,longname="trend_in_seasonal_precip_wet_frequency_Dec-Feb",prec="float",verbose=verbose)
    tdfwid.mam <- ncvar_def(name="summary_trend_wetfreq_MAM",dim=list(dimS), units="%/decade", 
                            missval=missval,longname="trend_in_seasonal_precip_wet_frequency_Mar-May",prec="float",verbose=verbose)
    tdfwid.jja <- ncvar_def(name="summary_trend_wetfreq_JJA",dim=list(dimS), units="%/decade", 
                            missval=missval,longname="trend_in_seasonal_precip_wet_frequency_Jun-Aug",prec="float",verbose=verbose)
    tdfwid.son <- ncvar_def(name="summary_trend_wetfreq_SON",dim=list(dimS), units="%/decade", 
                            missval=missval,longname="trend_in_seasonal_precip_wet_frequency_Sep-Nov",prec="float",verbose=verbose)
    lrid <- ncvar_def(name="summary_lastrains",dim=list(dimS), units="days", 
                      missval=missval,longname="number_of_dry_days_at_end_of_record",prec="float",verbose=verbose)
    ldid <- ncvar_def(name="summary_lastdry",dim=list(dimS), units="days", 
                      missval=missval,longname="number_of_wet_days_at_end_of_record",prec="float",verbose=verbose)
    sigma2id <- ncvar_def(name="summary_sigma2",dim=list(dimS), units="mm^2", 
                          missval=missval,longname="variance_daily_precip",prec="float",verbose=verbose)
    sigma2id.djf <- ncvar_def(name="summary_sigma2_DJF",dim=list(dimS), units="mm^2", 
                              missval=missval,longname="variance_daily_precip_Dec-Feb",prec="float",verbose=verbose)
    sigma2id.mam <- ncvar_def(name="summary_sigma2_MAM",dim=list(dimS), units="mm^2", 
                              missval=missval,longname="variance_daily_precip_Mar-May",prec="float",verbose=verbose)
    sigma2id.jja <- ncvar_def(name="summary_sigma2_JJA",dim=list(dimS), units="mm^2", 
                              missval=missval,longname="variance_daily_precip_Jun-Aug",prec="float",verbose=verbose)
    sigma2id.son <- ncvar_def(name="summary_sigma2_SON",dim=list(dimS), units="mm^2", 
                              missval=missval,longname="variance_daily_precip_Sep-Nov",prec="float",verbose=verbose)
    tsigma2id <- ncvar_def(name="summary_trend_sigma2",dim=list(dimS), units="mm^2", 
                           missval=missval,longname="variance_daily_precip_trend",prec="float",verbose=verbose)
    tsigma2id.djf <- ncvar_def(name="summary_trend_sigma2_DJF",dim=list(dimS), units="mm^2", 
                               missval=missval,longname="variance_daily_precip_trend_Dec-Feb",prec="float",verbose=verbose)
    tsigma2id.mam <- ncvar_def(name="summary_trend_sigma2_MAM",dim=list(dimS), units="mm^2", 
                               missval=missval,longname="variance_daily_precip_trend_Mar-May",prec="float",verbose=verbose)
    tsigma2id.jja <- ncvar_def(name="summary_trend_sigma2_JJA",dim=list(dimS), units="mm^2", 
                               missval=missval,longname="variance_daily_precip_trend_Jun-Aug",prec="float",verbose=verbose)
    tsigma2id.son <- ncvar_def(name="summary_trend_sigma2_SON",dim=list(dimS), units="mm^2", 
                               missval=missval,longname="variance_daily_precip_trend_Sep-Nov",prec="float",verbose=verbose)
    mwslid <- ncvar_def(name="summary_mean_wetdur",dim=list(dimS), units="day", 
                        missval=missval,longname="mean_wet-day-spell_length",prec="float",verbose=verbose)
    mdslid <- ncvar_def(name="summary_mean_drydur",dim=list(dimS), units="day", 
                        missval=missval,longname="mean_dry-spell_length",prec="float",verbose=verbose)
  } else {
    meanid <- ncvar_def(name="summary_mean",dim=list(dimS), units=attr(x,"unit")[1], 
                        missval=missval,longname=paste("mean_annual",varid(x)[1],sep='_'),prec="float",verbose=verbose)
    meanid.djf <- ncvar_def(name="summary_mean_DJF",dim=list(dimS), units=attr(x,"unit")[1], 
                            missval=missval,longname=paste("Dec-Feb_mean",varid(x)[1],sep='_'),prec="float",verbose=verbose)
    meanid.mam <- ncvar_def(name="summary_mean_MAM",dim=list(dimS), units=attr(x,"unit")[1], 
                            missval=missval,longname=paste("Mar-May_mean",varid(x)[1],sep='_'),prec="float",verbose=verbose)
    meanid.jja <- ncvar_def(name="summary_mean_JJA",dim=list(dimS), units=attr(x,"unit")[1], 
                            missval=missval,longname=paste("Jun-Aug_mean",varid(x)[1],sep='_'),prec="float",verbose=verbose)
    meanid.son <- ncvar_def(name="summary_mean_SON",dim=list(dimS), units=attr(x,"unit")[1], 
                            missval=missval,longname=paste("Sep-Nov_mean",varid(x)[1],sep='_'),prec="float",verbose=verbose)
    sdid <- ncvar_def(name="summary_sd",dim=list(dimS), units="degC", 
                      missval=missval,longname="annual_temperature_anomaly_standard_deviation",prec="float",verbose=verbose)
    sdid.djf <- ncvar_def(name="summary_sd_DJF",dim=list(dimS), units="degC", 
                          missval=missval,longname="temperature_anomaly_standard_deviation_Dec-Feb",prec="float",verbose=verbose)
    sdid.mam <- ncvar_def(name="summary_sd_MAM",dim=list(dimS), units="degC", 
                          missval=missval,longname="temperature_anomaly_standard_deviation_Mar-May",prec="float",verbose=verbose)
    sdid.jja <- ncvar_def(name="summary_sd_JJA",dim=list(dimS), units="degC", 
                          missval=missval,longname="temperature_anomaly_standard_deviation_Jun-Aug",prec="float",verbose=verbose)
    sdid.son <- ncvar_def(name="summary_sd_SON",dim=list(dimS), units="degC", 
                          missval=missval,longname="temperature_anomaly_standard_deviation_Sep-Nov",prec="float",verbose=verbose)
    tdid <- ncvar_def(name="summary_trend",dim=list(dimS), units=attr(x,"unit")[1], 
                      missval=missval,longname="annual_mean_temperature",prec="float",verbose=verbose)
    tdid.djf <- ncvar_def(name="summary_trend_DJF",dim=list(dimS), units=attr(x,"unit")[1], 
                          missval=missval,longname=paste("Dec-Feb_mean",varid(x)[1],sep='_'),prec="float",verbose=verbose)
    tdid.mam <- ncvar_def(name="summary_trend_MAM",dim=list(dimS), units=attr(x,"unit")[1], 
                          missval=missval,longname=paste("Mar-May_mean",varid(x)[1],sep='_'),prec="float",verbose=verbose)
    tdid.jja <- ncvar_def(name="summary_trend_JJA",dim=list(dimS), units=attr(x,"unit")[1], 
                          missval=missval,longname=paste("Jun-Aug_mean",varid(x)[1],sep='_'),prec="float",verbose=verbose)
    tdid.son <- ncvar_def(name="summary_trend_SON",dim=list(dimS), units=attr(x,"unit")[1], 
                          missval=missval,longname=paste("Sep-Nov_mean",varid(x)[1],sep='_'),prec="float",verbose=verbose)
    ## KMP 2018-11-02: devtools (run_examples) can only handle ASCII characters so I had to replace the 
    ## degree symbol with "\u00B0", but I'm not sure if it is going to work here.
    lelrid <- ncvar_def(name="last_element_lowest",dim=list(dimS), 
                        units=ifelse(attr(x,"unit")[1]=="\u00B0C", "degC",attr(x,"unit")[1]), 
                        missval=missval,longname="If_last_element_is_a_record",prec="short",verbose=verbose)
    nlrid <- ncvar_def(name="summary_lows",dim=list(dimS), units="degC", 
                       missval=missval,longname="fraction_of_low_records",prec="float",verbose=verbose)
  }
  
  if (param=='t2m') {
    ncid <- nc_create(file,vars=list(ncvar,lonid,latid,altid,locid,stid,cntrid, 
                                     fyrid,lyrid,nvid,meanid,meanid.djf,meanid.mam,meanid.jja,meanid.son,
                                     sdid,sdid.djf,sdid.mam,sdid.jja,sdid.son,maxid,minid,nhrid,nlrid,
                                     tdid,tdid.djf,tdid.mam,tdid.jja,tdid.son,lehrid,lelrid,doid))
  } else if (param=='precip') {
    ncid <- nc_create(file,vars=list(ncvar,lonid,latid,altid,locid,stid,cntrid, 
                                     fyrid,lyrid,nvid,meanid,meanid.djf,meanid.mam,meanid.jja,meanid.son,
                                     maxid,minid,nhrid,muid,muid.djf,muid.mam,muid.jja,muid.son,
                                     fwid,fwid.djf,fwid.mam,fwid.jja,fwid.son,
                                     tdid,tdid.djf,tdid.mam,tdid.jja,tdid.son,
                                     tdmuid,tdmuid.djf,tdmuid.mam,tdmuid.jja,tdmuid.son,
                                     tdfwid,tdfwid.djf,tdfwid.mam,tdfwid.jja,tdfwid.son,lrid,ldid,lehrid,
                                     sigma2id,sigma2id.djf,sigma2id.mam,sigma2id.jja,sigma2id.son,
                                     tsigma2id,tsigma2id.djf,tsigma2id.mam,tsigma2id.jja,tsigma2id.son,
                                     mwslid,mdslid,doid))
  } else {
    ncid <- nc_create(file,vars=list(ncvar,lonid,latid,altid,locid,stid,cntrid, 
                                     fyrid,lyrid,nvid,meanid,meanid.djf,meanid.mam,meanid.jja,meanid.son,
                                     sdid,sdid.djf,sdid.mam,sdid.jja,sdid.son,tdid,
                                     tdid.djf,tdid.mam,tdid.jja,tdid.son,maxid,minid,nhrid,nlrid,lehrid,
                                     lelrid,doid))
  }
  
  ## Set the global attributes
  ncatt_put( ncid, 0, 'class', class.x)
  ncatt_put( ncid, 0, 'title', paste(levels(factor(attr(x,"info"))),collapse="/"))
  ncatt_put( ncid, 0, 'source', src)
  ncatt_put( ncid, 0, 'history', paste(unlist(attr(x,"history")),collapse="/"))
  ncatt_put( ncid, 0, 'references', paste(levels(factor(attr(x,"reference"))),collapse="/"))
  ncatt_put( ncid, 0, "esd-version", attr(x,'history')$session$esd.version)
  ## Add global attributes suggested in https://adc.met.no/node/4
  ncatt_put( ncid, 0, 'id', doi)
  ncatt_put( ncid, 0, 'naming_authority', namingauthority)
  ncatt_put( ncid, 0, 'title', 'netCDF-file for station data')
  ncatt_put( ncid, 0, 'summary', 'Daily data organised as station ID number and time')
  ncatt_put( ncid, 0, 'keywords', 'Observational record, summary statistics, station data')
  ncatt_put( ncid, 0, 'geospatial_lat_min', min(lats) )
  ncatt_put( ncid, 0, 'geospatial_lat_max', max(lats) )
  ncatt_put( ncid, 0, 'geospatial_lon_min', min(lons) )
  ncatt_put( ncid, 0, 'geospatial_lon_max', max(lons) )
  ncatt_put( ncid, 0, 'time_coverage_start', as.character(min(time)) )
  ncatt_put( ncid, 0, 'time_coverage_end', as.character(max(time)) )
  ncatt_put( ncid, 0, 'Conventions', 'CF-inspired')
  ncatt_put( ncid, 0, 'processing_level', processinglevel)
  ncatt_put( ncid, 0, 'date_created', date())
  ncatt_put( ncid, 0, 'creator_type', creatortype)
  ncatt_put( ncid, 0, 'creator_email', creatoremail)
  ncatt_put( ncid, 0, 'creator_url', 'github.com/metno/esd')
  ncatt_put( ncid, 0, 'institution', institution)
  ncatt_put( ncid, 0, 'publisher_name', publishername)
  ncatt_put( ncid, 0, 'publisher_email', publisheremail)
  ncatt_put( ncid, 0, 'publisher_url', publisherurl)
  ncatt_put( ncid, 0, 'project', project)
  
  ## Now read the data, insert it into the correct time period and save
  start.is <- 1
  for (fname in fnames) {
    if (verbose) print(fname)
    ncadd <- nc_open(fname)
    ## Read the time information
    t.add <- ncvar_get(ncadd,'time')
    #if (max(diff(t.add))>1) print(paste('ncrcat - Warning time steps greater then one day',fname))
    tunit.add <- ncatt_get(ncadd,'time','units')$value
    t.add <- as.Date(t.add,orig=sub('days since ','',tunit.add))
    if (verbose) print(names(ncadd$var))
    ## Read the variables
    for (var in names(ncadd$var)) {
      if (verbose) print(var)
      x <- ncvar_get(ncadd,var)
      d <- dim(x)
      ## dimensions read is space - time
      stid <- ncvar_get(ncadd,'stationID')
      ns <- length(stid)
      ## For 2-D data matrices
      if (length(d)==2) { 
        scale <- try(ncatt_get(ncadd,var,'scale_factor')$value)
        offset  <- try(ncatt_get(ncadd,var,'add_offset')$value)
        missing <- try(ncatt_get(ncadd,var,'missing_value')$value)
        if (!inherits(scale,'try-error') & !inherits(offset,'try-error'))
          x <- (x - offset)/scale
        ## Locations
        y <- matrix(rep(NA,ns*nt),ns,nt)
        i1 <- is.element(t.add,t.new)
        i2 <- is.element(t.new,t.add)
        y[,i2] <- x[,i1]
        x <- y; rm('y')
        start <- c(start.is,1)
        count <- c(ns,nt)
        if (verbose) print(paste('Add',var,'from',fname,' start=c(',start[1],start[2],'), count=c(',
                                 count[1],count[2],')'))
      } else {
        start <- start.is
        count <- ns
        if (verbose) print(paste('Add',var,'from',fname,' start=',start, 'count=',count))
      }
      is <- ncvar_get(ncadd,'stid')
      if (!is.character(x)) ncvar_put(ncid,var,x,start=start,count=count) else
        ncvar_put(ncid,var,x,start=c(1,start),count=c(namelength,count))
      ## Read all the attributes
      if (start.is==1) {
        browser()
        vatts <- names(ncadd$var) 
      }
    }
    nc_close(ncadd)
    start.is <- start.is + ns
  }
  nc_close(ncid)
  if (verbose) print(paste('Finished - combined data is in',ncout))
}
