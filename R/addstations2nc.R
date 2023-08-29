#' Combine or append netCDF files for stations into one. Used with write2ncdf4.station to bypass memory demand/limitations.
#' The write2ncdf4.station with the argument append=TRUE also saves data chunk-wise. 
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
#' @param verbose used for debugging and diagnostics of executng the code
#'
#' @keywords utilities
#'
#' @examples 
#'
#' data(ferder)
#' plot(anomaly(ferder))
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
                                scale=0.1,missval=-99,torg='1899-12-31',it=NULL,stid=NULL,verbose=TRUE,
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
  if (verbose) {print('renamed statistics');print(statnames)}
  
  generate.station.ncfile(xval,out,stats,missval,offset,scale,torg,prec,it,verbose,stid_unlim,
                          namelength,calendar,doi,namingauthority,processinglevel,creatortype,
                          creatoremail,institution,publishername,publisheremail,
                          publisherurl,project,insufficient)

  
  if (verbose) print('nc4combine.stations successfuly finished')
}
