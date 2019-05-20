#' Read cyclone data
#' 
#' Methods for reading cyclone data.
#' 
#' \code{read.imilast} reads data from files following the standards of the
#' IMILAST project (Neu et al., 2013: IMILAST: A Community Effort to
#' Intercompare Extratropical Cyclone Detection and Tracking Algorithms. Bull.
#' Amer. Meteor. Soc., 94, 529–547, https://doi.org/10.1175/BAMS-D-11-00154.1).
#' 
#' \code{read.hurdat2} reads data from the Atlantic hurricane database
#' (http://www.nhc.noaa.gov/data/#hurdat).
#' 
#' 
#' @aliases read.imilast read.hurdat2
#' @param fname filename (for \code{read.hurdat}, filename can also be a url)
#' @param path path to file
#' @param verbose Logical value defaulting to FALSE. If FALSE, do not display
#' comments (silent mode). If TRUE, displays extra information on progress.
#' @return An "events" "data.frame" object containing the date, time, lon, and
#' lat, as well as additional information (e.g., trajectory number, slp or
#' other measure of storm strengt) of the cyclones.
#' @author K. Parding
#' @keywords cyclones storms IMILAST hurdat2
#' @examples
#' 
#' 
#' @export read.imilast
read.imilast <- function(fname,path=NULL,verbose=FALSE) {
  if(!is.null(path)) fname <- file.path(path,fname)
  # read and rearrange file header
  if(verbose) print(paste('reading file header:'))
  h <- tolower(readLines(fname,1))
  if(verbose) print(h)
  if(verbose) print('rearranging')
  h <- unlist(strsplit(h,','))
  if(length(h)==1) h <- unlist(strsplit(h,' '))
  h <- gsub('^\\s+|\\s+$','',h)
  h[grepl('cyclone',h) | grepl('trackn',h)] <- 'trajectory'
  h[grepl('lati',h) | grepl('latn',h)] <- 'lat'
  h[grep('long',h)] <- 'lon'
  h[grepl('99',h) | grepl('code',h)] <- 'code99'
  h[grepl('date',h) | grepl('yyyymmddhh',h)] <- 'datetime'
  h[grep('yyyy',h)] <- 'year'
  h[grep('mm',h)] <- 'month'
  h[grep('dd',h)] <- 'day'
  h[grepl('hh',h) | grepl('timestep\\.',h) | grepl('timestep\\_',h)] <- 'time'
  h[grepl('step',h) | grepl('ptn',h)] <- 'timestep'
  h <- unique(h)
  if(verbose) print(paste(h,collapse=', '))
  # check width of columns
  l <- readLines(fname,4)[4]
  l <- gsub('^\\s+|\\s+$','',l)
  blanks <- unlist(gregexpr(' ',l))
  breaks <- blanks[!(blanks %in% as.integer(blanks+1))]
  w <- c(blanks[1]-1,
         breaks[2:length(breaks)]-breaks[1:(length(breaks)-1)],
         nchar(l)-breaks[length(breaks)]+1)
  if(verbose) print(paste('width of columns:',paste(w,collapse=',')))                
  # read data
  if(verbose) print('reading data')
  x <- read.fwf(fname,widths=w,col.names=h,skip=1)
  x <- x[x$code99<90,]
  x <- x[!is.na(x$lon),]
  # rearrange date and time information
  dates <- round(x['datetime'][[1]]*1E-2)
  times <- x['datetime'][[1]] - round(dates)*1E2
  x['date'] <- dates
  x['time'] <- times
  x <- x[,-which(names(x) %in% c('datetime','year','month','day','timestep'))]
  #add attributes
  param <- 'storm tracks'
  longname <- 'mid-latitude storm trajectories'
  src <- 'IMILAST'
  url <- 'http://journals.ametsoc.org/doi/abs/10.1175/BAMS-D-11-00154.1'
  ref <- paste('Neu et al. 2013: IMILAST: A Community Effor to Intercompare Extratroptical Cyclone Detection and Tracking Algorithms.',
               'Bull. Amer. Meteor. Soc., 94, 529-547.')
  if (x$code99[1]>9) {
    method <- paste('M',as.character(x$code99[1]),sep='')
  } else {
    method <- paste('M0',as.character(x$code99[1]),sep='')
  }
  x <- as.events(x,longname=longname,param=param,method=method,src=src,
                 reference=ref,file=file.path(path,fname),url=url,verbose=verbose)
  attr(x, 'history')= history.stamp() 
  invisible(x)
}

