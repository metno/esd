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
#' @param verbose a boolean; If FALSE, do not display
#' comments (silent mode). If TRUE, display extra information on progress.
#' @return An "events" "data.frame" object containing the date, time, lon, and
#' lat, as well as additional information (e.g., trajectory number, slp or
#' other measure of storm strength) of the cyclones.
#' @author K. Parding
#' @keywords cyclones storms IMILAST hurdat2
#' @examples
#' 
#' 
#' @export
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

#fname <- 'http://www.aoml.noaa.gov/hrd/hurdat/Data_Storm.html'
read.hurdat2 <- function(fname='http://www.nhc.noaa.gov/data/hurdat/hurdat2-1851-2016-041117.txt',
                         path=NULL,verbose=FALSE,...) {
  if(verbose) print('read.hurdat2')
  if(verbose) print(paste('file:',fname))
  if(!is.null(path) & !is.url(fname)) {
    fname <- file.path(path,fname)
  } else if (is.url(fname)) {
    destfile <- sub('http://','',fname)
    destfile <- sub('.html','.txt',destfile)
    destfile <- sub('.*/','',destfile)
    if (!is.null(path)) destfile <- file.path(path,destfile)
    if(!file.exists(destfile)) download.file(url=fname, destfile, method='auto', 
                                             quiet=FALSE, mode='w', cacheOK=TRUE)
    fname <- destfile
  }
  hurdat2 <- readLines(fname)
  n <- as.vector(sapply(hurdat2,nchar))
  i.storm <- which(n>80)
  i.name <- which(diff(n)>50)
  i.start <- min(i.name)
  x <- strsplit(hurdat2[i.storm],',')
  x <- do.call(rbind, x)
  x <- x[,-21]
  d <- dim(x)
  ri <- sapply(x[,3],function(x) gsub(' ','',x))
  ri[nchar(ri)==0] <- 0
  ri <- sapply(ri,function(x) switch(x,'0'=0,'L'=1,'W'=2,'P'=3,'I'=4,'C'=5,'S'=6,'G'=7,'T'=8))
  ri <- as.numeric(unlist(as.character(ri)))
  dict.ri <- seq(0,8)
  names(dict.ri) <- c('-','L','W','P','I','C','S','G','T')
  hu <- unlist(sapply(x[,4],function(x) gsub(' ','',x)))
  hu <- sapply(hu,function(x) switch(x,'TD'=1,'TS'=2,'HU'=3,'EX'=4,'SD'=5,'SS'=6,'LO'=7,'WV'=8,'DB'=9))
  hu <- as.numeric(unlist(as.character(hu)))
  dict.hu <- seq(1,9)
  names(dict.hu) <- c('TD','TS','HU','EX','SD','SS','LO','WV','DB')
  x[,3] <- ri
  x[,4] <- hu
  x[,5] <- sapply(x[,5],function(x) {
    y <- as.numeric(gsub('N|S','',x))
    if(grepl('S',x)) y <- -y
    return(y) })
  x[,6] <- sapply(x[,6],function(x) {
    y <- as.numeric(gsub('E|W','',x))
    if(grepl('W',x)) y <- -y
    return(y) })
  x <- as.numeric(x)
  dim(x) <- d
  h <- strsplit(hurdat2[i.name],',')
  h <- do.call(rbind, h)
  h <- apply(h,c(1,2),function(x) gsub(' ','',x))
  h[,1] <- sapply(h[,1],function(x) gsub('AL','',x))
  i.x <- sapply(i.storm, function(i) max(seq(length(i.name))[i.name<i]))
  x <- cbind(as.numeric(h[i.x,3]),x)
  x <- cbind(as.numeric(h[i.x,1]),x)
  colnames(x) <- c('trajectory','n','date','time','record.id','status',
                   'lat','lon','max.speed','pcent',
                   'r34.ne','r34.se','r34.sw','r34.nw',
                   'r50.ne','r50.se','r50.sw','r50.nw',
                   'r64.ne','r64.se','r64.sw','r64.nw')
  if(verbose) print(colnames(x))
  dict.names <- h[,2]
  names(dict.names) <- factor2numeric(h[,1])
  x <- as.events(x,longname='tropical storms and hurricanes',
                 param='storm tracks',method='HURDAT 2nd generation',
                 src='National Hurricane Center (NHC) revised Atlantic hurricane database (HURDAT2)',
                 reference='Landsea, C. W., et al. 2004: The Atlantic hurricane database re-analysis project: Documentation for the 1851-1910 alterations and additions to the HURDAT database. Hurricanes and Typhoons: Past, Present and
Future, R. J. Murname and K.-B. Liu, Eds., Columbia University Press, 177-221.',
                 file=fname,url=fname,verbose=verbose)
  attr(x,'name') <- dict.names
  attr(x,'record.id') <- dict.ri
  attr(x,'status') <- dict.names
  invisible(x)
}



# read.otto <- function(fname,path=NULL,progress=TRUE,verbose=FALSE) {
#   #path <- '~/Dropbox/stormtracks/Otto'
#   #fname <- 'trkdat_pmsl_2001-2001.cfsr'
#   if(verbose) print('read.otto')
#   if(verbose) print(paste('file:',fname))
#   if(!is.null(path)) fname <- file.path(path,fname)
#   if(verbose) print('read data')
#   otto <- readLines(fname)
#   n <- as.vector(sapply(otto,nchar))
#   h <- hist(n,sort(c(unique(n)-0.5,unique(n)+0.5)),plot=FALSE)
#   n.storm <- h$mids[which.max(h$counts)]
#   i.header <- which(n==n.storm)[1]
#   i.storm <- which(n==n.storm & c(10,diff(n))<=0)
#   i.track <- which(diff(n)==diff(n)[i.header-1])-1
#   header <- unlist(strsplit(otto[i.header],'\\s+'))
#   header <- header[nzchar(header)]
#   header[header=='x'] <- 'lon'
#   header[header=='y'] <- 'lat'
#   header[header=='t'] <- 'timestep'
#   header[header=='da'] <- 'date'
#   header[header=='hr'] <- 'time'
#   x <- strsplit(otto[i.storm],'\\s+')
#   x <- lapply(x,function(x) as.numeric(x[nzchar(x)]))
#   x <- do.call(rbind, x)
#   if(verbose) print('arange data')
#   if (progress) pb <- txtProgressBar(style=3)
#   if (progress) setTxtProgressBar(pb,0/(length(i.track)))
#   x.track <- as.numeric(gsub(':.*','',gsub('.*Track','',otto[i.track])))
#   n.track <- sapply(2:length(i.track),function(i) {
#                     if(progress) setTxtProgressBar(pb,i/(length(i.track)))
#                     n.i <- sum(i.storm>i.track[i-1] & i.storm<i.track[i])
#                     return(n.i) })
#   n.track <- c(n.track,sum(i.storm>i.track[length(i.track)]))
#   tracks <- unlist(sapply(seq_along(x.track),function(i) rep(x.track[i],n.track[i])))
#   y <- cbind(tracks,x)
#   colnames(y) <- c('trajectory',header)
#   if(verbose) print('organize date and time')
#   da <- y[,'date']
#   da[y[,'date']<2E5] <- da[y[,'date']<2E5] + 2E7
#   da[y[,'date']>=2E5] <- da[y[,'date']>=2E5] + 1.9E7
#   hr <- y[,'time']*1E-2
#   y[,'date'] <- da
#   y[,'time'] <- hr
#   z <- as.events(y,longname='sub-tropical cyclones',
#                  param='storm tracks',method='Otto...',
#                  src='FMI',
#                  reference='Otto',
#                  file=fname,url=fname,verbose=verbose)
#   return(z)
# }

