
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

