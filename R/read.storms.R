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

