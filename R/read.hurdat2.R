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