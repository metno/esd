## R-script to read the annoying format of climate data from eKlima.
eklima <- function(fname) {
  test <- readLines(fname)
  Nnc <- as.numeric(table(nchar(test)))
  inc <- as.numeric(rownames(table(nchar(test))))
  nnc <- inc[is.element(Nnc,max(Nnc))]
  data <- test[is.element(nchar(test),nnc)]
  header <- test[grep('Date',test)][1]
  stid <- as.numeric(substr(data[1],1,11))
  t <- as.Date(paste(substr(data,18,21),substr(data,15,16),substr(data,12,13),sep='-'))
  varid <- strsplit(header[1],split=' ')[[1]]
  varid <- varid[!is.element(varid,c('','St.no','Date'))]
  for (i in 1:length(varid)) {
    ic <- regexpr(varid[i],header[1])[1]-2
    x <- as.numeric(substr(data,ic,ic+length(varid[i])+2))
    if (i==1) X <- x else X <- rbind(X,x)
  }
  z <- zoo(t(X),order.by=t)
  attr(z,'variable') <- varid
  ii <- grep('Stnr Name',test)+1
  attr(z,'longitude') <- as.numeric(substr(test[ii],70,80))
  attr(z,'latitude') <- as.numeric(substr(test[ii],62,70))
  attr(z,'altitude') <- as.numeric(substr(test[ii],50,60))
  iu <- grep('Unit',test) 
  attr(z,'unit') <- c(as.numeric(substr(test[iu+c(1:length(varid))],27,35)))
  attr(z,'longname') <- c(as.numeric(substr(test[iu+c(1:length(varid))],7,25)))
  attr(z,'location') <- rep(substr(test[ii],8,23),length(varid))
  class(z) <- c('station','daily','zoo')
  invisible(z)
  }
