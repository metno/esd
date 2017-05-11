## Function that reads HURDAT2 data
## Rasmus Benestad. Dubrovnik 2017-05-09

## This support function organises the data into lists with names elements
## for easy processing.
hur2list <- function(x) {
  south <- regexpr('S',x[7]) > 0
  west <- regexpr('W',x[8]) > 0
  Lon <- as.numeric(sub('W','',sub('E','',x[8])))
  Lat <- as.numeric(sub('S','',sub('N','',x[7])))
  if (south) Lat <- -Lat
  if (west) Lon <- -Lon
  #yyyymmdd <- paste(substr(x[3],1,4),substr(x[3],5,6),
  #                  substr(x[3],7,8),sep='-')
  yyyymmdd <- x[3]
  ##print(yyyymmdd)
  if ((nchar(x[3])==8) & (as.numeric(substr(x[3],1,4)) > 1700) &
      (as.numeric(substr(x[3],1,4)) < 2100) & 
      (as.numeric(substr(x[3],5,6)) > 0) & (as.numeric(substr(x[3],5,6)) < 13) &
      (as.numeric(substr(x[3],7,8)) > 0) & (as.numeric(substr(x[3],7,8)) < 32) ) {
     hurdat <- list(cyclone.no=as.numeric(x[1]),name=x[2],date= yyyymmdd,
                    type=x[6],lon=Lon,lat=Lat,speed=as.numeric(x[9]),
                    slp=as.numeric(x[10]))
  } else {print(x[3]); print(yyyymmdd); browser()}
  return(hurdat)
}

TC.track <- function(x,n) {
  ## TC.track collects all elements in a list that belongs to one cyclone
  ## into one matrix row
  ix <- 1:length(x$lon)
  if (length(ix) != length (x$lon)) {print(ix); print(x$lon)}
  start <- min(as.numeric(as.character(x$date)))
  end <- max(as.numeric(as.character(x$date)))
  if ( (sum(is.finite(x$lon)) > 3) & (sum(is.finite(x$lat)) > 3) &
       (sum(is.finite(x$slp)) > 3) ) {
    lons <- approx(x=ix,y=as.numeric(as.character(x$lon)),xout=1:n)$y
    lats <- approx(x=ix,y=as.numeric(as.character(x$lat)),xout=1:n)$y
    slps <- approx(x=ix,y=as.numeric(as.character(x$slp)),xout=1:n)$y
    row <- c(lons,lats,slps,start,end,n)
  } else row <- rep(NA,3*n+3)
  return(row)
}

read.hurdat2 <- function(url='http://www.aoml.noaa.gov/hrd/hurdat/Data_Storm.html',
                         n=20,verbose=FALSE) {
  ## This is an awkward file read first the lines
  aoml <- readLines(url)
  ## Sniff out the URL to the data from this website - the url seems to
  ## change with new updates.
  urltest <- aoml[grep('HURDAT 2',aoml)]
  urldata <- substr(urltest,regexpr('href=',urltest)+6,regexpr('txt',urltest)+2)[1]
  urldata <- paste('http://www.aoml.noaa.gov/hrd/hurdat/',urldata,sep='')
  datach <- readLines(urldata)
  ## Short rows contain storm names - select these for keeping track of each storm
  inm <- (1:length(datach))[nchar(datach)==37]
  ## Get the storm names
  stormnames <- unlist(lapply(strsplit(datach[inm],','),function(x) x[2]))
  stormnames <- gsub(' ','',stormnames)
  ## Identify the rows with the storm track data
  ist <- (1:length(datach))[nchar(datach)==120]
  ## Add the storm number and name to the rows of storm track data
  ns <- length(stormnames)
  ## Add a last one to capture all storms
  inm <- c(inm,length(datach)+1)
  for (i in 1:ns) {
    ii <- inm[i]:(inm[i+1]-1)
    print(c(i,stormnames[i],inm[i],(inm[i+1]-1)))
    datach[ii] <- paste(i,stormnames[i],datach[ii],sep=',')
  }
  
  X <- strsplit(datach[ist],',')
  
  ## Check the data structure
  xl <- unlist(lapply(X,function(x) length(x)))
  ci <- unlist(lapply(X,function(x) as.numeric(x[1])))
  if (verbose) print(table(xl))
  
  ## Stringsplit to create a list with storm tracks - the different elements
  ## contain different times, and one storm track includes several elements
  
  X <- lapply(X,hur2list)
  
  ## Restructure the data into a data.frame
  d <- c(length(X),length(X[[1]]))
  nms <- names(X[[1]])
  X <- as.data.frame(t(matrix(unlist(X),d[2],d[1])))
  names(X) <- nms
  
  ## Combine the elements belonging to the same storm track and structure as
  ## an trajectory object (as.trajectory). The structure is a matrix where 
  ## each row is [lon,lat,slp,start,end,n] and n is the number of lon/lat/slp points.
  ## Aggregate on name and use approx for interpolating between end points.
  
  Y <- matrix(rep(NA,(3*n+3)*ns),ns,3*n+3)
  for (i in 1:ns) Y[i,] <- TC.track(subset(X,cyclone.no==i),n=n)
  Y[Y <= -999] <- NA
  ok <- is.finite(rowMeans(Y))
  Y <- Y[ok,]
  colnames(Y) <- c(rep('lon',n),rep('lat',n),rep('slp',n),'start','end','n')
  
  # add attributes to trajectory matrix X
  attr(Y, "location")= NA
  attr(Y, "variable")= 'storm tracks'
  attr(Y, "longname")= 'Tropical cyclone storm tracks'
  attr(Y, "quality")= NA
  attr(Y, "calendar")= "gregorian"
  attr(Y, "source")= 'NOAA/USA'
  attr(Y, "URL")= url
  attr(Y, "unit")= NA
  attr(Y, "type")= "analysis"
  attr(Y, "aspect")= "interpolated"
  attr(Y, "reference")= ''
  attr(Y, "info")= info
  attr(Y, "method")= 'HURDAT2'
  attr(Y,"lon") <- NA
  attr(Y,"lat") <- NA
  attr(Y,"alt") <- NA
  attr(Y,"cntr") <- NA
  attr(Y,"stid") <- NA
  attr(Y, "history")= history.stamp()
  attr(Y,'storm name') <- stormnames[ok]
  class(Y) <- c('trajectory','matrix')
  invisible(Y)
}
