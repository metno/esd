#' @export
read.best.track <- function(url='http://www.usno.navy.mil/NOOC/nmfc-ph/RSS/jtwc/best_tracks/',domain='io',
                            start=1945,end=2014,n=20,verbose=TRUE) {
  ## domain = c('io','sh','wp')
  ## Format  	
  ## BASIN,CY, YYYYMMDDHH,TECHNUM, TECH,TAU, LatN/S,LonE/W, VMAX,MSLP
  ## TY,RAD, WINDCODE,RAD1, RAD2,RAD3, RAD4,RADP, RRP,MRD, GUSTS,EYE,
  ## SUBREGION,MAXSEAS, INITIALS,DIR, SPEED,STORMNAME, DEPTH,SEAS, SEASCODE,SEAS1, SEAS2,SEAS3, SEAS4 
  
  TC <- list(info='best-track cyclones')
  storms <- c()
  ii <- 1
  for (yr in start:end) {
    for (i in 1:30) {
      if (i < 10) ci <- paste(0,i,sep='') else ci <- as.character(i)
      urldata <- paste(url,yr,'/',yr,'s-b',domain,'/b',domain,ci,yr,'.txt',sep='')
      print(urldata)
      track <- try(read.table(urldata,sep=','),silent=TRUE)
      if (inherits(track,'try-error')) {
        urldata <- paste(url,yr,'/',yr,'s-b',domain,'/b',domain,ci,yr,'.dat',sep='')
        track <- try(read.table(urldata,sep=','),silent=TRUE)
      }
      if (!inherits(track,'try-error')) {
        tc <- paste(domain,yr,ii,sep='.')
        TC[[tc]] <- track
        storms <- c(storms,tc)
      }
    }
  }
  
  ## Create an esd trajectory-object from the tropical cyclone storm tracks
  ns <- length(storms)
  Y <- matrix(rep(NA,(3*n+3)*ns),ns,3*n+3)
  for (i in 1:ns) {
    z <- TC[[storms[i]]]
    south <- regexpr('S',z$V7) > 0
    west <- regexpr('W',z$V8) > 0
    Lon <- as.numeric(sub('W','',sub('E','',as.character(z$V8))))/10
    Lat <- as.numeric(sub('S','',sub('N','',as.character(z$V7))))/10
    if (south) Lat <- -Lat
    if (west) Lon <- -Lon
    np <- length(Lon)
    lons <- approx(1:np,Lon,1:n)$y
    lats <- approx(1:np,Lat,1:n)$y
    slps <- approx(1:np,z$V9,1:n)$y
    Y[i,] <- c(lons,lats,slps,start,end,n)
  }
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
  attr(Y, "source")= 'US NAVY JTWC'
  attr(Y, "URL")= url
  attr(Y, "unit")= NA
  attr(Y, "type")= "analysis"
  attr(Y, "aspect")= "interpolated"
  attr(Y, "reference")= ''
  attr(Y, "info")= 'http://www.usno.navy.mil/NOOC/nmfc-ph/RSS/jtwc/best_tracks/TC_bt_report.html'
  attr(Y, "method")= NA
  attr(Y,"lon") <- NA
  attr(Y,"lat") <- NA
  attr(Y,"alt") <- NA
  attr(Y,"cntr") <- NA
  attr(Y,"stid") <- NA
  attr(Y,'domain') <- switch(domain,
                             'io'='Northern Indian Ocean',
                             'sh'='Southern Hemisphere','wp'='Northwestern Pacific')
  attr(Y, "history")= history.stamp()
  attr(Y,'storm name') <- storms
  class(Y) <- c('trajectory','matrix')
  invisible(Y)
  return(Y)
}
