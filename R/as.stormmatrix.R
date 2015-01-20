as.stormmatrix <- function(x,...) UseMethod("as.stormmatrix")

as.stormmatrix.default <- function(x) {
  x$cyclone <- as.numeric(as.character(x$cyclone))
  x$Date <- as.numeric(as.character(x$Date))
  x$Year <- as.numeric(as.character(x$Year))
  x$Lon <- as.numeric(as.character(x$Lon))
  x$Lat <- as.numeric(as.character(x$Lat))
  x$Pressure <- as.numeric(as.character(x$Pressure))
  x <- x[x$Code99<90,]
  if (x$Code99[1]>9) {
    m <- paste("M",as.character(x$Code99[1]),sep="")
  } else m <- paste("M0",as.character(x$Code99[1]),sep="")
  x$Method <- m
  M <- imilast2matrix(x)
  invisible(M)
}

imilast2matrix <- function(x) {
  fn <- function(x) approx(x,n=10)$y
  x <- data.frame(x)
  aggregate(x$Lon, list(x$cyclone), fn)$x -> lon
  aggregate(x$Lat, list(x$cyclone), fn)$x -> lat
  aggregate(x$Pressure, list(x$cyclone), fn)$x -> slp
  aggregate(x$Date, list(x$cyclone), function(x) x[1])$x -> t1
  aggregate(x$Date, list(x$cyclone), function(x) x[length(x)])$x -> t2
  aggregate(x$Date, list(x$cyclone), length)$x -> n
  colnames(lon) <- rep('lon',10)
  colnames(lat) <- rep('lat',10)
  colnames(slp) <- rep('slp',10)
  #X <- list(lon=lon,lat=lat,slp=slp,start=t1,end=t2,n=n)
  X <- cbind(lon=lon,lat=lat,slp=slp,start=t1,end=t2,n=n)
  attr(X, "location")= NA
  attr(X, "variable")= 'storm tracks'
  attr(X, "unit")= NA
  attr(X, "longitude")= NA
  attr(X, "latitude")= NA
  attr(X, "altitude")= 0
  attr(X, "country")= NA
  attr(X, "longname")= "mid-latitude storm tracks"
  attr(X, "station_id")= NA
  attr(X, "quality")= NA
  attr(X, "calendar")= "gregorian"
  attr(X, "source")= "IMILAST"
  attr(X, "URL")= "http://journals.ametsoc.org/doi/abs/10.1175/BAMS-D-11-00154.1"
  attr(X, "type")= "analysis"
  attr(X, "aspect")= "interpolated"
  attr(X, "reference")= "Neu, et al. , 2013: IMILAST: A Community Effort to Intercompare Extratropical Cyclone Detection and Tracking Algorithms. Bull. Amer. Meteor. Soc., 94, 529–547."
  attr(X, "info")= NA
  attr(X, "method")= x$Method
  attr(X, "history")= history.stamp()
  class(X) <- 'stormmatrix'
  invisible(X)
}

#fname="/vol/fou/klima/IMILAST/ERAinterim_1.5_NH_M03_19890101_20090331_ST.txt"
#m03 <- read.fwf(fname,width=c(2,7,4,11,5,5,5,4,7,7,10),
#        col.names=c("Code99","cyclone","timeStep","Date","Year",
#        "Month","Day","Time","Lon","Lat","Pressure"))
#x <- as.stormtrack(m03)
