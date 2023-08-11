## Northeast and North Central Pacific hurricane database (HURDAT2) 1949-2022
## Ref. Landsea, C. W. and J. L. Franklin, 2013: Atlantic Hurricane Database 
## Uncertainty and Presentation of a New Database Format. Mon. Wea. Rev., 141, 3576-3592.

## Rasmus Benestad, 2023-08-08

library(esd)
url <- 'https://www.nhc.noaa.gov/data/hurdat/hurdat2-nepac-1949-2022-050423.txt'
rawdata <- readLines(url)

rows <- nchar(rawdata)
print(table(rows))

storm <- (1:length(rows))[rows==37]

hurdat2.stormtracks <- list()
attr(hurdat2.stormtracks,'url') <- url
attr(hurdat2.stormtracks,'time') <- Sys.time()
lon.rng <- c(); lat.rng <- c(); tim.rng <- c()
maxdist <- 0; maxlon <- 0

for (is in 1:length(storm)) {
  storminfo <- unlist(strsplit(rawdata[storm[is]],','))
  label <- paste(gsub(' ','',storminfo[2]),gsub('\t','',storminfo[1]),sep='.')
  nt <- as.numeric(storminfo[3])
  x <- strsplit(rawdata[(storm[is]+1):(storm[is]+nt)],',')
  date <- as.Date(unlist(lapply(x,function(x) paste(substr(x[1],1,4),substr(x[1],5,6),substr(x[1],7,8),sep='-'))))
  time <- unlist(lapply(x,function(x) x[2]))
  print(paste(label,date[1]))
  category <- unlist(lapply(x,function(x) x[4]))
  lat <- unlist(lapply(x,function(x) x[5]))
  NH <- grep('N',lat); SH <- grep('S',lat)
  lat <- as.numeric(sub('S','',sub('N','',lat)))
  if (length(SH)>0) lat[SH] <- -lat[SH]
  lon <- unlist(lapply(x,function(x) x[6]))
  WH <- grep('W',lon); EH <- grep('E',lat)
  lon <- as.numeric(sub('E','',sub('W','',lon)))
  if (length(WH)>0) lon[WH] <- -lon[WH]
  lon[lon > 0] <- lon[lon > 0] - 360
  stormtrack <- data.frame(date=date,time=time,lon=lon,lat=lat)
  hurdat2.stormtracks[[label]] <- stormtrack
  #print(paste0(paste(range(lon,collapse='-')),'E',range(lat,collapse='-'),'N'))
  lon.rng <- range(c(lon.rng,lon))
  lat.rng <- range(c(lat.rng,lat))
  tim.rng <- as.Date(range(c(tim.rng,date)))
  dist <- 0
  if (nt > 1) for (j in 2:nt) dist <- dist + distAB(lon[j-1],lat[j-1],lon[j],lat[j]) else dist <- 0
  maxdist <- round(max(maxdist,dist,na.rm=TRUE)/1000)
  maxlon <- max(c(maxlon,diff(range(lon))))
}

print(paste0('Space covered: ',paste(lon.rng,collapse='-'),'E, ',paste(lat.rng,collapse='-'),'N, ',
             'over the period ',paste(tim.rng,collapse=' to '),' number of storms= ',length(hurdat2.stormtracks),
             ', maximum distance=',maxdist,'km',', maximum longitude span= ',maxlon))
save(hurdat2.stormtracks,file='hurdat2.stormtracks.rda')

plot(lon.rng,lat.rng,type='n',main='Pacific tropical cyclones (HURDAT2)',xlab='',ylab='',
     sub=paste(tim.rng,collapse='to'))
data(geoborders)
lines(geoborders)
lines(geoborders$x - 360, geoborders$y)
grid()
lapply(hurdat2.stormtracks,function(x) {if (min(x$lon) < -200) col <- rgb(0,0,0) else 
  if (min(x$lon) < -180) col <- rgb(0,0,0.4) else col <- rgb(0.3,0.3,0.3,0.3); 
                                       lines(x$lon,x$lat,col=col)})

## 2023 Hurricane Dora
## <https://rammb-data.cira.colostate.edu/tc_realtime/storm.asp?storm_identifier=ep052023>
if (file.exists('~/Downloads/TC-Dora-2023.txt')) {
  dora <- read.table('~/Downloads/TC-Dora-2023.txt',header=TRUE)
  lines(dora$Longitude,dora$Latitude,col='red',lwd=2)
  nd <- length(dora$Longitude)
  text(dora$Longitude[nd],dora$Latitude[nd],'Dora (2023)',pos=1,col='red')
  duration.dora <- diff(range(as.Date(dora$Synoptic)))
  londiff.dora <- diff(range(dora$Longitude))
  print(paste('Dora has lived for',duration.dora,'days and trasversed',londiff.dora,'longitudes'))
} else {duration.dora <- NULL; londiff.dora <- NULL}

breaks <- seq(year(min(tim.rng)),year(max(tim.rng)),by=1)
nino3.4 <- subset(annual(anomaly(NINO3.4(),ref=1981:2010)),it=c(year(min(tim.rng)),year(max(tim.rng))))
col <- rep('grey',length(breaks))
col[nino3.4 > 0.5] <- 'red'
col[nino3.4 < -0.5] <- 'blue'
ntc.yr <- unlist(lapply(hurdat2.stormtracks,function(x) year(x$date[1])))
hist(ntc.yr,breaks=breaks,col=col,
     main='Number ot tropical cyclones per year in the western/central Pacific',xlab='year',ylab='count')
grid()


duration <- unlist(lapply(hurdat2.stormtracks,function(x) diff(range(x$date))))
hist(duration,main='Duration of TCs',xlab='days',ylab='count')
grid()
if (!is.null(duration.dora)) lines(rep(duration.dora,2),c(0,500),lwd=3,col="red")

londiff <- unlist(lapply(hurdat2.stormtracks,function(x) diff(range(x$lon))))
hist(londiff,main='Longitude travelled by TCs',xlab='days',ylab='count')
grid()
if (!is.null(londiff.dora)) lines(rep(londiff.dora,2),c(0,500),lwd=3,col="red")

