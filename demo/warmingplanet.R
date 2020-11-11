## Spinning globe with model topography

library(esd)
library(animation)

frame <- function(x,it) {
  cat('.')
  par(bg='black',col.axis='white',col.lab='white',col.main='white')
  y <- subset(x,it=it)
  i <- it %% 360
  j <- 0   # it %/% 120
  map(y,projection='sphere',main=year(y),lonR=-i,latR=j,style='night',
      colbar=list(breaks=seq(-6,6,by=0.1)),new=FALSE)
}

if (!file.exists('tas.nc'))
  download.file('http://climexp.knmi.nl/CMIP5/monthly/tas/tas_Amon_ens_rcp45_000.nc','tas.nc')
tas <- retrieve('tas.nc')
tas <- annual(tas)
tas <- anomaly(tas,ref=1961:1990)
n <- length(index(tas))

print('MP4:')
saveVideo({for (i in 1:n) frame(tas,i)},
          video.name="Warmingplanet.mp4",interval=0.1)

print('GIF:')
saveGIF({for (i in 1:n) frame(tas,i)},
        movie.name="Warmingplanet.gif",interval=0.03)