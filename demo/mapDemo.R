## Demo for testing map() of various sort
## map(), showmaps(), show(), display(), plot(), vis()

## showmaps()
t2m <- t2m.NCEP()
period <- paste(range(year(t2m)),collapse=' - ')
attr(t2m,'sub') <- paste0('Average over: ',period,'; source: NCEP')
showmaps(t2m,FUN='mean',colbar=list(breaks=seq(-60,20,by=5),show=FALSE))

## Ordinary lon-lat map
map(t2m,colbar=list(pal='t2m.IPCC')) -> mymap

## Make the same map with same colour palette but with a northern polar seterographic projection
map(mymap,colbar=attr(mymap,'colbar'),projection='np')

## This demo is used to test the various uses of map()

## A spherical projection
par(bg='black')
map(t2m,projection='sphere',style='night',colbar=list(show=FALSE))

par(bg='white')
## map of stations
data("t2m.NORDKLIM")
map(t2m.NORDKLIM,FUN='mean')

## Map EOFs
data("eof.t2m.NorESM.M")
map(eof.t2m.NorESM.M)
plot(eof.t2m.NorESM.M)

## Plot CCA
data("eof.slp.NCEP")
data("eof.t2m.NCEP")
cca <- CCA(eof.slp.NCEP,eof.t2m.NCEP)
plot(cca)

## Plot DS
data(Oslo)
ds <- DS(Oslo,eof.t2m.NCEP)
plot(ds)

## Map storm tracks and events

