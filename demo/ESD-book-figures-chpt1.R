## ESD book 2nd edition
## R-code for figures in Chapter 1.
## Rasmus Benestad, 2020-06-08

print("Figure 1.2")
library(esd)
## Fetch data on the mean sea-level pressure: select winter mean SLP
slp <- subset(as.4seasons(slp.NCEP()),it='djf')
## Fetch the NAO-index: select the winter mean
nao <- subset(as.4seasons(NAO()),it='djf')
## Show maps of correlation
corfield(nao,slp,new=FALSE)
  

z <- retrieve('~/Downloads/hgt.sfc.nc')
coredata(map(z,plot=FALSE)) -> Z
Z[Z <= 100] <- NA
map(Z)
