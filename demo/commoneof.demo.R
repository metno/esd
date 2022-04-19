## Demonstration of empirical orthogonal functions (EOFs) and common EOFs
## @RasmusBenestad 2021-12-09

## Activate the esd-package - open-source code freely available from https://github.com/metno/esd
library(esd)

## Get surface air temperature from two reanalysis data for an area over Africa
## NCEP reanalysis 1:
x <- retrieve('~/Downloads/air.mon.mean.nc',lon=c(-10,40),lat=c(-20,20))
## Take the annual mean temperature
x <- annual(x)
## Estimate EOFs of the 
eof <- EOF(x)
## Plot the EOFs
plot(eof)

## Get corresponding data from the ERA5 reanalysis for the same region
y <- retrieve('~/Downloads/ERA5_t2m_year.nc',lon=c(-10,40),lat=c(-20,20))
## Ensure that timestamp of y is year:
index(y) <- year(y)

## Combine the two reanalyses: take anomalies and regrid y to the grid of x
xy <- combine(x,y)
## Common EOFs are merely EOFs performed on the joint data matrix of anomalies
ceof <- EOF(xy)
## The common EOFs are great for comparing different reanalyses:
plot(ceof)
