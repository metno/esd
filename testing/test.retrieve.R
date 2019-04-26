# test.retrieve.R

library(esd)

path.lustre <- "/vol/lustre"
# Daily or sub-hourly
file1 <- file.path(path.lustre,"storeA/users/kajsamp/Data/CMIP5/6hr/psl_6hrPlev_HadGEM2-ES_historical_r1i1p1_195012010600-195112010000.nc")
#file2 <- file.path(path.lustre,"storeB/users/abdelkaderm/CixPag/wrfmonthly_d01_2047-08.nc")
file3 <- file.path(path.lustre,"storeA/users/abdelkaderm/cmip5/day/psl_day_MPI-ESM-LR_rcp85_r1i1p1_21400101-21491231.nc")
file4 <- file.path(path.lustre,"storeB/users/kajsamp/Data/ERAINT/ERAINT-full-slp_2016.nc")
# Monthly
file5 <- file.path(path.lustre,"storeB/users/kajsamp/Data/NCEP/slp.mon.mean.nc")
file6 <- file.path(path.lustre,"storeA/users/kajsamp/Data/CMIP5/KNMI/GCM106.tas.rcp45.nc")
file7 <- file.path(path.lustre,"storeB/users/kajsamp/Data/EOBS/pp_0.25deg_reg_v13.1.nc")

file <- file3
x <- retrieve(file, verbose=TRUE)
range(index(x))
class(x)
class(index(x))
map(x)
xsub <- subset(x, it=NULL, is=list(lon=c(-30,30),lat=c(35,75)))
xm <- as.monthly(xsub)
xeof <- EOF(xsub)
plot(xeof)

