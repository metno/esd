## Test time of netCDF retrieval for stations:

library(esd)

t0 <- Sys.time()
y1 <- retrieve.station('~/OpenClimateData/data/t2m.metnod.nc',stid=c(18700,39040,44640))
t1 <- Sys.time()
y2 <- retrieve.station('~/t2m.metnod2.nc',stid=c(18700,39040,44640))
t2 <- Sys.time()

print(paste(t1-t0,'for old ordering (time-space) and ',t2-t1,'for new ordering (space-time)'))