## Arctic amplification: plot trends in temperature at different pressure levels

library(esd)

layout(matrix(1:9,3,3),widths=rep(1,3),heights = rep(1,3))

for (plev in c(1000,850,700,500,250,100,30,1)) {
  ta <- annual(retrieve('~/Downloads/ERA5_TA_p-levels_mon.nc',lev=plev))
  map(ta,FUN='trend')
}

