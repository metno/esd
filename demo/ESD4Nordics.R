## Demo script for accessing and visualising netCDF files with multi-model CMIP6 ensemble mean results

url <- 'https://thredds.met.no/thredds/dodsC/metusers/rasmusb/'
x <- retrieve(file.path(url,'mu_Ayear_DSEns_Nordics_1850-2100_ssp370_.nc'))
par(mfcol=c(2,2))
map(x,main='Ensemble mean (SSP370)')
y <- retrieve(param='ens_sd_mu',file.path(url,'mu_Ayear_DSEns_Nordics_1850-2100_ssp370_.nc'))
map(y,main='Ensemble spread (SSP370)')
z <- retrieve(param='shapiro_test_normal',file.path(url,'mu_Ayear_DSEns_Nordics_1850-2100_ssp370_.nc'))
map(z,main='Test ensemble spread for normal distribution')
w <- retrieve(param='ks_test_uniform',file.path(url,'mu_Ayear_DSEns_Nordics_1850-2100_ssp370_.nc'))
map(w,main='Test ensemble statistics against observations')
              