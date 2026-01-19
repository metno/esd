## Demo script showing how to make use of downscaled wet-day frequency and mean precipitation
## and WG to produce daily rainfall data

## Retrieve gridded key annual rainfall statistics for the Nordic countries from thredds/OpenDap
url <- 'https://thredds.met.no/thredds/dodsC/metusers/rasmusb/'
mMU <- retrieve(file.path(url,'mu_Ayear_DSEns_Nordics_1850-2100_ssp370_.nc'))
mFW <- retrieve(file.path(url,'fw_Ayear_DSEns_Nordics_1850-2100_ssp370_.nc'))
sMU <- retrieve(param='ens_sd_mu',file.path(url,'mu_Ayear_DSEns_Nordics_1850-2100_ssp370_.nc'))
sFW <- retrieve(param='ens_sd_fw',file.path(url,'fw_Ayear_DSEns_Nordics_1850-2100_ssp370_.nc'))
## Use bjornholt (near Oslo, Norway) as an example
data(bjornholt)
## Extract daily annual statistics for the coordinates corresponding to selected site
## using bi-linear interpolation: the ensemble mean
mmu <- regrid(mMU,is=bjornholt)
mfw <- regrid(mFW,is=bjornholt)
## There was some missing data in fw
ok <- is.finite(mfw)
coredata(mfw) <- approx(year(mfw[ok]),coredata(mfw)[ok],xout = year(mfw))$y
## The ensemble spread
smu <- regrid(sMU,is=bjornholt)
sfw <- regrid(sFW,is=bjornholt)
ok <- is.finite(sfw)
coredata(sfw) <- approx(year(sfw[ok]),coredata(sfw)[ok],xout = year(sfw))$y
## Create annual statistics based on mean and standard deviation
mu <- zoo(rnorm(length(mmu),mean=mmu,sd=smu),order.by=year(mmu))
fw <- zoo(rnorm(length(mfw),mean=mfw,sd=sfw),order.by=year(mfw))
## Here
z <- WG(bjornholt,mu=mu,fw=fw)
yz <- combine.stations(bjornholt,z)
plot(yz)
