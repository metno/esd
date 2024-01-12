## Global mean temperature: models & observations

## Activate the esd-packages
#library(esd)
## If the jpeg-package is not installed, then install it.
install.jpeg <- ("jpeg" %in% rownames(installed.packages()) == FALSE)
if (install.jpeg) install.packages('jpeg')
library(jpeg)

url.earth="http://eoimages.gsfc.nasa.gov/images/imagerecords/36000/36019/AS8-16-2593_lrg.jpg"
image.earth <- 'earth.jpg'
if (!file.exists(image.earth)) {
  download.file(url.earth,image.earth)
}
img.earth <- readJPEG(image.earth)

## Function used for extracting one calendar month from list objects
Mon <- function(x,mon='Sep') {y <- subset(x,it=mon); return(y)}

## Read the CSV data downloaded from Copernicus C3S and organise the data:
HadCRUT5 <- read.csv("https://climate.metoffice.cloud/formatted_data/gmt_HadCRUT5.csv")
NOAA <- read.csv("https://climate.metoffice.cloud/formatted_data/gmt_NOAAGlobalTemp.csv")
gistemp <- read.csv("https://climate.metoffice.cloud/formatted_data/gmt_GISTEMP.csv")
era5 <- read.csv("https://climate.metoffice.cloud/formatted_data/gmt_ERA5.csv")
jra55 <- read.csv("https://climate.metoffice.cloud/formatted_data/gmt_JRA-55.csv")
berkeley <- read.csv("https://climate.metoffice.cloud/formatted_data/gmt_Berkeley%20Earth.csv")

## CMIP6 GCMs
#load('~/R/ssp245.mon.gmst.rda')
data('ssp245.mon.gmst')
ssp245 <- lapply(results,annual)
ns <- length(ssp245); nt <- length(ssp245[[1]]); t <- year(ssp245[[1]])
ssp245 <- unlist(ssp245); dim(ssp245) <- c(nt,ns)
ssp245 <- zoo(ssp245,order.by=t)

## Set plotting options, such ascolours
col <- c('black','red','blue','darkgreen','cyan','purple')
lwd <- c(5,rep(1,4))
it <- c(1850,1900)

par(bty="n",mar=rep(0,4),bg="black",col.main="white",col.lab="white",col.axis="white",col='white')
plot(c(1850,2030),c(-0.5,2),type="n",xlab="",ylab="")
rasterImage(img.earth,1950,-0.7,2050,1.3)

## Plot the GMST
par(new=TRUE,mar=c(7,4,4,1),las=1)
plot(subset(anomaly(ssp245,ref=1850:1900),it=c(1900,2050)),col=rgb(1,0.5,0.5,0.1),
     lwd=3,plot.type = 'single',xlab='',ylab=expression(degree*C))
legend(1910,3.5,c('ERA5','HadCUT5','CMIP6 SSP245'),col=c('white','white','pink'),
       lty=1,lwd=c(3,3,1),text.col='white')
for (i in seq(0.5,2,by=0.01)) {
  points(era5$Year,era5$ERA5..degC.,col=rgb(i/2,i/2,i/2,0.5),pch=19,lty=3,cex=1/(i+0.2))
  points(HadCRUT5$Year,HadCRUT5$HadCRUT5..degC.,col=rgb(i/2,i/2,i/2,0.5),pch=19,lty=3,cex=1/(i+0.2))
}
grid(col='white')
