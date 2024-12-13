station.YR <- function(is,url="https://api.met.no/weatherapi/locationforecast/2.0/",verbose=TRUE) { 
# Load jsonlite package
library(jsonlite)

if (inherits(is,'station')) is <- list(lon=lon(is),lat=lat(is))
## If field, then need to repeat so that there is equal lons and lats
# URL
#url <- "https://api.met.no/weatherapi/locationforecast/2.0/compact?lat=60.10&lon=9.58"

URL <- paste0(url,'compact?lat=',lat,'&lon=',lon)
# Read JSON data from the URL
data <- fromJSON(url)

# View the data
View(data)
}