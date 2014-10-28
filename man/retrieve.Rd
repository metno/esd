\name{retrieve}
\alias{retrieve}
\title{Retrieve field data from a netcdf file.}
\description{Retrieve data from a netcdf file and return a zoo field object with attributes. 
  \code{\link{retrieve}}.
}
\usage{
%retrieve.default(...)
retrieve.default(ncfile,...)
retrieve.ncdf(ncfile = ncfile, path = path , param = "auto",
               lon = NULL, lat = NULL, lev = NULL, time = NULL,
               miss2na = TRUE, greenwich = TRUE , ncdf.check = TRUE,
               verbose = FALSE , plot = TRUE)
retrieve.ncdf4(ncfile = ncfile, path = path , param = "auto",
               lon = NULL, lat = NULL, lev = NULL, time = NULL,
               miss2na = TRUE, greenwich = TRUE , ncdf.check = TRUE,
               verbose = FALSE , plot = TRUE) 

summary.ncdf4(ncfile, verbose=FALSE)

check.ncdf4(ncid = ncid, param = "auto" , verbose = FALSE)
check.ncdf(ncid = ncid, param = "auto" , verbose = FALSE)
}

\arguments{
  \item{ncfile}{A character string of full path netcdf file name
  (include the path if necessary) or any object of class 'ncdf' or 'ncdf4'.}
  \item{ncid}{An object of class 'ncdf' or 'ncdf4'.}
  \item{lon}{Numeric value of longitude for the reference point (in decimal degrees East) or a vector containing the range of longitude values in the form of c(lon.min,lon.max)}
  \item{lat}{Numeric value of latitude for the reference point (in decimal degrees North) or a vector containing the range of latitude values in the form of c(lat.min,lat.max)}
  \item{lev}{Numeric value of pressure levels or a vector containing the range of pressure level values in the form of c(lev.min,lev.max) }
  \item{time}{Numerical year values or date values of time or a vector containing the range of values in the form of c(time.min,time.max). Date format should be in the form of "YYYY-MM-DD".}
  \item{param}{Parameter or element type. There are several core parameters or elements as well as a number of additional parameters. The parameters or elements are:
           auto         = automatic selection.
           prcp, pr     = Precipitation (mm)
   	   tas, tavg    = 2m-surface temperature (in degrees Celcius)
           tmax, tasmax = Maximum temperature (in degrees Celcius)
           tmin, tasmin = Minimum temperature (in degrees Celcius)
           }
  \item{plot}{Logical value. if, TRUE provides a map.}
  \item{greenwich}{Logical value. If FALSE, convert longitudes to -180E/180E or centre maps on Greenwhich meridian (0 deg E).} 
  \item{ncdf.check}{Logical value. If TRUE, performs a quick check of the ncfile contents} 
  \item{miss2na}{Logical value. If TRUE missing values are converted to "NA"}
  \item{verbose}{Logical value defaulting to FALSE. If FALSE, do not
  display comments (silent mode). If TRUE, displays extra information on progress.}
}
\value{A "zoo" "field" object with additional attributes used for further processing.}

\seealso{
       \code{\link{test.retrieve.ncdf4}}.
}

\author{A. Mezghani}
\examples{
#\dontrun{
# Download air surface temperature (tas) for the 'NorESM1-ME' model
# output prepared for 'CMIP5 RCP4.5' and for run 'r1i1p1' from the climate
# explorer web portal (http://climexp.knmi.nl) and store the file into the temporary folder '/tmp' (Size ~96Mb)
download.file(url="http://climexp.knmi.nl/CMIP5/monthly/tas/tas_Amon_NorESM1-ME_rcp45_000.nc",
 destfile="/tmp/tas_Amon_NorESM1-ME_rcp45_000.nc", method = "auto", quiet = FALSE,mode = "w",cacheOK = TRUE)
# Retrieve the data into "gcm" object
gcm <-retrieve.ncdf4(ncfile="/tmp/tas_Amon_NorESM1-ME_rcp45_000.nc",param="tas",plot=TRUE)
# Download the air surface temperature (tas) for RCP 4.5 scenarios and
# NorESM1-ME model from the climate explorer and store it in destfile. 
# Compute the anomalies
gcm.a <- as.anomaly(gcm,ref=c(1960:2001))
map(gcm.a,projection="sphere")
#}
}
\keyword{?}
