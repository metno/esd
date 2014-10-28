\name{station}
\alias{select.station}
\alias{station}
\alias{station.default}
\alias{station.ecad}
\alias{station.nacd}
\alias{station.narp}
\alias{station.nordklim}
\alias{station.metnod}
\alias{station.metnom}
\alias{station.ghcnd}
\alias{station.ghcnm}

\title{Retrieve meta data and data from observational weather stations.}

\description{Retrieve station record from a given data source.}

\usage{

# Retrieve meta data 
select.station(loc=NULL , param = NULL,  stid = NULL ,lon = NULL,
               lat = NULL, alt = NULL, cntr = NULL, src = NULL ,
               it = NULL , nmin = NULL , verbose=FALSE) 

# Retrieve data
station(stid=NULL,...)

station.default(loc=NULL, param="t2m", src = c("GHCNM", "ECAD", "GHCND",
                "NACD","NARP","NORDKLIMA"), stid=NULL, lon=NULL, lat=NULL,
                alt=NULL, cntr=NULL, it= NULL,nmin=30,
                path=NULL, plot=FALSE,verbose=FALSE)

# Monthly weather stations
station.nordklim(...)
station.narp(...)
station.nacd(...)
station.ghcnm(...,path="data.GHCNM")
station.metnom(...)

#Daily weather stations
station.ecad(..., path="data.ECAD")
station.ghcnd(..., path="data.GHCND")
station.metnod(...)

}

\arguments{
  \item{loc}{A string of characters as the name of the location
    (weather/climate station) or an object of class "stationmeta".}
  \item{param}{Parameter or element type or variable identifier. There are several core parameters or elements as well as a number of additional parameters. The parameters or elements are:
    prcp, pr     = Precipitation (mm)
    tas, tavg    = 2m-surface temperature (in degrees Celcius)
    tmax, tasmax = Maximum temperature (in degrees Celcius)
    tmin, tasmin = Minimum temperature (in degrees Celcius)
  }
  \item{src}{Source: limit the downscaling to a specific data set ("NARP",
    "NACD", "NORDKLIMA", "GHCNM", "METNOM", "ECAD", "GHCND" and "METNOD")}
  \item{stid}{A string of characters as an identifier of the weather/climate station.}
  \item{lon}{Numeric value of longitude (in decimal degrees East) for
    the reference point (e.g. weather station) as a single value or a
    vector containing the range of longitude values in the form of c(lon.min,lon.max)}
  \item{lat}{Numeric value of latitude for the reference point (in
    decimal degrees North) or a vector containing the range of latitude
    values in the form of c(lat.min,lat.max)}
  \item{alt}{Numeric value of altitude (in meters a.s.l.) used for selection. Positive value,
    select all stations above this altitude; for negative values, select
    all stations below this latitude.}
  \item{cntr}{A string or a vector of strings of the full name of the
    country: Select the stations from a specified country or a
    set of countries.}
  \item{it}{A single integer or a vector of integers or Dates. An integer in the
  range of [1:12] for months, an integer of 4 digits for years
  (e.g. 2014), or a vector of Dates in the form "2014-01-01").}
  \item{nmin}{Select only stations with at least nmin number of
  years, months or days depending on the class of object x (e.g. 30 years).}
  \item{plot}{Logical value. If, TRUE provides a plot.}
  \item{verbose}{Logical value defaulting to FALSE. If FALSE, do not
    display comments (silent mode). If TRUE, displays extra information on
    progress.}
  \item{path}{The path where the data are stored. Can be a symbolic link.}
}
\value{A time series of "zoo" "station" class with additional attributes
  used for further processing.}

\seealso{
  \code{\link{station.meta}} and \code{\link{map.station}}.
}

\author{A. Mezghani}
\examples{
# \dontrun{
# Get daily and monthly mean temperature for "Oslo" station ("18700") from METNO data source
t2m.dly <- station.metnod(stid="18700",param="t2m")
t2m.mon <- station.metnom(stid="18700",param="t2m")

# Get daily data from the ECA&D data source:
# If called for the first time, the script will download a huge chunk of
# data and store it locally.
# select meta for "De Bilt" station into ss,
ss <- select.station(loc = "de bilt",param="t2m",src="ECAD")
# Retrieve the data from the local directory specified in path based on
# previous selected station 
t2m.dly <- station.ecad(loc=ss,path="data.ECAD")
# or directly retrieve the data without a prior selection 
t2m.dly <- station.ecad(loc = "oslo - blindern",param="t2m",path="data.ECAD")
plot(t2m.dly)
# Aggregate to monthly and annual mean temperature values and plot the results
t2m.mon <- as.monthly(t2m.dly, FUN="mean") ; plot(t2m.mon)
t2m.ann <- as.annual(t2m.mon, FUN = "mean") ; plot(t2m.ann)
# specify one station from ECAD, and this time get daily mean precipitation
precip.dly <- station.ecad(loc="Oxford",param="precip") ; plot(precip.dly)
# Aggregate to annual accumulated precipitation values and plot the result
precip.ann <- as.annual(precip.dly,FUN="sum") ; plot(precip.ann)

# Get daily data from the GHCND data source
# select the meta data of a subset of stations in ss and map the geographical location of the selected stations with a zoom on Norway.
ss <- select.station(cntr="NORWAY",param="precip",src="GHCND")
map(ss, xlim = c(-10,30) , ylim = c(50,70))
# Select a subset of stations across Norway with a minimum number of
# 130 years using "GHCND" as a data source, retrieve the data and show its
# structure.
ss <- select.station(cntr="NORWAY",param="precip",src="GHCND",nmin=130)
y <- station.ghcnd(loc=ss , path="data.GHCND",plot=TRUE)
str(y)
# Subselect one station and display the geographical location of both selected stations and highlight the subselected station (is=2).
y1 <- subset(y,is=2)
map(y, xlim = c(-10,30), ylim = c(50,70), cex=1, select=y1, cex.select=2, showall=TRUE)
#}
}
\keyword{select.station}
