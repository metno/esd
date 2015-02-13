\name{test.retrieve.ncdf4}
\alias{test.retrieve.ncdf4}
\alias{test.cmip3.arctic}
\alias{test.cmip5.arctic}
\alias{test.cmip5.global}
\alias{test.cmip3.global}
\alias{test.cmip5.scandinavia}
\alias{test.cmip3.scandinavia}
\alias{test.cmip35.global}
\alias{}
\alias{}
\alias{}
\title{Test functions}
\description{This routine contains a series of test functions which compute the global mean 2m-temperature anomalies and other predefined regions based on both CMIP3 and CMIP5 experiments. The main function is \link{test.retrieve.ncdf4}. The others are specific cases for the "Scandinavian" and "Arctic" regions for both CMIP3 and CMIP5 projects}
\usage{
test.retrieve.ncdf4(path="CMIP3.monthly/20c3m-sresa1b",param="auto",cntr= "Global",
                    lon=NULL,lat=NULL,time=NULL,lev=NULL,experiment="CMIP3",
                    climatology=c(1986,2005),fig=TRUE,saveinfile=TRUE,verbose=FALSE)

test.cmip3.arctic(...)
test.cmip5.arctic(...)
test.cmip3.global(...)
test.cmip5.global(...)
test.cmip3.scandinavia(...)
test.cmip5.scandinavia

test.cmip35.global()
}

\arguments{
  \item{path}{Character vector of the path to the CMIP data.}
  \item{param}{Parameter or element type. There are several core parameters or elements as well as a number of additional parameters. The parameters or elements are #(Abdelkader : We need to update this list):
           auto         = automatic selection.
           prcp, pr     = Precipitation (mm)
   	   tas, tavg    = 2m-surface temperature (in degrees Celcius)
           tmax, tasmax = Maximum temperature (in degrees Celcius)
           tmin, tasmin = Minimum temperature (in degrees Celcius)
                }
  \item{lon}{Numeric value of longitude for the reference point (in decimal degrees East) or a vector containing the range of longitude values in the form of c(lon.min,lon.max)}
  \item{lat}{Numeric value of latitude for the reference point (in decimal degrees North) or a vector containing the range of latitude values in the form of c(lat.min,lat.max)}
  \item{lev}{Numeric value of pressure levels or a vector containing the range of pressure level values in the form of c(lev.min,lev.max) }
  \item{time}{Numerical year values or date values of time or a vector containing the range of values in the form of c(time.min,time.max). Date format should be in the form of "YYYY-MM-DD".}
  \item{saveinfile}{Logical or a character value. The output filename is set automatically if not specified. Default value is "TRUE".}
  \item{verbose}{Logical value. If TRUE, do not display comments (silent mode)}
}
\value{A field object}

\seealso{
}

\author{A.  Mezghanil}
\examples{
# Eg.1 : 
# Compute the global mean surface temperature anomalies from CMIP3 experiment
test.cmip3.global(...)
# Compute the global mean surface temperature anomalies from CMIP5 experiment
test.cmip5.global(...)
# Plot the Global t2m-temperature anomalies from both CMIP3 and CMIP5 experiments.
plot.cmip35.global()

# Eg.2 Compute the mean surface temperature anomalies from CMIP3 experiment over Arctic
test.cmip3.arctic(...) 

# Eg.3 Compute the mean surface temperature anomalies from CMIP5 experiment over Scandinavia.
test.cmip5.scandinavia(...)

}
\keyword{utilities}
