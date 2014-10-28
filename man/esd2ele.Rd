\name{ele2param}
\alias{esd2ele}

\title{Dictionary and conversion tools between esd element identifier and variables names and specifications.}
\description{Converts between esd element/parameter identifier and names from
  different data sources.
}

\usage{
ele2param(ele = '101', src = 'GHCND')
esd2ele(param='t2m',src='GHCND',verbose=FALSE)
}

\arguments{
  \item{param,ele}{Parameter or element identifier. There are several core parameters or elements as well as a number of additional parameters. The parameters or elements are :
    PARAMETER                LONGNAME                                   ELE ID           
    auto                     Automatic selection.
    prcp, pr, rr, precip     Precipitation (mm)                         '601'
    tas, tavg, t2m, t2       2m-surface temperature (in degrees Celcius)'101'
    tmax, tasmax             Maximum temperature (in degrees Celcius)   '111'
    tmin, tasmin             Minimum temperature (in degrees Celcius)   '121'
    slp, mslp                Mean sea level pressure (hPa)              '401'
    cloud                    Cloud cover (\%)                           '801'
    ...
  }
  \item{src}{A character string for the acronym of the data source. The
    data sources are :
    NACD     : North Atlantic Climatological Dataset - Monthly
    NARP     : North Atlantic Research Programme ? - Monthly
    NORDKLIM : Nordic co-operation within climate activities - Monthly
    GHCNM    : Global historical climate network - Monthly
    ECAD     : European Climate Assessement & Dataset- Daily
    GHCND    : Global Historical Climate Network - Daily
    METNO    : MET Norway climate network - Daily
  }
}
\value{A meta data matrix object with the glossary of the different variables or
  element identifiers as originally defined by each data source}

\seealso{
  \code{\link{}}
}
\author{A. Mezghani, MET Norway}
\examples{
\dontrun{
# Eg.1 # Display the glossary of paramerters or element identifiers for 'GHCND' data source.
print(ele2param(ele=NULL,src='GHCND'))
# Eg.2 # Display the glossary of parameters or element identifiers for all data sources. 
print(ele2param())
# Eg.3 # Convert mean temperature parameter (param) to esd element (ele).
ele <- esd2ele(param='t2m')
print(ele)
}
}
\keyword{parameter,element}