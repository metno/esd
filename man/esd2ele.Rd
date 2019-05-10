\name{ele2param}
\alias{ele2param}
\alias{esd2ele}

\title{Dictionary and conversion tools between esd element identifier and variables names and specifications.}

\description{Converts between esd element/parameter identifier and
  variable names from different data sources.
}

\usage{
ele2param(ele = '101', src = 'GHCND')
esd2ele(param='t2m',src='GHCND',verbose=FALSE)
}

\arguments{
  \item{param,ele}{Parameter or element identifier. There are several
  core parameters or elements as well as a number of additional
  parameters. The parameters or elements are :}

\tabular{lll}{
  PARAMETER   \tab LONGNAME   \tab    ELE ID     \cr      
  auto         \tab  Automatic selection.  \tab  \cr
  prcp, pr, rr, precip   \tab   Precipitation (mm)      \tab   '601' \cr
  tas, tavg, t2m, t2   \tab  2m-surface temperature (in degrees
  Celcius)  \tab '101' \cr
  tmax, tasmax \tab           Maximum temperature (in degrees Celcius)
  \tab  '111' \cr
  tmin, tasmin \tab            Minimum temperature (in degrees
  Celcius) \tab  '121' \cr
  slp, mslp    \tab            Mean sea level pressure (hPa)
  \tab  '401' \cr
  cloud        \tab            Cloud cover (\%)
  \tab  '801' \cr
}

\item{src}{A character string for the acronym of the data source. The
  data sources are :}
\tabular{lll}{
  Data Source  \tab Full name   \tab Type \cr 
  NACD     \tab North Atlantic Climatological Dataset \tab  Monthly \cr
  NARP     \tab North Atlantic Research Programme ? \tab Monthly \cr
  NORDKLIM \tab Nordic co-operation within climate activities \tab
  Monthly \cr
  METNOM   \tab MET Norway climate network \tab Monthly \cr
  GHCNM    \tab Global historical climate network \tab Monthly \cr
  ECAD     \tab European Climate Assessement & Dataset \tab Daily \cr
  GHCND    \tab Global Historical Climate Network \tab Daily \cr
  METNOD   \tab MET Norway climate network \tab Daily \cr
}
}

\value{A meta data matrix object with the glossary of the different variables or
  element identifiers as originally defined by each data source}

\author{A. Mezghani, MET Norway}
\examples{
# Eg.1 # Display the glossary of paramerters or element identifiers for 'GHCND' data source.
print(ele2param(ele=NULL,src='GHCND'))
# Eg.2 # Display the glossary of parameters or element identifiers for all data sources. 
print(ele2param())
# Eg.3 # Convert mean temperature parameter (param) to esd element (ele).
ele <- esd2ele(param='t2m')
print(ele)

}
\keyword{parameter,element}