\name{Oslo}
\alias{Oslo}
\docType{data}
\title{
%%   ~~ data name/kind ... ~~
Oslo monthly mean temperature time series
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
Oslo temperature monthy record from 1837 up to now. 
}
\usage{data(Oslo)}
\format{
  The format is:
  'zoo' series from 1837-01-01 to 2014-02-01 \cr 
  Data: atomic [1:2126] NaN NaN NaN 3.1 8.9 13.5 16 15.7 11.2 7.6 ... \cr 
  - attr(*, "location")= chr "Oslo" \cr
  - attr(*, "variable")= chr "T2m" \cr
  - attr(*, "unit")= chr "deg C" \cr
  - attr(*, "longitude")= num 10.7 \cr
  - attr(*, "latitude")= num 59.9 \cr
  - attr(*, "altitude")= num 94 \cr
  - attr(*, "country")= chr "Norway" \cr
  - attr(*, "longname")= chr "temperature at 2m" \cr
  - attr(*, "station_id")= num 18700 \cr
  - attr(*, "quality")= chr "homogenised" \cr
  - attr(*, "calendar")= chr "gregorian" \cr
  - attr(*, "source")= chr "Dr. Nordli, 2013, met.no" \cr
  - attr(*, "URL")= logi NA \cr
  - attr(*, "type")= chr "observation"\cr 
  - attr(*, "aspect")= chr "original" \cr
  - attr(*, "reference")= chr  "Nordli et  al. (in progress). 'The Oslo Temperature series 1837-2012: Homogeneity testing and Climate Analysis'" \cr
  - attr(*, "info")= logi NA \cr
  - attr(*, "method")= chr "Blended recnostruction (1877-1936) and instrumental data (1937-)" \cr
  - attr(*, "history")=List of 3 \cr
  ..$ call       :length 19 as.station.data.frame(oslo, loc = "Oslo",
  param = "T2m", unit = "deg C",      lon = 10.7, lat = 59.9, alt = 94,
  cntr = "Norway", longname = "temperature at 2m",  ... \cr
  .. ..- attr(*, "srcref")=Class 'srcref'  atomic [1:8] 85 1 90 104 1
  104 85 90 \cr
  .. .. .. ..- attr(*, "srcfile")=Classes 'srcfilecopy', 'srcfile'
  <environment: 0x23b1980>  \cr
  ..$ timestamp  : chr "Tue Dec 10 14:53:51 2013" \cr
  ..$ sessioninfo:List of 3 \cr
  .. ..$ R.version  : chr "R version 3.0.2 (2013-09-25)" \cr
  .. ..$ esd.version: chr "esd_0.2-1" \cr
  .. ..$ platform   : chr "x86_64-pc-linux-gnu (64-bit)" \cr
  Index:  Date[1:2126], format: "1837-01-01" "1837-02-01" "1837-03-01"
  "1837-04-01" ... \cr
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
Oslo surface temperature recorded on a monthly basis from 1837 up to
2012. It corresponds to a blended recnostruction (1837-1936) and
instrumental data (1937-2012). An homogenisation
procedure has been carried out on the data by Øyvind Nordli at MET Norway. 
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
MET Norway
}
\references{
%%  ~~ possibly secondary sources and usages ~~
Nordli et  al. (in progress). 'The Oslo Temperature series 1837-2012: Homogeneity testing and Climate Analysis'
}
\examples{
data(Oslo)
## maybe str(Oslo) ; plot(Oslo) ...
}
\keyword{datasets}
