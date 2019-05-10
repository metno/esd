\name{station.meta}
\alias{station.meta}
\docType{data}
\title{Weather station metadata}

\description{
Meta datasets based on different data sources or datasets included in
esd package. Mainly weather stations' coordinates and other meta data. 
}

\usage{data(station.meta)}

\format{
  The format is: \cr
List of 12 \cr
 $ station_id: chr [1:284239] "6447" "6193" "21100" "25140" ... \cr
 $ location  : chr [1:284239] "UCCLE" "HAMMERODDE_FYR" "VESTERVIG"
 "NORDBY" ... \cr
 $ country   : chr [1:284239] "BELGIA" "DENMARK" "DENMARK" "DENMARK"
 ... \cr
 $ longitude : num [1:284239] 4.35 14.78 8.32 8.4 10.6 ... \cr
 $ latitude  : num [1:284239] 50.8 55.3 56.8 55.4 55.9 ... \cr
 $ altitude  : num [1:284239] 100 11 18 5 11 9 4 51 85 105 ... \cr
 $ element   : chr [1:284239] "101" "101" "101" "101" ... \cr
 $ start     : chr [1:284239] "1833" "1853" "1873" "1872" ... \cr
 $ end       : chr [1:284239] NA NA NA NA ... \cr
 $ source    : chr [1:284239] "NACD" "NACD" "NACD" "NACD" ... \cr
 $ wmo       : num [1:284239] 6447 6193 -999 -999 -999 ... \cr
 $ quality   : int [1:284239] 2 2 2 2 2 1 1 2 5 5 ... \cr
 - attr(*, "history")= chr [1:7] \cr
 "meta2esd.R - data taken from the clim.pact package and consolidated for NACD and NARP" "nordklim.meta.rda" "ecad.meta.rda" "ghcnd.meta.rda" ... \cr
 - attr(*, "date")=function ()  \cr
 - attr(*, "call")= language meta2esd(save = TRUE) \cr
 - attr(*, "author")= chr "R.E. Benestad & A. Mezghani" \cr
 - attr(*, "URLs")= chr [1:6] "www.dmi.dk/dmi/sr96-1.pdf" \cr
 "http://www.norden.org/en/publications/publikationer/2005-450" \cr
 "http://www.smhi.se/hfa_coord/nordklim/" "http://eca.knmi.nl/" ... \cr
}
\examples{
data(station.meta)
str(station.meta)
map(station.meta)
}

\seealso{\code{\link{select.station}},\code{\link{station}}}
\keyword{datasets}
