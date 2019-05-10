\name{read.imilast}
\alias{read.imilast}
\alias{read.hurdat2}

\title{
Read cyclone data
}
\description{
Methods for reading cyclone data. 

\code{read.imilast} reads data from files following the standards of the IMILAST project (Neu et al., 2013: IMILAST: A Community Effort to Intercompare Extratropical Cyclone Detection and Tracking Algorithms. Bull. Amer. Meteor. Soc., 94, 529–547, https://doi.org/10.1175/BAMS-D-11-00154.1).

\code{read.hurdat2} reads data from the Atlantic hurricane database (http://www.nhc.noaa.gov/data/#hurdat).
}

\usage{
read.imilast(fname, path = NULL, verbose = FALSE)
read.hurdat2(fname='http://www.aoml.noaa.gov/hrd/hurdat/hurdat2-1851-2014-022315.html',
             path = NULL, verbose = FALSE)
}
\arguments{
  \item{fname}{filename (for \code{read.hurdat}, filename can also be a url)}
  \item{path}{path to file}
  \item{verbose}{Logical value defaulting to FALSE. If FALSE, do not
  display comments (silent mode). If TRUE, displays extra information on progress.}
}
\value{
 An "events" "data.frame" object containing the date, time, lon, and lat, as well as additional information (e.g., trajectory number, slp or other measure of storm strengt) of the cyclones. 
}
\author{K. Parding }

\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{cyclones}% use one of  RShowDoc("KEYWORDS")
\keyword{storms}% __ONLY ONE__ keyword per line
\keyword{IMILAST}
\keyword{hurdat2}
