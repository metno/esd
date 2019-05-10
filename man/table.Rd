\name{summary.dsensemble}
\alias{summary.dsensemble}
\alias{summary.station}
\alias{summary.ds}
\alias{summary.eof}
\alias{summary.cca}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Summary showing summary of objects
%%  ~~function to do ... ~~
}
\description{
Produce a summary table
}
\usage{
summary.dsensemble(x, years = seq(1990, 2090, by = 20))
summary.station(x,im=1:12)
summary.ds(x)
summary.eof(x)
summary.cca(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{an object of type 'DSensemble', 'station', 'DS', 'EOF' or 'CCA'}
  \item{years}{A set of years for which to produce summary statistics}
  \item{im}{The order of months in the table. Use \code{im=c(10:12,1:9)} for Oct-Sep.}
}
\value{A matrix containing summary statistics}

\examples{
data("Oslo")
summary(Oslo)
}
