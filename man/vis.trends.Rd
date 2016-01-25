\name{vis.trends}
\alias{vis.trends}

\title{Visualise trends for multiple overlapping periods}
\description{Produce a plot showing trends for multiple periods within
  a time series. The strength of the trend is represented by the color
  scale and significant trends are marked with black borders.}
\usage{
vis.trends(x, unitlabel = "unit", varlabel = "", pmax = 0.01, minlen =
15, lwd = NA, vmax = NA, new = TRUE, show.significance = TRUE, verbose =
FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{the 'x' argument provides the time series for which the trend
    analysis is performed. Only zoo objects are accepted.}
  \item{minlen}{minimum time interval to calculate trends for in units
    of years.}
  \item{unitlabel}{unit of x.}
  \item{varlabel}{name of x.}
  \item{vmax}{upper limit of trend scale.}
  \item{show.significance}{TRUE to mark statistically significant trends.}
  \item{pmax}{maximum p-value of trends marked as significant.}
}

\author{Kajsa Parding}

\examples{

t <- seq(as.Date("1955-01-01"),as.Date("2004-12-31"),by=1)
x <- zoo(sample(seq(-30,30,1e-1),length(t),rep=TRUE),order.by=t)
vis.trends(x,show.significance=FALSE)

data(Oslo)
vis.trends(Oslo, unitlabel="oC", varlabel = "Temperature",
  pmax = 1e-2, minlen = 40)
vis.trends(subset(Oslo,it='jja'), unitlabel="oC",
  varlabel = "Temperature JJA",
  pmax = 1e-3, vmax=0.5, minlen = 40)
vis.trends(subset(Oslo,it='mam'), unitlabel="oC",
  varlabel = "Temperature MAM",
  pmax = 1e-3, vmax=0.5, minlen = 40)
}

\keyword{ trend }
