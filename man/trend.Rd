\name{trend}
\alias{trend.station}
\alias{trend.eof}
\alias{trend.field}
\title{Trending and detrending data}
\description{
  Trend analysis and de-trending of data.ls 
}
\usage{
trend(x,result="trend",model="y ~ t",...)
trend.default(x,result="trend",model="y ~ t",...)
trend.one.station(x,result="trend",model="y ~ t",...)
trend.station(x,result="trend",model="y ~ t",...)
trend.eof(x,result="trend",model="y ~ t",...)
trend.field(x,result="trend",model="y ~ t",...)
trend.zoo(x,result="trend",model="y ~ t",...) 
}
\arguments{
  \item{x}{The data object}
  \item{result}{"trend" returns the trend; "residual" returns the residual}
  \item{model}{The trend model used by \code{\link{lm}}.}
}
\value{Similar type object as the input object}

\seealso{\code{link{climatology}}, \code{link{anomaly}}
}

\author{R.E. Benestad}
\examples{
data(ferder)
tr <- trend(annual(ferder,'max'))
attr(tr,'coefficients')
}
\keyword{utilities}
