\name{trend}
\alias{trend.one.station}
\alias{trend.station}
\alias{trend.eof}
\alias{trend.field}
\alias{trend.zoo}
\alias{trend.zoo.multi}
\alias{trend.coef}
\alias{trend.err}
\alias{trend.pval}
\title{Trending and detrending data}
\description{
  Trend analysis and de-trending of data. The three methods \code{trend.coef},
  \code{trend.err} and \code{trend.pval} are somewhat different to the
  other trend methods and designed for the use in \code{apply}
  operations, as reflected in the different sets of arguments. They are
  used in the other methods if the \code{result} argument is set to one
  of ["coef","err","pval"]. 
}
\usage{
trend(x,result="trend",model="y ~ t",verbose=FALSE,...)
trend.default(x,result="trend",model="y ~ t",verbose=FALSE,...)
trend.one.station(x,result="trend",model="y ~ t",verbose=FALSE,...)
trend.station(x,result="trend",model="y ~ t",verbose=FALSE,...)
trend.eof(x,result="trend",model="y ~ t",verbose=FALSE,...)
trend.field(x,result="trend",model="y ~ t",verbose=FALSE,...)
trend.zoo(x,result="trend",model="y ~ t",verbose=FALSE,...)
trend.zoo.multi(x,result="trend",model="y ~ t",verbose=FALSE,...)
trend.coef(x,...)
trend.err(x,...)
trend.pval(x,...)
}
\arguments{
  \item{x}{The data object}
  \item{result}{"trend" returns the trend; "residual" returns the
    residual; "coef" returns the trend coefficient; "err" the error
    estimate; "pval" the p-value.}
  \item{model}{The trend model used by \code{\link{lm}}.}
}
\value{Similar type object as the input object}

\seealso{\code{link{climatology}}, \code{link{anomaly}}
}

\author{R.E. Benestad}
\examples{
data(ferder)

plot(annual(ferder,'max'))
tr <- trend(annual(ferder,'max'))
lines(tr)
grid()
print(attr(tr,'coefficients'))
print(trend(ferder,results='pval'))
}
\keyword{utilities}
