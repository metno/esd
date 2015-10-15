\name{crossval}
\alias{crossval}
\alias{crossval.ds}
\alias{crossval.list}
\title{Cross-validation}
\description{
Applies a cross-validation of DS results, using the same strategy as in
the DS exercise. Any step-wise screening is applied for each iteration
independently of that used to identify the subset of skillful predictors
in the original analysis. The model coeffiecients (beta) is saved for
each iteration, and both correlation and root-mean-squared-error are
returned as scores.
}
\usage{
crossval(x, m=5, ...)
crossval.ds(x, m=5, ...)
crossval.list(x, m=5, ...)
}
\arguments{
    \item{x}{The results from \code{\link{DS}}.}
    \item{m}{window with - leave m-out for each iteration. There are
      also some pre-set options: 'cordex-esd-exp1' and 'value-exp1' for
      experiments defined at CORDEX-ESD and COST-VALUE.}
 }
\value{Cross-validation object.}
\author{R.E. Benestad}
\examples{
data(Oslo)
t2m <- t2m.NCEP(lon=c(-20,40),lat=c(45,65))
eof <- EOF(t2m,1)

ds <- DS(Oslo,eof)
xv <- crossval(ds)
plot(xv)
}
\keyword{manip}
