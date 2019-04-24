\name{crossval}
\alias{crossval}
\alias{crossval.ds}
\alias{crossval.list}
\alias{crossval.dsensemble}
\title{Cross-validation}
\description{
Applies a cross-validation of DS results, using the same strategy as in
the DS exercise. Any step-wise screening is applied for each iteration
independently of that used to identify the subset of skillful predictors
in the original analysis. The model coeffiecients (beta) is saved for
each iteration, and both correlation and root-mean-squared-error are
returned as scores.

\code{crossval.dsensemble} will make use of the \code{evaluation} attribute with cross-validation results and returns the correlation.
}
\usage{
crossval(x, m=5, ...)
crossval.ds(x, m=5, ...)
crossval.list(x, m=5, ...)
crossval.dsensemble(x,plot=TRUE)
}
\arguments{
    \item{x}{The results from \code{\link{DS}}.}
    \item{m}{window with - leave m-out for each iteration. There are
      also some pre-set options: 'cordex-esd-exp1', 'value-exp1', and 'loo' for
      experiments defined at CORDEX-ESD, COST-VALUE, and leave-one-out
      ('loo') cross-validation.}
 }
\value{Cross-validation object.}
\author{R.E. Benestad}
\examples{
data(Oslo)
t2m <- t2m.DNMI(lon=c(-20,40),lat=c(45,65))
eof <- EOF(t2m)

ds <- DS(Oslo,eof)
xv <- crossval(ds)
plot(xv)
}
\keyword{manip}
