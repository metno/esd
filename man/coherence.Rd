\name{coherence}
\alias{coherence}
\title{Coherence spectrum - cross-spectrum analysis}
\description{
 Based on: \url{http://en.wikipedia.org/wiki/Wiener-Khinchin_theorem};
 Press et al. (1989) 'Numerical Recipes in Pascal', Cambridge, section 12.8
                     'Maximum Entropy (All Poles) Method';  
 von Storch & Zwiers (1999) 'Statistical Analysis in climate Research',
 Cambridge, section 11.4, eq 11.67, p. 235;
 
 A test with two identical series the original equation (eq 11.67) from von Storch
 & Zwiers (1999) gave uniform values: 1. The denominator was changed from
  \eqn{( \Gamma_{xx} * \Gamma_{yy} )}{(Gamxx * Gamyy )} to \eqn{(
    \sqrt{\Gamma_{xx} * \Gamma_{yy} )}}{sqrt(Gamxx * Gamyy )}.
 
}
\usage{
coherence(x,y,dt=1,M=NULL,plot=TRUE)
}
\arguments{
    \item{x}{A vector (time series).}
    \item{y}{A vector (time series).}
    \item{dt}{time incremet - for plotting.}
    \item{M}{Window length - default= half series length}
    \item{plot}{Flag: plot the diagnostics.}
 }
\value{A complex vector .}
\author{R.E. Benestad}
\examples{
\dontrun{
# Add example of usage
}
}
\keyword{manip}
