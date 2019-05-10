\name{SSA}
\alias{Singular Spectrum Analysis}
\alias{SSA}
\alias{plot.ssa}
\title{Singular Spectrum Analysis}
\description{After von Storch & Zwiers (1999), Statistical Analysis in Climate Research, p. 312}
\usage{
SSA(x, m, plot=TRUE, main="SSA analysis", sub="",
    LINPACK=TRUE,param="t2m", anom=TRUE, ip=1,
    verbose=FALSE)
plotSSA(ssa, main="SSA analysis",sub="")
}
\arguments{
    \item{x}{A station or eof object.}
    \item{m}{Window length.}
    \item{plot}{Flag: plot the diagnostics.}
    \item{LINPACK}{'TRUE': svd; 'FALSE':La.svd}
    \item{main}{main title (see \code{link{plot}}).}
    \item{sub}{subtitle (see \code{link{plot}}).}
    \item{ssa}{An 'SSA' object returned by SSA().}
    \item{param}{Which parameter ("daily.station.record") to use: "precip", "t2m" or other.}
    \item{anom}{TRUE if analysis on anomalies}
    \item{ip}{If x is an eof-object, which PC to use.}
    \item{verbose}{Print out diagnostics.}
 }
\value{A SSA object: An \code{link{svd}} object with additional parameters: m (window length), nt (original length of series), Nm (effective length of series= nt - m), anom (FLAG for use of anomaly), param (name of parameter, typically 'precip' or 't2m'), station (the station object to which SSA is applied).
}
\author{R.E. Benestad}
\examples{
\dontrun{
data(DNMI.t2m)
eof.1 <- EOF(DNMI.t2m,mon=1)
pop <- POP(eof.1)
}
}
\keyword{manip}
