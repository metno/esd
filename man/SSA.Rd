\name{SSA}
\alias{Singular Spectrum Analysis}
\alias{SSA}
\alias{plotSSA}
\title{Singular Spectrum Analysis}
\description{After von Storch & Zwiers (1999), Statistical Analysis in Climate Research, p. 312}
\usage{
SSA(x,m,plot=TRUE,main="SSA analysis",sub="",LINPACK=TRUE,
    param = "t2m", anom = TRUE,i.eof=1)
plotSSA(ssa,main="SSA analysis",sub="")
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
    \item{i.eof}{If x is an eof-object, which PC to use.}
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
