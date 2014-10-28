\name{diagnose}
\alias{diagnose}
\alias{diagnose.default}
\alias{diagnose.comb}
\alias{diagnose.eof}
\alias{diagnose.comb.eof}
\alias{diagnose.mvr}
\alias{diagnose.cca}
\alias{diagnose.ds}
\alias{diagnose.station}
\alias{diagnose.dsensemble}
\title{Diagnose}
\description{Diagnose and examine combined fields, MVR, and CCA
  results. applies some tests to check for consistency.

  The method \code{diagnose.comb.eof} which estimates the difference in
  the mean for the PCs of the calibration data and GCMs over a common
  period in addition to the ratio of standard deviations and lag-one
  autocorrelation. This 'bias correction' is described in Imbert and
  Benestad (2005), \emph{Theor. Appl. Clim.}
  \url{http://dx.doi.org/10.1007/s00704-005-0133-4}.

  \code{climvar} estimates the climatological variance, e.g. how the
  inter-annual variance varies with seasons.
}
\usage{
diagnose(x,...)
diagnose.default(x,...)
diagnose.comb(x)
diagnose.eof(x)
diagnose.comb.eof(x)
diagnose.mvr(x)
diagnose.cca(x)
diagnose.ds(x,plot=FALSE)
diagnose.station(x,it=NULL,...)
diagnose.dsensemble(x,plot=TRUE,plot.type='target',...)
}
\arguments{
  \item{x}{data object}
  \item{it}{teporal selection - see \code{\link{subset}}}
  \item{plot}{}
  \item{plot.type}{}
 }
\value{A 'diag' object containing test results}

\seealso{
}

\author{R.E. Benestad}
\examples{
t2m <- t2m.NCEP(lon=c(-40,40),lat=c(30,70))
T2m <- t2m.NorESM.M(lon=c(-40,40),lat=c(30,70))
# Combine in time to compute common EOFs:
X <- combine(t2m,T2m)
diagnose(X)

ceof <-EOF(X,it=1)
plot(diagnose(ceof))

slp <- slp.NCEP(lat=c(-40,40),anomaly=TRUE)
sst <- sst.NCEP(lat=c(-40,40),anomaly=TRUE)
eof.1 <- EOF(slp,it=1)
eof.2 <- EOF(sst,it=1)
cca <- CCA(eof.1,eof.2)
diagnose(cca)


}
\keyword{utilities}
