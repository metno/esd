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
\alias{diagnose.distr}
\alias{diagnose.matrix}
\alias{diagnose.dsensemble}
\alias{diagnose.ds.pca}
\title{Diagnose}
\description{Diagnose and examine combined fields, MVR, and CCA
  results. applies some tests to check for consistency.

  The method \code{diagnose.comb.eof} which estimates the difference in
  the mean for the PCs of the calibration data and GCMs over a common
  period in addition to the ratio of standard deviations and lag-one
  autocorrelation. The x-axis shows the difference in the mean of the
  segments in the PCs representing the different data souces, the y-axis
  shows difference in standard deviation and the size of the symbols the
  difference in the autocorrelation (open symbols if the autocorrelation
  have different signs).

  \code{climvar} estimates the climatological variance, e.g. how the
  inter-annual variance varies with seasons.

 
}
\usage{
diagnose(x,...)
diagnose.default(x,...)
diagnose.comb(x,...)
diagnose.eof(x,...)
diagnose.comb.eof(x,...)
diagnose.mvr(x,...)
diagnose.cca(x,...)
diagnose.ds(x,...,plot=FALSE,verbose=FALSE,new=TRUE)
diagnose.ds.pca(x,...,plot=FALSE,verbose=FALSE,new=TRUE)
diagnose.station(x,...,it=NULL,...)
diagnose.matrix(x,...,xlim=NULL,ylim=NULL,verbose=FALSE)
diagnose.dsensemble(x,...,plot=TRUE,map=TRUE,plot.type='target')
diagnose.dsensemble.list(x,...,plot=FALSE,is=NULL,ip=NULL,
                 map.show=TRUE,alpha=0.6,xrange=NULL,yrange=NULL,
                 main=NULL,verbose=FALSE,new=TRUE)
diagnose.distr(x,...,main='Mean - quantile',
               xlab='mean',ylab=expression(q[p]),
               sub=src(x),probs=0.95)
}
\arguments{
  \item{x}{data object}
  \item{it}{teporal selection - see \code{\link{subset}}}
  \item{plot}{if TRUE, plot results}
  \item{plot.type}{type of plot}
  \item{verbose}{Logical value defaulting to FALSE. If FALSE, do not
  display comments (silent mode). If TRUE, displays extra information on progress.}
  \item{new}{if TRUE plot in new window}
  \item{xlim}{range of x-axis}
  \item{xlim}{range of y-axis}
  \item{alpha}{factor modifying the opacity alpha; typically in [0,1]}
  \item{map.show}{if TRUE show map}
  \item{xrange}{longitude range to display in map}
  \item{yrange}{latitude range to display in map}
  \item{main}{main label in plot}
  \item{sub}{smaller label (subtitle) in plot}
  \item{xlab}{label of x-axis}
  \item{ylab}{label of y-axis}
  \item{probs}{quantile to display in plot, e.g., probs=0.95 
               gives a diagnosis of the 95th percentile of the data.}
 }
\value{A 'diag' object containing test results}

\author{R.E. Benestad}
\examples{
t2m <- t2m.NCEP(lon=c(-40,40),lat=c(30,70))
T2m <- t2m.NorESM.M(lon=c(-40,40),lat=c(30,70))
# Combine in time to compute common EOFs:
X <- combine(t2m,T2m)
diagnose(X)

ceof <- EOF(X,it="jan")
plot(diagnose(ceof))

slp <- slp.NCEP(lat=c(-40,40),anomaly=TRUE)
sst <- sst.NCEP(lat=c(-40,40),anomaly=TRUE)
eof.1 <- EOF(slp,it="jan")
eof.2 <- EOF(sst,it="jan")
cca <- CCA(eof.1,eof.2)
diagnose(cca)

}
\keyword{utilities}
