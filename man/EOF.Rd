\name{EOF}
\alias{EOF}
\alias{EOF.default}
\alias{EOF.field}
\alias{EOF.comb}
\alias{eof2field}
\alias{PCA}
\alias{PCA.default}
\alias{PCA.station}
\alias{pca2station}
\alias{Empirical orthogonal Functions}

\title{Empirical Orthogonal Functions (EOFs).}
\description{
  Computes EOFs (a type of principal component
  analysis) for combinations of data sets, typically from
  gridded data, reanalysis and climate models results.
  
  [ref: Benestad (2001), "A comparison between two empirical downscaling
  strategies", \emph{Int. J. Climatology}, \bold{vol 21}, Issue 13,
  pp.1645-1668. DOI 10.1002/joc.703]. and  \code{\link{mixFields}} prepares for
  mixed-field EOF analysis [ref. Bretherton et al. (1992) "An
  Intercomparison of Methods for finding Coupled Patterns in Climate
  Data", \emph{J. Climate}, \bold{vol 5}, 541-560;
  Benestad et al. (2002), "Empirically downscaled temperature scenarios
  for Svalbard", \emph{Atm. Sci. Lett.}, doi.10.1006/asle.2002.0051]. 

  Uncertainty estimates are computed according to North et al. (1982),
  "Sampling Errors in the Estimation of Empirical Orthogonal Functions",
  \emph{Mon. Weather Rev.}, \bold{vol 110}, 699-706. 

  The EOFs are based on \code{\link{svd}}.

  See the course notes from Environmental statistics for climate
  researchers \url{http://www.gfi.uib.no/~nilsg/kurs/notes/course.html}
  for a discussion on EOF analysis.

  The method \code{PCA} is similar to EOF, but designed for parallel
  station series (e.g. grouped together with \code{\link{merge}}). PCA
  does not assume gridded values and hence does not weigth according to
  grid area. PCA is useful for downscaling where the spatial
  covariance/coherence is important, e.u involving different variables
  from same site, same variable from different sites, or a mix between
  these. For instance, PCA can be applied to the two wind components
  from a specific site and hence extract the most important wind
  directions/speeds. 
}

\usage{
EOF(X,it=NULL,n=20,lon=NULL,lat=NULL,verbose=FALSE,...)
EOF.default(X,it=NULL,n=20,lon=NULL,lat=NULL,
            area.mean.expl=TRUE,verbose=FALSE,...)
EOF.field(X,it=NULL,n=20,lon=NULL,lat=NULL,
          area.mean.expl=TRUE,verbose=FALSE)
EOF.comb(X,it=NULL,n=20,lon=NULL,lat=NULL,area.mean.expl=TRUE,verbose=FALSE)
eof2field(x,it=NULL,is=NULL,anomaly=FALSE,verbose=FALSE)
PCA(X,...) 
EOF.default(X,...)
PCA.station(X,na.action='fill',verbose=FALSE)
pca2station(X,lon=NULL,lat=NULL,anomaly=FALSE,what='pca',verbose=FALSE)
}
\arguments{
  \item{X}{a 'field' or 'pca' object}
  \item{it}{see \code{\link{subset}}}
  \item{n}{number of EOFs}
  \item{is}{Spatial subsetting - see \code{\link{subset.eof}}}
  \item{lon}{set longitude range - see \code{\link{t2m.ERAINT}}}
  \item{lat}{set latitude range}
  \item{verbose}{TRUE - clutter the screen with messages}
  \item{area.mean.expl}{When TRUE, subtract the area mean for the domain
  and use as a the first co-variate before the PCs from the EOF
  analysis.}
  \item{na.action}{'fill' uses \code{\link{approx}} to interpolate the
    NA-values before the PCA.}
  \item{what}{Default set to 'pca' for convertin the PCA-results, but
    also works on downscaled PCA results. Option 'xval' returns the
    station values for cross-validation results and 'test' returns the
    station values used in the cross-validation (same interval as for 'xval').}
}
\value{File containing an 'eof' object which is based on the 'zoo' class.

}
\author{R.E. Benestad}
\examples{
# Simple EOF for annual mean SST:
sst <- sst.NCEP(lon=c(-90,20),lat=c(0,70))
SST <- aggregate(sst, year, mean, na.rm = FALSE)
eof.sst <- EOF(SST)
plot(eof.sst)

# EOF of July SST:
eof.sst7 <- EOF(sst,it=7)
plot(eof.sst7)

# common EOF for model
# Get some sample data, extract regions:
GCM <- t2m.NorESM.M()
gcm <- subset(GCM,is=list(lon=c(-50,60),lat=c(30,70)))
t2m.eraint <- t2m.ERAINT()
eraint <- subset(t2m.eraint,is=list(lon=c(-50,60),lat=c(30,70)))

OBS <-  aggregate(eraint, by=year, mean, na.rm = FALSE)
GCM <- aggregate(gcm, by=year, mean, na.rm = FALSE)
OBSGCM <- combine(OBS,GCM,dimension='time')

ceof <- EOF(OBSGCM)
plot(ceof)

# Example for using PCA in downscaling
## nacd <- station(src='nacd')
## X <- annual(nacd)
X <- station(src='nacd')
nv <- function(x) sum(is.finite(x))
ok <- (1:dim(X)[2])[apply(X,2,nv) == dim(X)[1]]
X <- subset(X,is=ok)
pca <- PCA(X)
map(pca)

slp <- slp.NCEP(lon=c(-20,30),lat=c(30,70))
eof <- EOF(slp,it=1)
ds <- DS(pca,eof)
# ds is a PCA-object
plot(ds)

# Recover the station data:
Z <- pca2station(pca)
plot(Z,plot.type='multiple')


}
\keyword{spatial}
\keyword{ts}
\keyword{multivariate}

