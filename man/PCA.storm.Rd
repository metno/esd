\name{PCA.storm}
\alias{PCA.storm}
\alias{pca2storm}
\alias{plot.PCA.storm}

\title{Principle component analysis of storm tracks.}
\description{
  Computes principal component analysis for storm track data.
  Add some reference and details about the method.
  The PCA is based on \code{\link{svd}}.
}

\usage{
PCA.storm(X,neofs=20,param=c('lon','lat','slp'),anomaly=TRUE,verbose=FALSE)
pca2storm(X)
}
\arguments{
  \item{X}{a 'storm' object}
  \item{verbose}{TRUE - clutter the screen with messages}
  \item{anomaly}{logical. If TRUE, subtract the first latitude/longitude
    from each storm track.}
  \item{param}{parameters to include in principle component analysis.}
}

\author{K. Parding}
\examples{
# Simple EOF for annual mean SST:
data(imilast.M03)
x <- subset(imilast.M03,is=list(lon=c(-20,20),lat=c(50,70)))
pca <- PCA.storm(x)
plot.pca.storm(pca)
map.pca.storm(pca,projection='latlon')
}

\keyword{spatial}
\keyword{ts}
\keyword{multivariate}

