\name{PCA.trajectory}
\alias{PCA.trajectory}
\alias{pca2trajectory}
\alias{plot.PCA.trajectory}

\title{Principle component analysis of trajectory objects.}
\description{
  Computes principal component analysis for trajectory data, e.g., storm tracks.
  Add some reference and details about the method.
  The PCA is based on \code{\link{svd}}.
}

\usage{
PCA.trajectory(X,...,neofs=20,param=c('lon','lat'),anomaly=TRUE,verbose=FALSE)
pca2trajectory(X)
}
\arguments{
  \item{X}{a 'trajectory' object}
  \item{verbose}{TRUE - clutter the screen with messages}
  \item{anomaly}{logical. If TRUE, subtract the first latitude/longitude
    from each trajectory.}
  \item{param}{parameters to include in principle component analysis.}
}

\author{K. Parding}
\examples{
# Simple EOF for annual mean SST:
data(imilast.M03)
x <- subset(imilast.M03,is=list(lon=c(-20,20),lat=c(50,70)))
# PCA of longitude and latitude
pca <- PCA.trajectory(x,param=c('lon','lat'))
plot.pca.trajectory(pca)
map.pca.trajectory(pca,projection='latlon')

# latitude only
pca <- PCA.trajectory(x,param=c('lat'))
plot.pca.trajectory(pca)
}

\keyword{spatial}
\keyword{ts}
\keyword{multivariate}

