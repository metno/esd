\name{CCA}
\alias{Canonical correlation analysis}
\alias{CCA}
\alias{CCA.default}
\alias{CCA.eof}
\alias{CCA.pca}
\alias{CCA.field}
\alias{test.cca}
\alias{predict.cca}
\title{Canonical correlation analysis}
\description{
Applies a canonical correlation analysis (CCA) to two data sets. The CCA here
can be carried out based on an \code{\link{svd}} based approach (after 
Bretherton et al. (1992), J. Clim. Vol 5, p. 541, also documented in Benestad (1998): 
"Evaluation of Seasonal Forecast Potential for Norwegian  Land Temperatures and Precipitation using CCA", 
DNMI KLIMA Report 23/98 at \url{http://met.no/english/r_and_d_activities/publications/1998.html}) or ii) a 
covariance-eigenvalue approach (after Wilks, 1995, "Statistical methods in the Atmospheric 
Sciences", Academic Press, p. 401). 

The analysis can also be applied to either EOFs or fields.

Note: the analysis has sometimes been somewhat unstable, returning inconsistent results. The recommendation is to use EOFs and SVD option.

The CCA analysis can be used to develope statistical models according to:

\deqn{Y = \Psi X}{Y = Psi X}

Where Y is the predictand and X the predictor.
\code{plotCCA} plots the CCA results, \code{testCCA} is for code
verification, and \code{Psi} returns the matrix \deqn{\Psi}{Psi}.

\code{stations2field} turns a group of station objects into a field by
the means of a simple and crude interpolation/gridding. \code{check.repeat} is a quality-control
function that eliminates repeated years in the station objects.

Try the same type of argument as in lm (' y ~ x, data= ')
}
\usage{
CCA(Y,X,...)
CCA.default(Y,X,...) 
CCA.eof(Y,X,...,ip=1:8,verbose=FALSE)
CCA.pca(Y,X,...,ip=1:8,verbose=FALSE)
CCA.field(Y,X,...,ip=1:8,verbose=FALSE)
test.cca(method="CCA",reconstr=FALSE,mode=1,test=TRUE,LINPACK=TRUE,
         SVD=TRUE,n.pc=4,synthetic=TRUE)

predict.cca <- function(object, newdata=NULL, ...)
}
\arguments{
    \item{Y}{An object with climate data: field, eof, pca.}
    \item{X}{Same as Y.}
    \item{SVD}{Use a singular value decomposition as a basis for the PCA.}
    \item{ip}{Which EOFs to include in the CCA.}
    \item{LINPACK}{an option for \code{\link{svd}}.}
    \item{object}{The result from CCA.}
    \item{newdata}{The same as X.}
 }
\value{A CCA object: a list containing a.m, b.m, u.k, v.k, and r, describing the Canonical Correlation variates, patterns and correlations. a.m and b.m are the patterns and u.k and v.k the vectors (time evolution).}
\author{R.E. Benestad}
\examples{
# CCA with two eofs
slp <- slp.NCEP(lat=c(-40,40),anomaly=TRUE)
sst <- sst.NCEP(lat=c(-40,40),anomaly=TRUE)
eof.1 <- EOF(slp, it='Jan')
eof.2 <- EOF(sst, it='Jan')
cca <- CCA(eof.1,eof.2)
plot(cca)

# CCA with PCA and EOF:
\dontrun{
NACD <- station.nacd()
plot(annual(NACD))
map(NACD,FUN="sd")
pca <- PCA(NACD)
plot(pca)
naslp <- slp.NCEP(lon=c(-30,40),lat=c(30,70),anomaly=TRUE)
map(naslp)
eof <- EOF(naslp,it='Jan')
nacca <- CCA(pca,eof)
plot(nacca)
cca.pre <- precit.cca(nacca)
}
}
\keyword{manip}
