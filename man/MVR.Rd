\name{MVR}
\alias{Multi-variate regression}
\alias{MVR}
\alias{predict.MVR}
\alias{plot.MVR}
\title{Multi-variate regression}
\description{MVR solves the equation
\deqn{Y = \Psi X}{Y = Psi X}
and estimates \deqn{\Psi}{Psi} by inverting the equation. Predictions
give the varlue of Y, given this matrix and some input. MVR is useful
for data where Y contains several time series where the spatial
coherence/covariance is important to reproduce. For instance, Y may be a
combination of stations, the two wind components from one station, or a
set of different elements from a group of stations. 

}
\usage{
MVR(Y,X,...) 
MVR.default(Y,X,...) 
MVR.field(Y,X,SVD=TRUE,LINPACK=FALSE) 
MVR.pca(Y,X,SVD=TRUE,LINPACK=FALSE) 
MVR.eof(Y,X,SVD=TRUE,LINPACK=FALSE) 
predict.MVR(object, newdata=NULL, ...)
}
\arguments{
    \item{Y}{An object with climate data: field, eof, or pca.}
    \item{X}{Same as Y or any zoo object.}
    \item{SVD}{Use a singular value decomposition as a basis for the PCA.}
    \item{LINPACK}{an option for \code{\link{svd}}.}
    \item{object}{The result from CCA.}
    \item{newdata}{The same as X.}
 }
\value{A CCA object: a list containing a.m, b.m, u.k, v.k, and r, describing the Canonical Correlation variates, patterns and correlations. a.m and b.m are the patterns and u.k and v.k the vectors (time evolution).}
\author{R.E. Benestad}
\examples{
\dontrun{
# Example for using EOF and MVR
slp <- slp.NCEP(lat=c(-40,40),anomaly=TRUE)
sst <- sst.NCEP(lat=c(-40,40),anomaly=TRUE)
eof.1 <- EOF(slp,mon=1)
eof.2 <- EOF(sst,mon=1)
mvr <- MVR(eof.1,eof.2)
plot(mvr)

# Example for using PCA and MVR
oslo <- station(src="NACD",loc="Oslo")
bergen <- station.nacd("Bergen")
stockholm <- station.nacd("Stockholm")
copenhagen <- station.nacd("Koebenhavn")
helsinki <- station.nacd("Helsinki")
reykjavik <- station.nacd("Stykkisholmur")
edinburgh <- station.nacd("Edinburgh")
debilt <- station.nacd("De_Bilt")
uccle <- station.nacd("Uccle")
tromso <- station.nacd("Tromsoe")
falun <- station.nacd("Falun")
stensele <- station.nacd("Stensele")
kuopio <- station.nacd("Kuopio")
valentia <- station.nacd("Valentia")
X <- combine(oslo,bergen,stockholm,copenhagen,helsinki,reykjavik,
           edinburgh,debilt,uccle,tromso,falun,stensele,kuopio,valentia)
pca <- PCA(X)
slp <- slp.NCEP(lon=c(-20,30),lat=c(30,70))
eof <- EOF(slp)
mvr <- MVR(pca,eof)
plot(mvr)

# Find the teleconnection pattern to the NAO 
data("NAOI")
data("sunspots")
data("NINO3.4")
X <- merge(NAOI,sunspots,NINO3.4,all=FALSE)

mvr <- MVR(pca,X)

# Find the pattern for NAOI:
teleconnection <- predict(mvr,newdata= c(1,0,0))
map(teleconnection,cex=2)
}
}
\keyword{manip}
