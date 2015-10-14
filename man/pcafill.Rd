\name{pcafill}
\alias{pcafill}
\alias{pcafill.test}
\alias{fitpc}
\title{PCA-based missing-value filling}
\description{
Fills missing (station) values by predicting their values using multiple
regression. The regression uses as input pcincipal components from PCA
from the same (group of station) data, but where series with missing
data have been excluded. This makes sense for (station) data where most
of the variability is accounted for by a few leading modes. This method
is not expected to be useful when there are many large data gaps.
}
\usage{
pcafill(X, insertmiss = 0, eofs = 1:4, test = FALSE, verbose = FALSE)
pcafill.test(X,N=100,max.miss=100, eofs = 1:4, verbose=FALSE)
fitpc(y,x,eofs=4)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{station data (group of stations)}
  \item{insertmiss}{Used for testing and evaluating. Missing data are
    introduced to test the predictive capability}
  \item{eofs}{Number of EOFs/PCAs to include in filling in. In many
    cases, it may be useful to keep this to a small set of values.}
  \item{test}{Extra test - debugging}
  \item{verbose}{Print diagnostics - debugging}
  \item{N}{Number of runs in Monte-Carlo simulation}
  \item{max.miss}{Maximum NAs to insert (insertmiss) in Monte-Carlo
    simulations}
  \item{x}{time series for calibrating regression analysis}
  \item{y}{PC input for regression analysis}
}
\details{
This function is handy for the downscaling of PCAs. See Benestad, R.E., D. Chen, A. Mezghani, L. Fan, K. Parding, On using principal components to represent stations in empirical-statistical downscaling, Tellus A 28326,  accepted.
}
\value{The same as the input - station object with filled-in values}
\references{

}
\author{Rasmus Benestad}

\seealso{
\code{\link{PCA}}, \code{\link{allgood}}
}
\examples{
download.file('http://files.figshare.com/2073466/Norway.Tx.rda',
              'Norway.Tx.rda')
load('Norway.Tx.rda')
X <- annual(Tx,FUN='mean',nmin=200)
ok<- apply(X,1,nv)
X <- subset(X,it=ok > 0)
Y <- pcafill(X)

plot(PCA(Y))
plot(c(coredata(Y)),c(coredata(X)))

## Monte-Carlo test with random selection of data points set to NA:
Y.test <- pcafill.test(X,max.miss=10,eofs=1:3)
cor(Y.test)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{PCA}
\keyword{missing data}
