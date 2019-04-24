\name{dX}
\alias{dX}
\alias{dY}
\alias{dT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Derivatives
%%  ~~function to do ... ~~
}
\description{Functions to estimate derivatives for gridded field objects based on a fit to truncated Fourier series. The three functions give the x-, y- and time derivatives respectively.  See Benestad & Chen (2006) 'The use of a Calculus-based Cyclone Identification method for generating storm statistics' (Tellus A 58A, 473-486, doi:10.1111/j.1600-0870.2006.00191) for more details. 
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
dX(Z, m = 10, mask.bad = TRUE, plot = FALSE, r = 6378000, chk.conf = 1, 
   accuracy = NULL, progress = TRUE, verbose = FALSE)
dY(Z, m = 10, mask.bad = TRUE, plot = FALSE, r = 6378000, chk.conf = 1, 
   accuracy = NULL, progress = TRUE, verbose = FALSE)
dT(y, m = NULL, plot = FALSE, verbose = FALSE, chk.conf = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Z}{A field object
%%     ~~Describe \code{Z} here~~
}
  \item{m}{number of harmonics for fitting the Fourier series
%%     ~~Describe \code{m} here~~
}
  \item{mask.bad}{mask missing data
%%     ~~Describe \code{mask.bad} here~~
}
  \item{plot}{for 
%%     ~~Describe \code{plot} here~~
}
  \item{r}{radius of the Earth (m)
%%     ~~Describe \code{r} here~~
}
  \item{chk.conf}{
%%     ~~Describe \code{chk.conf} here~~
}
  \item{accuracy}{
%%     ~~Describe \code{accuracy} here~~
}
  \item{progress}{show the progress intereactively
%%     ~~Describe \code{progress} here~~
}
  \item{verbose}{show diagnostics of the progress
%%     ~~Describe \code{verbose} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{a list with several comonents:

\item{Z}{original data}
\item{a}{Fourier coefficients for cosine}
\item{b}{Fourier coeffieicnes for sine}
\item{z0}{defunct?}
\item{dZ}{The component contains the first derivative.}
\item{dZ2}{The component contains the second derivative (quicker to do both in one go).}
\item{lon}{longitude}
\item{lat}{latitude}
\item{dx}{spatial resolution}
\item{span}{spatial extent}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
data(slp.ERA5)
slp.dx <- dX(slp.ERA5,verbose=TRUE)
map(slp.dx$Z) # map of SLP 
map(slp.dx$dZ) # map of first derivative in longitude direction
map(slp.dx$dZ2) # map of second derivative in longitude direction
\dontrun{
u10 <- retrieve('~/Downloads/Jan2018_ERAINT_uvp.nc',param='u10')
v10 <- retrieve('~/Downloads/Jan2018_ERAINT_uvp.nc',param='v10')
## Estimate the vorticity
zeta <- dX(v10)$dZ - dY(u10)$dZ
zeta <- attrcp(u10,zeta)
class(zeta) <- class(u10)
attr(zeta,'variable') <- 'vorticity'
attr(zeta,'unit') <- '1/s'
map(subset(zeta,it=1),projection='np')
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
