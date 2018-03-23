\name{rainequation}
\alias{rainequation}
\alias{fract.gt.x}
\alias{test.rainequation}
\alias{scatterplot.rainequation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Rain equation
%%  ~~function to do ... ~~
}
\description{The rain equation  \eqn{Pr(X>x)=fw exp(-x/mu)}{Pr(X>x)==f_w exp(-x/\mu)} estimates the likelihood of 24-hr precipitation exceedint a threshold value x. It is analogous to the normal distribution used to describe the statistical distribution of e.g. daily temperature over a season, but applied to precipitation. It has two parameters, the wet-day frequency fw and the wet-day mean mu. It assumes that the distribution for wet-day precipitation can be approximated with an exponential distribution, and has one tail.  
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
rainequation(x, x0 = 10, threshold = NULL)
test.rainequation(loc='DE BILT',src='ecad',nmin=150,x0=20,threshold=1)
scatterplot.rainequation(src='ecad',nmin=150,x0=c(10,20,30,40))
fract.gt.x(x,x0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{x0}{
%%     ~~Describe \code{x0} here~~
}
  \item{threshold}{
%%     ~~Describe \code{threshold} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{a station object
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
Benestad R. and  A. Mezghani (2015), On downscaling probabilities for heavy 24-hr precipitation events at seasonal-to-decadal scales, Tellus A 2015, 67, 25954, http://dx.doi.org/10.3402/tellusa.v67.25954
}
\author{Rasmus Benestad
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
data(bjornholt)
plot(rainequation(bjornholt))
test.rainequation(bjornholt,threshold=30)
\dontrun{scatterplot.rainequation()}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
