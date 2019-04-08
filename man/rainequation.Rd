\name{rainequation}
\alias{rainequation}
\alias{fract.gt.x}
\alias{test.rainequation}
\alias{scatterplot.rainequation}
\alias{rainvar}
\alias{rainvartrend}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Rain equation
%%  ~~function to do ... ~~
}
\description{The rain equation  \eqn{Pr(X>x)=fw exp(-x/mu)}{Pr(X>x)==f_w exp(-x/\mu)} estimates the likelihood of 24-hr precipitation exceedint a threshold value x. It is analogous to the normal distribution used to describe the statistical distribution of e.g. daily temperature over a season, but applied to precipitation. It has two parameters, the wet-day frequency fw and the wet-day mean mu. It assumes that the distribution for wet-day precipitation can be approximated with an exponential distribution, and has one tail.  

The function \code{rainvar} returns the dail variance of the 24-hr precipitation according to \eqn{sigma^2 = 2 fw * mu^3}{\sigma^2 = 2 f_w \mu^3} and \code{rainvartrend} calculates the first derivative accroding to \eqn{d sigma^2/dt = 2 mu^3 dfw/dt + 6 fw mu^2 dmu/dt}{d \frac{\sigma^2}{dt}= 2 \mu^3 \frac{d f_w}{dt} + 6 f_w * \mu^2 * \frac{d \mu}{dt}}.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
rainequation(x, x0 = 10, threshold = NULL)
test.rainequation(loc='DE BILT',src='ecad',nmin=150,
                  x0=20,threshold=1,verbose=FALSE,plot=TRUE)
scatterplot.rainequation(src='ecad',nmin=150,x0=c(10,20,30,40),
                  colour.by='x0',col=NULL,verbose=FALSE)
fract.gt.x(x,x0)
rainvar(x,x0=1,na.rm=FALSE)
rainvartrend(x,x0=1,na.rm=TRUE,mean=TRUE,nmin=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A station object - single station
%%     ~~Describe \code{x} here~~
}
  \item{x0}{The threshold value defining an event.
%%     ~~Describe \code{x0} here~~
}
  \item{threshold}{The threshold defining a 'wet day'.
%%     ~~Describe \code{threshold} here~~
}
  \item{src}{Data source}
  \item{nmin}{Minimum number of years with data}
  \item{verbose}{TRUE for printing out diagnostics}
  \item{colour.by}{Used to plot the data points with differentcolours according to e.g. 'x0', 'stid', 'alt', 'lon', or 'lat'}
  \item{col}{colour palette}
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
