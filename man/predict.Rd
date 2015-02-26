\name{predict.ds}
\alias{predict.ds}
\alias{predict.ds.eof}
\alias{predict.ds.com}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Prediction based on DS model
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~

}
\usage{

predict.ds(x, newdata = NULL, addnoise = FALSE, n = 100)
predict.ds.eof(x, newdata = NULL, addnoise = FALSE, n = 100)
predict.ds.comb(x, newdata = NULL, addnoise = FALSE, n = 100)
predict.mvr(object, newdata = NULL, ...) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A ds object}
%%     ~~Describe \code{x} here~~
  \item{newdata}{An eof object containing the new data sets on which
    the prediction is made.
%%     ~~Describe \code{newdata} here~~
}
  \item{addnoise}{If TRUE, will add an attribute called "noise" to the
    ouput based on WG
%%     ~~Describe \code{addnoise} here~~
}
  \item{n}{Number of runs to be generated, used only if addnoise is set to TRUE
%%     ~~Describe \code{n} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Predicted ds values.
  %%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
A. Mezghani
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{DS}} %% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# Get predictor
## Get reanalysis
X <- t2m.ERA40(lon=c(-40,50),lat=c(40,75))
## Get Gcm output
Y <- t2m.NorESM.M(lon=c(-40,50),lat=c(40,75))
## Combine
XY <- combine(X,Y)
# Compute common eof for January
ceof <- EOF(XY,it='jan')
# Get predictand
data(Oslo)
# Do the downscaling
ds <- DS(Oslo,ceof)
# Plot ds results
plot(ds)
# Do the prediction
ds.pre <- predict.ds(ds)
#Plot predicted results based on ds object
plot(ds.pre)
# Display the attribute "aspect"
attr(ds.pre, "aspect")
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
