\name{provenance}
\alias{provenance}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Provenance
%%  ~~function to do ... ~~
}
\description{A function that extracts the history/provenance of an object. Eg. the list of steps in a chain of analysis.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
provenance(x, what = "call")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{any esd-object
%%     ~~Describe \code{x} here~~
}
  \item{what}{What to return from the history attribute
%%     ~~Describe \code{what} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{rasmus.benestad-at-met.no
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{history.stamp
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x, what = "call") 
return(unlist(attr(x, "history")[[what]]))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
