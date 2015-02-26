\name{predict.ds}
\alias{predict.ds}
\alias{predict.ds.eof}
alias{predict.ds.com}
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
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{newdata}{
%%     ~~Describe \code{newdata} here~~
}
  \item{addnoise}{
%%     ~~Describe \code{addnoise} here~~
}
  \item{n}{
%%     ~~Describe \code{n} here~~
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
\author{
%%  ~~who you are~~
A. Mezghani
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x, newdata = NULL, addnoise = FALSE, n = 100) 
{
    stopifnot(!missing(x), inherits(x, "ds"))
    if ((inherits(x, "eof")) & (is.null(newdata))) {
        if (inherits(x, "comb")) 
            y <- predict.ds.comb(x, newdata = newdata, addnoise = addnoise, 
                n = n)
        else if (inherits(x, "field")) 
            y <- predict.ds.eof(x, newdata = newdata, addnoise = addnoise, 
                n = n)
    }
    else if (inherits(x, "eof")) {
        if (inherits(newdata, "comb")) 
            y <- predict.ds.comb(x, newdata = newdata, addnoise = addnoise, 
                n = n)
        else if (inherits(newdata, "eof")) 
            y <- predict.ds.eof(x, newdata = newdata, addnoise = addnoise, 
                n = n)
    }
    y <- attrcp(attr(x, "original_data"), y)
    attr(y, "history") <- history.stamp(x)
    class(y) <- class(x)
    invisible(y)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
