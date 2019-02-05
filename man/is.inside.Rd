\name{is.inside}
\alias{is.inside}
\alias{test.is.inside}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{is.inside checks whether a point or a set of points is inside a polygon}
\description{
\code{is.inside} returns \code{TRUE} if the point(s) \code{x} is/are inside a polygon \code{y}. 
}
\usage{
is.inside(x, y, verbose = FALSE, plot = FALSE)
test.is.inside((N=5000,verbose=TRUE,plot=TRUE))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{an esd-object or a list/data.frame with the elements \code{x$x} and \code{x$y} containing the coordinates.}
  \item{y}{A polygon in the shape of a list/data.frame with the elements \code{x$x} and \code{x$y} containing the coordinates.}
  \item{verbose}{
\code{TRUE} prints out diagnostics for the code.
}
  \item{plot}{
\code{TRUE} provides a graphical disgnostic.
}
\item{N}{Number of tests with random coordinates}
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
function (x, y, verbose = FALSE, plot = FALSE) 
{
    if (verbose) 
        print("is.inside")
    if (plot) 
        plot(y, type = "l")
    if ((is.list(x)) | is.data.frame(x)) {
        els <- names(x)
        nn <- length(els)
        ixy <- (1:nn)[is.element(substr(tolower(els), 1, 3), 
            c("x", "lon"))][1]
        iyy <- (1:nn)[is.element(substr(tolower(els), 1, 3), 
            c("y", "lat"))][1]
        if (verbose) 
            print(c(els[ixy], els[iyy]))
        X <- data.frame(x = x[[els[ixy]]], y = x[[els[ixy]]])
    }
    else if ((!is.null(lon(x))) & (!is.null(lat(x)))) 
        X <- data.frame(x = lon(x), y = lat(x))
    if ((is.list(y)) | is.data.frame(y)) {
        els <- names(y)
        nn <- length(els)
        ixy <- (1:nn)[is.element(substr(tolower(els), 1, 3), 
            c("x", "lon"))][1]
        iyy <- (1:nn)[is.element(substr(tolower(els), 1, 3), 
            c("y", "lat"))][1]
        if (verbose) 
            print(c(els[ixy], els[iyy]))
        Y <- data.frame(x = y[[els[ixy]]], y = y[[els[ixy]]])
    }
    else if ((!is.null(lon(y))) & (!is.null(lat(y)))) 
        Y <- data.frame(x = lon(y), y = lat(y))
    d <- dim(x)
    if (verbose) 
        print(d)
    inside <- rep(FALSE, d[1])
    for (i in 1:d[1]) {
        if (verbose) 
            print(i)
        inside[i] <- is.inside.1(x[i, ], y, verbose = verbose, 
            plot = plot)
    }
    if (verbose) 
        print(table(inside))
    if (plot) {
        col <- rep("lightblue", dim(x)[1])
        col[inside] <- "darkgreen"
        points(x, pch = 19, cex = 1.2, col = col)
        text(x$x, x$y, 1:length(x$x), sub = 1, cex = 0.6, font = 2)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
