\name{gridmap}
\alias{gridmap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Creates a griddded map
%%  ~~function to do ... ~~
}
\description{A function that uses \code{LAtticeKrieg} and elevation data to grid data and present a map.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
gridmap(Y, FUN = "mean", colbar = NULL, project = "lonlat", xlim = NULL, ylim = NULL, zlim = NULL, verbose = FALSE, plot = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
%%     ~~Describe \code{Y} here~~
}
  \item{FUN}{
%%     ~~Describe \code{FUN} here~~
}
  \item{colbar}{
%%     ~~Describe \code{colbar} here~~
}
  \item{project}{
%%     ~~Describe \code{project} here~~
}
  \item{xlim}{
%%     ~~Describe \code{xlim} here~~
}
  \item{ylim}{
%%     ~~Describe \code{ylim} here~~
}
  \item{zlim}{
%%     ~~Describe \code{zlim} here~~
}
  \item{verbose}{
%%     ~~Describe \code{verbose} here~~
}
  \item{plot}{
%%     ~~Describe \code{plot} here~~
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
function (Y, FUN = "mean", colbar = NULL, project = "lonlat", 
    xlim = NULL, ylim = NULL, zlim = NULL, verbose = FALSE, plot = FALSE) 
{
    if (verbose) 
        print(paste("gridmap", FUN))
    if (!requireNamespace("LatticeKrig", quietly = TRUE)) {
        stop("Package 'LatticeKrig' needed to use 'gridmap'. Please install it.")
    }
    else {
        if (is.null(xlim)) 
            xlim <- range(lon(Y))
        if (is.null(ylim)) 
            ylim <- range(lat(Y))
        if (!is.null(dim(Y))) {
            y <- apply(Y, 2, FUN, na.rm = TRUE)
        }
        else {
            y <- Y
        }
        if (verbose) 
            print("Use etopo5 elevation data")
        data(etopo5, envir = environment())
        etopo5 <- subset(etopo5, is = list(lon = range(lon(Y)) + 
            c(-1, 1), lat = range(lat(Y)) + c(-1, 1)))
        etopo5[etopo5 <= -1] <- NA
        if (!is.null(zlim)) {
            etopo5[(etopo5 < min(zlim)) | ((etopo5 > max(zlim)))] <- NA
        }
        if (verbose) 
            print("Use same structure as etopo5")
        grid <- structure(list(x = lon(etopo5), y = lat(etopo5)), 
            class = "gridList")
        if (verbose) 
            print("Check for duplicates")
        ok <- !(duplicated(lon(Y)) & duplicated(lat(Y)))
        if (verbose) 
            print(paste("Apply kriging to", sum(ok), "locations"))
        obj <- LatticeKrig::LatticeKrig(x = cbind(lon(Y)[ok], 
            lat(Y)[ok]), y = y[ok], Z = alt(Y)[ok])
        if (verbose) 
            print("Predict surface")
        w <- fields::predictSurface(obj, grid.list = grid, Z = etopo5)
        w$z[is.na(etopo5)] <- NA
        W <- w$z
        attr(W, "variable") <- varid(Y)[1]
        attr(W, "unit") <- unit(Y)[1]
        attr(W, "longitude") <- w$x
        attr(W, "latitude") <- w$y
        class(W) <- class(etopo5)
        if (verbose | plot) 
            print("make the map")
        map(W, xlim = xlim, ylim = ylim, zlim = zlim, colbar = colbar, 
            project = project)
        invisible(W)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
