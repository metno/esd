\name{validate.dsensemble}
\alias{validate.dsensemble}
\alias{validate.eof.comb}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Validate
%%  ~~function to do ... ~~
}
\description{
  %%  ~~ A concise (1-5 lines) description of what the function does. ~~
  The method \code{validate} 
}
\usage{
validate.eof.comb(x, ...)
validate.dsensemble(x, conf.int=c(0.05, 0.95), 
                    colbar=list(breaks=seq(0,1,by=0.1), cex=1.5, 
                                col=colscal(11, col="t2m", alpha=0.5)),
                    plot=TRUE, verbose=FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{esd object to be validated
%%     ~~Describe \code{x} here~~
}
  \item{conf.int}{confidence interval
%%     ~~Describe \code{conf.int} here~~
}
  \item{colbar}{for plotting. See \code{colbar}
%%     ~~Describe \code{colbar} here~~
}
  \item{plot}{
%%     ~~Describe \code{plot} here~~
}
  \item{verbose}{
%%     ~~Describe \code{verbose} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
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
function (x, conf.int = c(0.05, 0.95), colbar = list(breaks = seq(0, 
    1, by = 0.1), cex = 1.5, col = colscal(11, col = "t2m", alpha = 0.5)), 
    plot = TRUE, verbose = FALSE, ...) 
{
    ranktest <- function(x) {
        obs <- zoo(attr(x, "station"))
        res <- zoo(x)
        comb <- coredata(merge(obs, res, all = FALSE))
        rs <- apply(comb, 1, function(x) order(x)[1])
        medH0 <- median(1:dim(comb)[2])
        p.val <- wilcox.test(rs, alternative = "two.sided", mu = medH0, 
            conf.int = TRUE)$p.value
        return(p.val)
    }
    if (verbose) 
        print("validate.dsensemble")
    if (inherits(x, "pca")) 
        x <- as.station(x)
    ro <- unlist(lapply(x, ranktest))
    attr(ro, "longitude") <- unlist(lapply(x, function(x) lon(attr(x, 
        "station"))))
    attr(ro, "latitude") <- unlist(lapply(x, function(x) lat(attr(x, 
        "station"))))
    attr(ro, "variable") <- "Wilcox-test score"
    attr(ro, "unit") <- "probability"
    attr(ro, "longname") <- "Wilcox-test of how observation ranks amongst model results"
    attr(ro, "history") <- history.stamp(x)
    if (is.null(colbar)) 
        colbar <- colbar.ini(ro, verbose = verbose)
    cols <- colbar$col[round(ro * length(colbar$col))]
    if (plot) {
        par0 <- par()
        par(bty = "n", fig = c(0, 1, 0, 0.85))
        plot(lon(ro), lat(ro), pch = 19, col = cols, cex = colbar$cex, 
            xlab = "", ylab = "", main = "Validation of downscaled ensemble", 
            sub = "p-values from a Wilcox-test: how observation ranks amongst model results")
        data(geoborders)
        lines(geoborders, col = "grey")
        points(lon(ro), lat(ro), cex = colbar$cex, col = "grey")
        good <- (ro > conf.int[1]) & (ro < conf.int[2])
        points(lon(ro)[good], lat(ro)[good], cex = colbar$cex, 
            lwd = 2)
        colbar(colbar$breaks, colbar$col, fig = c(0.12, 0.15, 
            0.75, 0.9))
        par(bty = "n", fig = c(0.7, 1, 0.8, 0.95), new = TRUE, 
            xaxt = "n")
        hist(ro, breaks = colbar$breaks, col = colbar$col, main = "")
        par(bty = "n", fig = c(0.25, 0.7, 0.8, 0.95), new = TRUE, 
            xaxt = "n", yaxt = "n")
        plot(c(0, 1), c(0, 1), type = "n", xlab = "n", ylab = "n")
        legend(0, 1, c(paste("ci=[", paste(conf.int, collapse = ","), 
            "]"), expression(p \%in\% ci), expression(p \%notin\% 
            ci)), pch = 21, lty = 0, lwd = c(1, 2, 1), col = c("white", 
            "black", "grey"), bty = "n")
        par(par0)
    }
    invisible(ro)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
