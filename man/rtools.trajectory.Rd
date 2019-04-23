\name{anomaly.trajectory}
\alias{anomaly.trajectory}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{R-tools for trajectories
%%  ~~function to do ... ~~
}
\description{Functions to process and analyse storm trajectories
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
anomaly.trajectory(x, ..., type = "first", param = c("lon", "lat"), verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{trajectory object.
%%     ~~Describe \code{x} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
  \item{type}{
%%     ~~Describe \code{type} here~~
}
  \item{param}{
%%     ~~Describe \code{param} here~~
}
  \item{verbose}{
%%     ~~Describe \code{verbose} here~~
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
\author{Kajsa M. PArding
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
function (x, ..., type = "first", param = c("lon", "lat"), verbose = FALSE) 
{
    if (verbose) 
        print("anomaly.trajectory")
    stopifnot(!missing(x), inherits(x, "trajectory"))
    if (is.null(param)) {
        param <- unique(colnames(x)[!(colnames(x) \%in\% c("start", 
            "end", "n"))])
    }
    else if (any(!param \%in\% colnames(x))) {
        print(paste("Warning! Input error: param", paste(param[!param \%in\% 
            colnames(x)], collapse = " "), "missing"))
    }
    if (verbose) 
        print(param)
    y <- x
    y <- attrcp(x, y)
    if (any("anomaly" \%in\% attr(x, "aspect")) & any(param \%in\% 
        attr(x, "mean"))) 
        x <- anomaly2trajectory(y)
    m <- vector(mode = "list", length = length(param))
    names(m) <- param
    for (p in param) {
        i <- which(colnames(x) == p)
        xi <- x[, i]
        if (p == "lon") 
            xi <- t(apply(x[, i], 1, lontrack))
        if (type == "first") {
            m[param == p] <- list(apply(xi, 1, function(x) x[1]))
            p.anomaly <- apply(xi, 1, function(x) x - x[1])
        }
        else if (type == "mean") {
            m[param == p] <- list(mean(xi, na.rm = TRUE))
            p.anomaly <- apply(xi, 1, function(x) x - m[param == 
                p][[1]])
        }
        if (verbose) 
            print(p)
        if (verbose) 
            print(dim(x[, i]))
        if (verbose) 
            print(dim(p.anomaly))
        if (verbose) 
            print(x[1, i])
        if (verbose) 
            print(t(p.anomaly)[1, ])
        y[, i] <- t(p.anomaly)
    }
    if (any("anomaly" \%in\% attr(y, "aspect")) & any(!names(attr(y, 
        "mean")) \%in\% param)) {
        p <- names(attr(y, "mean"))
        p <- p[!p \%in\% c(param)]
        m0 <- vector(mode = "list", length = length(p))
        names(m0) <- p
        for (q in p) {
            eval(parse(text = paste("m0[param=", q, "] <- attr(y,'mean')$", 
                q, sep = "")))
        }
        m <- c(m, m0)
    }
    attr(y, "mean") <- m
    attr(y, "aspect") <- unique(c("anomaly", attr(x, "aspect")))
    attr(y, "history") <- history.stamp(x)
    invisible(y)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}% use one of  RShowDoc("KEYWORDS")
\keyword{spatial}% __ONLY ONE__ keyword per line
