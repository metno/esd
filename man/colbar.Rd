\name{col.bar}
\alias{col.bar}
\alias{colbar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Color bar
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Display a color bar object on an existing plot.
}
\usage{
col.bar(breaks,horiz=TRUE,pch=21,v=1,h=1,col=col,cex=2,cex.lab=0.6,
        cex.axis=0.9,type="r",verbose=FALSE,vl=0.5,border=FALSE,...)
colbar(breaks,col,fig=c(0.15,0.2,0.15,0.3),horiz=FALSE,
       mar=c(1,0,0,0),new=TRUE,las=1,cex.axis=0.6,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{breaks}{A numeric vector of breakpoints for the colours
%%     ~~Describe \code{breaks} here~~
}
  \item{horiz}{
%%     ~~Describe \code{horiz} here~~
}
  \item{pch}{\code{\link{par}}
%%     ~~Describe \code{pch} here~~
}
  \item{v}{Vertical space between color bar points
%%     ~~Describe \code{v} here~~
}
  \item{h}{horizontal space between color bar points
%%     ~~Describe \code{h} here~~
}
  \item{col}{\code{\link{par}}
%%     ~~Describe \code{col} here~~
}
  \item{pal}{palette: c("bwr",rwb","faint.bwr","faint.rwb","rainbow","gray.colors","heat.colors","terrain.colors","topo.colors","cm.colors","grmg","brbu","budor","budrd","bugr","bugy","buor","buorr","bu","rd","cat","cold","warm")}
  \item{cex}{\code{\link{par}}
%%     ~~Describe \code{cex} here~~
}
  \item{cex.lab}{\code{\link{par}}
%%     ~~Describe \code{cex.lab} here~~
}
  \item{type}{r : rectangular shape , p : for points
%%     ~~Describe \code{type} here~~
}
  \item{verbose}{
%%     ~~Describe \code{verbose} here~~
}
  \item{vl}{Vertical lines
%%     ~~Describe \code{vl} here~~
}
  \item{border}{Color bar borders
%%     ~~Describe \code{border} here~~
}
  \item{\dots}{More graphical parameters to be passed
%%     ~~Describe \code{\dots} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
Insert a color bar or color points into an exisiting plot or map
}
\author{
A. Mezghani
  %%  ~~who you are~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (breaks, horiz = TRUE, pch = 21, v = 1, h = 1, col = col, 
    cex = 2, cex.lab = 0.6, type = "r", verbose = FALSE, vl = 0.5, 
    border = FALSE, ...) 
{
    par0 <- par()
    xleft <- par()$usr[1]
    xright <- par()$usr[2]
    ybottom <- par()$usr[4] - 1 - h
    ytop <- par()$usr[4] - 1
    by <- (xright - xleft - v * (length(col)))/(length(breaks))
    steps <- seq(0, (xright - xleft - v * (length(col))), by = by)
    nsteps <- length(steps)
    if (verbose) 
        print(steps)
    if (verbose) 
        print(breaks)
    if (verbose) 
        print(nsteps)
    k <- 1/2
    for (i in 1:(nsteps - 2)) {
        if (!is.null(v)) 
            if (i == 1) 
                k <- k + v/2
            else k <- k + v
        if (type == "r") {
            rect(xleft = k + xleft + steps[i], xright = k + xleft + 
                steps[i + 1], ybottom = ybottom, ytop = ytop, 
                col = col[i], border = border)
        }
        else if (type == "p") {
            points(x = k + xleft + (steps[i] + steps[i + 1])/2, 
                y = (ybottom + ytop)/2, pch = pch, bg = col[i], 
                cex = cex, ...)
        }
        text(x = k + xleft + (steps[i] + steps[i + 1])/2, y = ybottom - 
            vl, labels = levels(cut(breaks, breaks))[i], col = "grey50", 
            cex = cex.lab)
    }
    par(fig = par0$fig)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
