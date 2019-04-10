\name{scatter}
\alias{scatter}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Advanced scatter plots
%%  ~~function to do ... ~~
}
\description{Various functions that display bi-variate data in different ways. The default presents the number of points falling into pixles of a 2D grid. 
}
\usage{
scatter(x,y,type='heat', verbose=FALSE,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{x-variable}
  \item{y}{y-variable}
  \item{type}{Type of scatter cplot \code{c('heat','sunflower','hexbin'}
  \item{verbose}{TRUE for diagnosing the internals of the function.}
}
\details{The 'heat' type produces a heat map type display whereas 'sunflower' and 'hexbin' present the 
}
\value{
Only graphics.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Kajsa M. Parding and Rasmus E. Benestad
}
\seealso{
\code{\link{vis},\link{map},\link{plot}}
}
\examples{
scatter(rnorm(10000),rnorm(10000),dig=1,log=TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{graphics}% use one of  RShowDoc("KEYWORDS")
\keyword{infographics}% __ONLY ONE__ keyword per line
