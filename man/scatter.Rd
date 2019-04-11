\name{scatter}
\alias{scatter}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Advanced scatter plots}
\description{Various functions that display bi-variate data in different ways. The default presents the number of points falling into pixles of a 2D grid. 
}
\usage{
scatter(x,y,type='heat', verbose=FALSE,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{x-variable}
  \item{y}{y-variable}
  \item{type}{Type of scatter cplot \code{c('heat','sunflower','hexbin'}}
  \item{verbose}{TRUE for diagnosing the internals of the function.}
}
\details{The 'heat' type produces a heat map type display whereas 'sunflower' and 'hexbin' present alternative infographics}
\value{
Graphics and visualisation only.
}
\author{
Kajsa M. Parding and Rasmus E. Benestad
}
\seealso{
\code{\link{vis}},\code{\link{map}},\code{\link{plot}}
}
\examples{
scatter(rnorm(10000),rnorm(10000),dig=1,log=TRUE)
scatter(rnorm(10000),rnorm(10000),type='sunflower',petalsize=7,dx=NULL,dy=NULL,
                              xgrid=NULL,ygrid=NULL,xlim=NULL,ylim=NULL,
                              xlab=NULL,ylab=NULL,main=NULL,leg=TRUE,rotate=TRUE,
                              alpha=0.6,leg.loc=2,new=TRUE, verbose=FALSE)
scatter(rnorm(10000),rnorm(10000),type='hexbin',new=TRUE,Nmax=NULL,
                           dx=NULL,dy=NULL,xgrid=NULL,ygrid=NULL,
                           xlim=NULL,ylim=NULL,
                           xlab=NULL,ylab=NULL,main=NULL,
                           leg=TRUE,col='blue',border='white',
                           colmap='heat.colors',
                           scale.col=TRUE,scale.size=FALSE, verbose=FALSE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{graphics}% use one of  RShowDoc("KEYWORDS")
\keyword{infographics}% __ONLY ONE__ keyword per line
