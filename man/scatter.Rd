\name{scatter}
\alias{scatter}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Advanced scatter plots}
\description{Various functions that display bi-variate data in different ways. The default presents the number of points falling into pixles of a 2D grid. 
}
\usage{
scatter(x,y,type='heat', verbose=FALSE,...)
}

\arguments{
  \item{x}{x-variable}
  \item{y}{y-variable}
  \item{type}{Type of scatter cplot \code{c('heat','sunflower','hexbin'}}
  \item{verbose}{TRUE for diagnosing the internals of the function.}
  \item{breaks}{('heat' option) Set the colour scaling.}
  \item{ignorezero}{('heat' option - default = TRUE) Zeros are blank.}
  \item{log}{('heat' option - default = FALSE) Use logarithmic colour scaling.}
  \item{dig}{('heat' option) Resolution in number of decimal points/digits.}
  \item{fig}{('heat' option) Figure region for the colour bar.}
}
\details{The 'heat' type produces a heat map type display whereas 'sunflower' and 'hexbin' present alternative infographics. Scatter also takes arguments similar to plot such as 'xlim', 'ylim', 'main', 'sub', 'xlab', and 'ylab.}
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
