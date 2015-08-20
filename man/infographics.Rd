\name{InfoGraphics}
\alias{wheel}
\alias{wheel.station}
\alias{wheel.spell}
\alias{diagram}
\alias{diagram.dsensemble}
\alias{diagram.station}
\alias{seasevol}
\alias{seasevol.station}
\alias{vis}
\alias{diagram}
\alias{cumugram}
\alias{climvar}
\title{InfoGraphics.}
\alias{colscal}
\alias{prob}
\description{
  Wheel

  Risk

  prob - boxes with forseen outcomes - area proportional to probability

  conf - confidence intervals and uncertainty - clouds...

  vis

  diagram

  cumugram
  
  prob
}
\usage{
vis(x,...)
diagram(x,...)
diagram.dsensemble(x,it=0,...)
diagram.station(x,...)
wheel(x,...)
wheel.station(x,y=NULL,new=TRUE,lwd=2,col=NULL,bg="grey90",...)
wheel.spell(x,y=NULL,new=TRUE,lwd=2,col=NULL...)
cumugram(x,it=NULL,prog=FALSE,...)
climvar(x,FUN='sd',plot=TRUE,...)
colscal(n=30,col="bwr",alpha=NULL,test=FALSE)
seasevol(x,...)
seasevol.station(x,it=NULL,nv=25,...)
prob(x,...) 
prob.default(x,...)
prob.station(x,y=NULL,is=1,...)
prob.station.precip(x,y=NULL,is=1,threshold=1,...)
}
\arguments{
 \item{x}{a data object}
 }
\value{A field object}

\seealso{
  \code{\link{map}}, \code{\link{plot.station}}, \code{\link{hist.spell}}  
}

\author{R.E. Benestad and A.  Mezghanil}
\examples{
data(bjornholt)
wheel(bjornholt)

z <- spell(bjornholt,threshold=1)
wheel(z)

}
\keyword{utilities}
