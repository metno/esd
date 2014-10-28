\name{InfoGraphics}
\alias{wheel}
\alias{wheel.station}
\alias{wheel.spell}

\alias{seasevol}
\alias{vis}
\alias{diagram}
\alias{cumugram}
\alias{climvar}
\title{InfoGraphics.}
\alias{colscal}
\description{
  Wheel

  Risk

  prob - boxes with forseen outcomes - area proportional to probability

  conf - confidence intervals and uncertainty - clouds...

  vis

  diagram

  cumugram
}
\usage{
vis(x,...)
diagram(x,...)
diagram.dsensemble(x,it=0,...)
wheel(x,...)
wheel.station(x,new=TRUE,lwd=2,col=NULL,bg="grey90",...)
wheel.spell(x,new=TRUE,lwd=2,col=NULL...)
cumugram(x,it=NULL,...)
climvar(x,FUN='sd',plot=TRUE,...)
colscal(cols="bwr",n=100,test=FALSE)
seasevol(x,...)
seasevol.station(x,it=NULL,nv=25,...)
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

z <- spell(bjornholt)
wheel(spell,dry=TRUE)

}
\keyword{utilities}
