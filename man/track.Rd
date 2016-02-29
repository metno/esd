\name{track}
\alias{Cylcone tracking algorithm}
\alias{track.default}
\alias{track.events}
\alias{track.trajectory}
\alias{Track}
\alias{Trackstats}
\title{3-step cyclone tracking algorithm}
\description{
Applies a tracking algorithm to a set of cyclones (\code{\link{CCI}}).

The algorithm connects events in three subsequent time steps, chosing
the path that minimizes the change in angle and displacement between them. 

The analysis can be applied to 'events' objects as well as 'trajectory'
objects, in the latter case producing new trajectories.

Note: the analysis has been developed for tracking midlatitude cyclones
in the northern hemisphere. The algorithm may not work as well for other
regions or 'events' of different types, e.g., anti.cyclones. Try
tweeking the input parameters 'dmin', 'dmax', 'amax', etc.

}
\usage{
track(X,...)
track.default(X,...) 
track.events(X,...)
track.trajectory(X,...)

}
\arguments{
    \item{X}{An object with cyclones, anticyclones or other events.}
  
 }
\value{An events object containing trajectory number ('trajectory') and
  statistical properties ('trackcount - number of events in path;
  tracklen - distance between start and end point of path').}
\author{R.E. Benestad}
\examples{

}
\keyword{track}
