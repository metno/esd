\name{track}
\alias{Cylcone tracking algorithm}
\alias{track.default}
\alias{track.events}
\alias{Track}
\alias{Trackstats}
\title{3-step cyclone tracking algorithm.}
\description{
Applies a tracking algorithm to a set of cyclones (\code{\link{CCI}}).

The algorithm connects events in three subsequent time steps, chosing
the path that minimizes the change in angle and displacement between them. 
The analysis can be applied to 'events' objects.

The maximum displacement ('dmax') between two time steps depends on the 
angle of direction. The circle of maximum displacement around a cyclone can be 
shifted in the east-west ('dE') and north-south ('dN') direction. For example, 
'dE' = 0.3 means that the maximum displacement is 1.3*dmax in the eastward direction 
and only 0.7*dmax in the westward direction, which is in line with the general 
westerly flow in the storm track regions. 

The maximum change in the angle of direction ('amax') is enforced for medium to long displacement, 
but for short displacements (< dmax/4) there is no restriction of the change of direction.

Note: The algorithm has been developed for tracking midlatitude cyclones
in the northern hemisphere and may not work as well for other
regions or 'events' of different types, e.g., anti-cyclones. 
}
\usage{
track(X,x0=NULL,it=NULL,is=NULL,dmax=1E6,amax=90,
             nmax=31*24,nmin=5,dmin=5E5,dE=0.3,dN=0.2,
             lplot=FALSE,progress=TRUE,verbose=FALSE)
track.events(X,...) 

}
\arguments{
    \item{X}{An 'events' object containing temporal and spatial information about a set of cyclones or anticyclones.}
    \item{x0}{A tracked 'events' object from previous time steps, used as a starting point for the tracking 
              of X so that trajectories can continue from x0 to X.}
     \item{it}{A list providing time index, e.g. month.}
    \item{is}{A list providing space index, lon and/or lat.}
    \item{dmax}{Maximum displacement of events between two time steps in a trajectory. Unit: m.}
    \item{amax}{Maximum change in the angle of direction. Unit: degrees.}
    \item{nmax}{Maximum lifetime of a trajectory. Unit: hours.}
    \item{nmin}{Minimum lifetime of a trajectory. Unit: hours.}
    \item{dmin}{Minimum displacement between the start and end point of a trajectory. Unit: m.}
    \item{dE}{Shift of the maximum displacement in the eastward direction.}
    \item{dN}{Shift of the maximum displacement in the northward direction.}
    \item{lplot}{TRUE: Produce plots of trajectories for selected time steps.}
    \item{progress}{TRUE: Show progress bar.}
    \item{verbose}{TRUE: Print out diagnosics.}  
 }
\value{An 'events' object containing the original information as well as 
  the trajectory number ('trajectory') of each event and statistical properties 
  of the trajectories ('trackcount' - number of events in path;
  'tracklen' - distance between start and end point of path').}
\author{K. Parding}
\examples{

}
\keyword{track}
