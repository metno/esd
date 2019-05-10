\name{anomaly.trajectory}
\alias{anomaly.trajectory}
\alias{polyfit.trajectory}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{R-tools for trajectories
%%  ~~function to do ... ~~
}
\description{Functions to process and analyse storm trajectories
}
\usage{
anomaly.trajectory(x, ..., type = "first", param = c("lon", "lat"), verbose = FALSE)
anomaly2trajectory(x, verbose=FALSE)
polyfit.trajectory(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A trajectory object}
  \item{\dots}{Other arguments}
  \item{type}{type of anomaly: 'first' gives you the spatial anomaly with regards to the first time step of the trajectory and 'mean' centers the trajectories around the mean longitude and latitude}
  \item{param}{parameters to calculate anomaly of}
  \item{verbose}{if TRUE print information about progress}
 }
  \value{A trajecory object}
\author{Kajsa M. Parding}

\examples{
# Load trajectory example data
data('imilast.M03')
# Calculate anomaly
x <- anomaly.trajectory(imilast.M03)
# Show maps of original trajectories and spatial anomalies
map(imilast.M03)
map(x)
# Transform trajectory anomalies back to regular trajectories
y <- anomaly2trajectory(x)
# Print longitudes of first trajectory
imilast.M03[1,1:12]
x[1,1:12]
y[1,12]

# Fit polynomial to trajectories
y <- polyfit.trajectory(x)
# Show coefficients of first trajectory
print(attr(y,"coefs")[,1])
# Plot original trajectory and polynomial fit
ilon <- colnames(x)=="lon"
ilat <- colnames(x)=="lat"
plot(x[1,ilon], x[1,ilat], col="black", pch=1)
lines(y[1,ilon], y[1,ilat], col="blue", lty=2)
}