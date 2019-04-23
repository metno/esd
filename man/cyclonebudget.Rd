\name{cyclonebudget}
\alias{calculate.cyclonebudget}
\alias{plot.cyclonebudget}

\title{Calculate and plot the cyclone budget}

\description{Calculate and plot the cyclone budget: 
total - tracks all grid boxes visited, can be visited many times
system - tracks all grid boxes visited, but only the first visit   
genesis - record the position of the first step
lysis - record the position of the last step
outN,E,S,W - from where did the cyclone come to the present grid box?
inN,E,S,W - to where will the cyclone go from the present grid box?
}

\usage{
  calculate.cyclonebudget(X,is=NULL,it=NULL,
                          resolution.lon=12,resolution.lat=6,
                          progress=TRUE,verbose=FALSE)
}

\arguments{
  \item{X}{A 'trajectory' or 'event' object of cyclone trajectories}
  \item{it}{A list or data.frame providing time index, e.g. months, season, year range}
  \item{is}{A list providing space index, e.g., list(lon=c(-50,50),lat=c(45,70)}
  \item{resolution.lon}{Longitudinal resolution}
  \item{resolution.lat}{Latitudinal resolution}
  \item{progress}{Show progress bar. TRUE or FALSE.}
  \item{verbose}{Print out diagnostics. TRUE or FALSE.}
}

\value{A 'cyclonebudget' object: a list of various aspects of the cyclone budget.}

\author{K. Parding, MET Norway}
\examples{
\dontrun{
data(storms)
storms.deep <- trackfilter(storms,param="pcent",pmax=970,FUN="any")
storms.deep <- trackfilter(storms.deep,param="max.gradient",pmin=2.5e-2,FUN="any")
bud <- calculate.cyclonebudget(storms.deep)
plot(bud,col=colscal(n=9,name= "bu"))
}
}
\keyword{cyclonebudget}