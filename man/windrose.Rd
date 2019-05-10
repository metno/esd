\name{windrose}
\alias{windrose}
\alias{TGW}
\alias{geostrophicwind}
\alias{geostrophicwind.station}
\alias{geostrophicwind.field}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Wind analysis
}
\description{A function that plots windroses from station objects which
  contain both the zonal and meridional components. These are stored
  seperately as if they were two different station records. 
}
\usage{
TGW(triangle,f=1.25e-4,rho=1.25,verbose=FALSE)
geostrophicwind.station(x,f=1.25e-4,rho=1.25,verbose=FALSE,progressbar=TRUE)
geostrophicwind.field(x,f=1.25e-4,rho=1.25,verbose=FALSE,...)
windrose(x, saw = 10, max.scale = NULL, 
         cols = c("grey90", "yellow", "green", "red", "blue", 
                  "darkgreen", "darkred", "magenta", "black"), 
         param = c("u", "v"), simple = TRUE, verbose = FALSE)
}
\arguments{
  \item{x}{station object
}
  \item{thrshld}{Threshold for windspeed
}
  \item{saw}{Directional resolution in degrees
}
  \item{max.scale}{scaling factor for windrose
}
  \item{cols}{a vector defining colors for plot
}
  \item{param}{Name of the variables representing zonal and meridional wind
}
  \item{simple}{Only plot the windrose, not an additional historgram for windspeed
}
  \item{verbose}{if TRUE print information about progress
}
}
\author{R.E. Benestad
}
\note{Adapted from clim.pact
}

\seealso{\code{link{geostrophicwind}}
}

\examples{
\dontrun{
slp <- station(param='slp',cntr='Denmark',src='ecad')
uv <- TGW(subset(slp,is=c(1,3,11)))
UV <- geostrophicwind(slp)
windrose(uv)
map(UV,FUN='q95')
}
}
