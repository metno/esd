\name{scatter.hexbin}
\alias{scatter.hexbin}

\title{Produce a binned scatter plot with a hexagon grid}
\description{
  Multiple points are plotted as hexagons of different sizes such that
    overplotting is visualized instead of accidental and invisible.
    The size of the hexagonal markers represents the count in each grid point. 
}
\usage{
scatter.hexbin(x, y, new = TRUE, Nmax = NULL, dx = NULL, dy = NULL, xgrid = NULL, ygrid = NULL, xlim = NULL, ylim = NULL, xlab = NULL, ylab = NULL, leg = TRUE, col = "blue", border = "black")
}
\arguments{
  \item{x, y}{the ‘x’ and ‘y’ arguments provide the x and y coordinates for
          the plot.  Any reasonable way of defining the coordinates is
          acceptable.  See the function ‘xy.coords’ for details.}
  \item{Nmax}{the minimum count represented by the largest hexagonal marker.}
  \item{dx, dy}{width and height of the grid points.}
  \item{xgrid, ygrid}{the first row and column of the hexagonal grid. 
          If 'xgrid' is of length 2, it is interpreted as the first and last point of a grid with distance 'dx' between the intermediate points.}
  \item{xlim}{limits of the x-axis}
  \item{ylim}{limits of the y-axis}
  \item{xlab}{a label for the x-axis}
  \item{ylab}{a label for the y-axis}
  \item{leg}{logical. If TRUE, shows legend.}
  \item{col}{face color of markers}
  \item{border}{border color of markers}
}

\author{Kajsa Parding}

\seealso{\code{\link{scatter.sunflower}}, \code\link{{scatter.smooth}}
}

\examples{
x <- sample(seq(1,5,1e-2),5000,rep=T)
y <- sample(seq(5,10,1e-2),5000,rep=T)
scatter.hexbin(x,y)

x <- station.metnod(stid=39100,param="Tmax")
y <- station.metnod(stid=39100,param="Tmin")
OK <- (!is.na(x) & !is.na(y))
x <- x[OK]; y <- y[OK]
scatter.hexbin(x,y,Nmax=500,
  dx=2,dy=1.5,xlim=c(-20,35),ylim=c(-25,25),
  col="red",border="white",xlab="Tmax",ylab="Tmin")
}

\keyword{ scatter }
\keyword{ hexbin }

