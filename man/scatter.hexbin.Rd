\name{scatter.hexbin}
\alias{scatter.hexbin}

\title{Produce a binned scatter plot with a hexagon grid}
\description{
  Multiple points are plotted as hexagons of different sizes or colors such that
    overplotting is visualized instead of accidental and invisible.
}
\usage{
scatter.hexbin(x, y, new = TRUE, scale.col=TRUE, scale.size=TRUE, Nmax = NULL, dx = NULL, dy = NULL, xgrid = NULL, ygrid = NULL, xlim = NULL, ylim = NULL, xlab = NULL, ylab = NULL, leg = TRUE, colmap='gray.colors', col = "blue", border = "black")
}
\arguments{
  \item{x, y}{the ‘x’ and ‘y’ arguments provide the x and y coordinates for
          the plot.  Any reasonable way of defining the coordinates is
          acceptable.  See the function ‘xy.coords’ for details.}
  \item{scale.col}{logical. If TRUE, a color scale represents the count in each grid point. If FALSE, all markers are the same color.}
  \item{scale.size}{logical. If TRUE, the marker size represent the count in each grid point. If FALSE, all markers are the same size.}
  \item{colmap}{color scale of markers if scale.col is TRUE}
  \item{col}{face color of markers if scale.col is FALSE}
  \item{border}{border color of markers if scale.col is FALSE}
  \item{Nmax}{minimum count in grid points with the largest hexagonal marker. used only if scale.size is TRUE.}
  \item{dx, dy}{circumradius of the hexagonal markers in the x and y direction.}
  \item{xgrid, ygrid}{the first row and column of the hexagonal grid. 
          If 'xgrid' is of length 2, it is interpreted as the edges of a grid of hexagons of circumradius 'dx' and 'dy'.}
  \item{leg}{logical. If TRUE, shows legend.}
  \item{xlim}{limits of the x-axis}
  \item{ylim}{limits of the y-axis}
  \item{xlab}{a label for the x-axis}
  \item{ylab}{a label for the y-axis}
}

\author{Kajsa Parding}

\seealso{\code{\link{scatter.sunflower}}}

\examples{
x <- sample(seq(1,5,1e-2),5000,rep=T)
y <- sample(seq(5,10,1e-2),5000,rep=T)
scatter.hexbin(x,y)

x <- station.metnod(stid=39100,param="Tmax")
y <- station.metnod(stid=39100,param="Tmin")
OK <- (!is.na(x) & !is.na(y))
x <- x[OK]; y <- y[OK]
scatter.hexbin(x,y,
  scale.size=T,scale.col=T,Nmax=500,colmap="heat.colors",
  dx=2,dy=1.5,xlim=c(-20,35),ylim=c(-25,25),xlab="Tmax",ylab="Tmin")
}

\keyword{scatter}
\keyword{hexbin}

