\name{scatter.sunflower}
\alias{scatter.sunflower}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Produce a sunflower scatter plot with a hexagonal grid}
\description{Multiple points are plotted as ‘sunflowers’ with multiple leaves
     (‘petals’) such that overplotting is visualized instead of
     accidental and invisible.
}
\usage{
scatter.sunflower(x, y, petalsize = 7, dx = NULL, dy = NULL, xgrid = NULL, ygrid = NULL, xlim = NULL, ylim = NULL, xlab = NULL, ylab = NULL, leg = TRUE, rotate = TRUE, alpha = 0.6, leg.loc = 2, new = TRUE)
}
\arguments{
  \item{x, y}{the ‘x’ and ‘y’ arguments provide the x and y coordinates for
          the plot.  Any reasonable way of defining the coordinates is
          acceptable.  See the function ‘xy.coords’ for details.}
  \item{petalsize}{counts per petal of the larger sunflowers}
  \item{xgrid, ygrid}{the first row and column of the hexagonal grid. 
          If 'xgrid' is of length 2, it is interpreted as the first and last point of a grid with distance 'dx' between the intermediate points.}
  \item{xlim}{limits of the x-axis}
  \item{ylim}{limits of the y-axis}
  \item{xlab}{a label for the x-axis}
  \item{ylab}{a label for the y-axis}
  \item{leg}{logical. If TRUE, shows legend.}
  \item{rotate}{logical. If TRUE, randomly rotates petals in sunflowers}
  \item{alpha}{factor modifying the opacity alpha; typically in [0,1]}
  \item{leg.loc}{location of legend; "upper left", "upper right", "lower
    left", or "lower right"}
}

\references{
%% ~put references to the literature/web site here ~
}
\author{
 Kajsa Parding
}

\seealso{
 \code{\link{scatter.smooth}},\code{\link{scatter.hexbin}}
}

\examples{

x <- sample(seq(1,5,1e-2),5000,replace=TRUE)
y <- sample(seq(5,10,1e-2),5000,replace=TRUE)
scatter.sunflower(x,y)

x <- station.metnod(stid=39100,param="Tmax")
y <- station.metnod(stid=39100,param="Tmin")
OK <- (!is.na(x) & !is.na(y))
x <- x[OK]; y <- y[OK]
scatter.sunflower(x,y,petalsize=40,
  dx=2,dy=2,xlim=c(-20,35),ylim=c(-25,25),
  xlab="Tmax",ylab="Tmin")

}

\keyword{ scatter }
\keyword{ sunflower }% __ONLY ONE__ keyword per line
