\name{hotsummerdays}
\alias{hotsummerdays}
\alias{coldwinterdays}
\alias{coldspells}
\alias{heatwavespells}
\alias{nwetdays}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Projection of hot and cold day statistics}
\description{
The functions \code{hotsummerdays}, \code{heatwavespells},
\code{coldwinterdays}, and \code{coldspells} estimate statistics for
heatwaves/hot days or cold spells based on seasonal mean
temperatures. The estimations are based on a regression analysis (GLM) between
observed number of events or spell lengths and seasonal mean from
station data. \code{nwetdays} estimates the number of days per year with
precipitation amount exceeding a threshold values.
}
\usage{
coldwinterdays(x,y=NULL,dse=NULL,it='djf',threshold=0,
               verbose=FALSE,plot=TRUE,nmin=90,new=TRUE,...)
coldspells(x,y=NULL,dse=NULL,it='jja',threshold=0,
               verbose=FALSE,plot=TRUE,nmin=90,new=TRUE,...)
hotsummerdays(x, y=NULL,dse = NULL, it = "jja", threshold = 30,
               verbose = FALSE, plot = TRUE)
heatwavespells(x,y=NULL,dse=NULL,it='jja',threshold=30,
               verbose=FALSE,plot=TRUE,new=TRUE,...)
nwetdays(x,y=NULL,dse=NULL,threshold=10,
         FUN='mu',verbose=FALSE,plot=TRUE,new=TRUE,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{station object, e.g. the temperature. Matches the element
    used in the dsensemble object 'dse'}
  \item{y}{station object which may be some statistics with dependency
    to x, e.g. snow depth.}
%%     ~~Describe \code{x} here~~
  \item{dse}{a dsensembel object. If NULL, then run DSensemble}
%%     ~~Describe \code{dse} here~~
  \item{it}{Default season set for northern hemisphere.
%%     ~~Describe \code{it} here~~
}
  \item{threshold}{Temperature threshold}
%%     ~~Describe \code{threshold} here~~
  \item{verbose}{TRUE for trouble shooting, debugging etc.}
%%     ~~Describe \code{verbose} here~~
  \item{plot}{TRUE - produce graphics}
%%     ~~Describe \code{plot} here~~
}
\details{
The estimation of these statistics makes use of general linear models
(GLMs) and take the counts to follow the 'Poisson family' whereas the
spall lengths belong to the geometric distribution. The seasonal mean
temperature or annual wet-mean precipitation are used as independent
variable.
}
\author{
R.E. Benestad
}
\examples{
data(ferder)
data(dse.ferder)
hw <- hotsummerdays(ferder,dse.ferder,threshold=20)
}
