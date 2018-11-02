\name{kge}
\alias{KGE}
\alias{Kling-Gupta}
\alias{nse}
\alias{NSE}
\alias{Nash-Sutcliffe}
\title{Model efficiency evaluation}
\description{The KGE function returns the 2012 Kling-Gupta model efficiency using CV for the variability of x in predicting y. 
The NSE function returns the dimensionless Nash Sutcliffe model efficiency, evalutaing x as a prediction of y.
An efficiency of one corresponds to a perfect match, while the lowest score is -infinity.}
\usage{
  KGE(x,y)
  NSE(x,y)
  kge(x,y,return_all=FALSE, scc=1., sv=1.,sm=1.)
  nse(x,y)
}
\arguments{
  \item{x}{a 1-D zoo object corresponding to the prediction}
  \item{y}{a 1-D zoo object corresponding to the reference}
  \item{return_all}{Directional resolution in degrees}
  \item{scc}{weight of (0-1)}
  \item{sv}{weight of (0-1)}
  \item{sm}{weight of (0-1)}
}
\details{The NSE function returns the multi-objective Nash-Sutcliffe efficiency metric, 
  which corresponds to the unbiased \eqn{R^2}. The metric is scaled by the observed variance. 
  The metric is given as nse = 1 - mean_square_error/observed_variance. 
  It is calculated for pairwise complete observations

  Warning: Keep in mind that in regions with higher variance (e.g. seasonality) equal absolute deviations will 
  be penalized less than for regions with a lower variance.

  The KGE function returns the dimensionless, multi-objective Kling-Gupta efficiency metric, 
  which is based on the correlation, variability ratio, and mean ratio
  between the prediction and reference objects and is a modified version of the
  Nash-Sutcliffe efficiency. 
  The weight of the sub-metrics (scc, sv, sm) may be adjusted from 1 in the arguments
  kge = 1 - sqrt( (scc(cc-1)^2 + sv(CV_m/CV_o - 1)^2 + sm(mean(mod)/mean(obs)-1)^2 )).
  
  Warning: The mean ratio may not be approriate for Celcius (can blow up around 0)
  or other variables that may have a mean within (-1,1). Though the KGE etric is dimensionless, 
  keep in mind that variables with higher averages (e.g. temperature in Kelvin, longwave radiation) 
  will be less penalized for mean deviations than variables with lower averages 
  (temp. in Celsius, winter SW radiation).
}
\references{
Gupta HV, Kling H, Yilmaz KK and Martinez GF (2009)., "Decomposition of the mean squared error and NSE performance criteria: Implications for improving hydrological modelling.", Journal of Hydrology, 377(1), 80-91.
Nash JE and Sutcliffe JV (1970), River flow forecasting through conceptual models part I-A discussion of principles.???, Journal of hydrology, 10(3), pp. 282-290.
Gupta HV, Kling H, Yilmaz KK and Martinez GF (2009)., "Decomposition of the mean squared error and NSE performance criteria: Implications for improving hydrological modelling.", Journal of Hydrology, 377(1), pp. 80-91.}
\author{H.B. Erlandsen}

\examples{
  \dontrun{
    data('Oslo')
    x <- subset(Oslo,it=!is.na(Oslo))
    mydata <- data.frame(x=x,t=index(x))
    fit <- lm(x ~ t, data=mydata)
    y <- zoo(fitted(fit,index(x)),order.by=index(x))
    KGE(x,y)
  }
}

\keyword{"validation"}
\keyword{"model-efficiency"}