\name{Models}
\alias{mu.eq.f.tx}

\title{Calibrated models.}
\description{
\code{mu.eq.f.tx} contains a regression model
\eqn{mu=f(tx)}{\mu==f(t_x)}, relating the wet-day mean (mu) to the mean
daily maximum temperature - however, the input used in this model is the
saturation vapour pressure according to the Clausius-Clapeyron equation,
and 'f()' here is 'beta * C.C.eq()'.

}
\usage{
data(mu.eq.f.tx)
}

\arguments{
}

\seealso{
}


\value{lm() object.}
\author{R.E. Benestad}
\examples{
## Retrieve the model
data(mu.eq.f.tx)
## Sample data - temperature from Ferder Lighthouse
data(ferder)
pre <- data.frame(x=mean(C.C.eq(ferder),na.rm=TRUE))
## Predict the wet-day mean based on mean temperatures
predict(mu.eq.f.tx,newdata=pre)
}
\keyword{datasets}
