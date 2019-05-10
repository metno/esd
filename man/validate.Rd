\name{validate.dsensemble}
\alias{validate.dsensemble}
\alias{validate.eof.comb}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Validate
}
\description{
  The method \code{validate} 
}
\usage{
validate.eof.comb(x, ...)
validate.dsensemble(x, conf.int=c(0.05, 0.95), 
                    colbar=list(breaks=seq(0,1,by=0.1), cex=1.5, 
                                col=colscal(11, col="t2m", alpha=0.5)),
                    plot=TRUE, verbose=FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{esd object to be validated
}
  \item{conf.int}{confidence interval
}
  \item{colbar}{for plotting. See \code{colbar}
}
  \item{plot}{if TRUE produce plot
}
  \item{verbose}{if TRUE print progress
}
  \item{\dots}{other arguments
}
}

\examples{
  slp1 <- slp.DNMI(lon=c(-50,50),lat=c(30,70))
  slp2 <- slp.DNMI(lon=c(-50,50),lat=c(30,70))
  slpcomb <- combine(slp1,slp2)
  eofcomb <- EOF(slpcomb)
  validate(eofcomb)
}
