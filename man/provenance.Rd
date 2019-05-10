\name{provenance}
\alias{provenance}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Provenance
%%  ~~function to do ... ~~
}
\description{A function that extracts the history/provenance of an object. Eg. the list of steps in a chain of analysis.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
provenance(x, what = "call")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{any esd-object
%%     ~~Describe \code{x} here~~
}
  \item{what}{What to return from the history attribute
%%     ~~Describe \code{what} here~~
}
}

\seealso{history.stamp}

\examples{
data(Oslo)
provenance(Oslo)
}