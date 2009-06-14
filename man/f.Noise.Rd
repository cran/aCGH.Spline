\name{f.Noise}
\alias{f.Noise}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Noise calculating function. }
\description{
  Noise estimation function.
}
\usage{
f.Noise(r1, fact = 4.5, p = 0.68, typ = "percentile")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{r1}{	- numeric vector of log ratio values. }
  \item{fact}{		- numeric value (factor by which to rise the noise value). }
  \item{p}{		- value from 0 to 1, the percentile to evaluate. }
  \item{typ}{	- noise estimation type c("percentile", "difference", "combined"). }
}

\value{
An estimation of noise.
}
\author{ Tomas William Fitzgerald }

\examples{
 r = rnorm(100, 0, 3)
 f.Noise(r, fact=4.5, p=0.68, typ="percentile")
 f.Noise(r, fact=4.5, typ="derivative")
 }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ optimize }
