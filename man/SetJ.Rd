\name{SetJ}
\alias{SetJ}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Function to set the java class path and heap size. }
\description{
  Sets the java class path and adjusts the heap size. 
}
\usage{
SetJ()
}
%- maybe also 'usage' for other objects documented here.

\author{ Tomas William Fitzgerald }
\note{ 
 R seems to be limited at -Xmx1600, this allows FE files containing 1 million rows to be read (much more may start to cause problems). 
}

\examples{
SetJ()
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ environment }
