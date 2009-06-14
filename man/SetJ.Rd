\name{SetJ}
\alias{SetJ}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Function to set the java class path and heap size. }
\description{
  Sets the java class path and adjusts the heap size. 
}
\usage{
SetJ(JavaHeapSize = "-Xmx1000m")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{JavaHeapSize}{		- argument to adjust the heap size - defualts to 1000Mb. }
 }
\author{ Tomas William Fitzgerald }
\note{ 
 R seems to be limited at -Xmx1600, this allows FE files containing 1 million rows to be read (much more may start to cause problems).\cr
 This function needs to be run first to adjust the heap (heap size can only be set once within an R session).\cr
 NB. if either "fe.read", "fe.write" or "Jspline" are run before "SetJ()" the heap will be set to default (1000Mb).
}

\examples{
SetJ()
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ environment }
