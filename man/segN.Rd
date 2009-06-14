\name{segN}
\alias{segN}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Segmentation method. }
\description{
Java method that performs segmentation using a 'RandomWalk' algorithm.
}
\usage{
segN(ddd)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ddd}{		- vector of ratio values (sorted by genomic postion). }
}
\details{
This method segments and assesses the difference between consecutive data points using a 'RandomWalk' approach.
}
\value{
A vector of length (ddd) containing the segment medians. 
}
\author{ Tomas William Fitzgerald }
\note{ 
This method is static and has no paramethers available for the user.\cr
NB. The full method will be available as an R package soon.
}
\examples{

v = seq(1,100000,0.5)
d = sin(2*pi/500 * v) 
red = d  +  rnorm(length(d),0,100) + 1000
dd = sin(2*pi/1000 * v) 
green = dd  +  rnorm(length(dd),0,120) + 1000
 rat = log2(red / green) - median(log2(red / green), na.rm=TRUE)
 rat[20000:30000] = abs(rat[20000:30000] * 2)
 rat[60000:70000] = -abs(rat[60000:70000] * 2)
 seg  = segN(rat)
 par(mfrow=c(2,1))
 plot(rat, pch=46, ylim=c(-2,2), main="Before_segN")
 plot(rat, pch=46, ylim=c(-2,2), main="After_segN")
 matplot(seg, pch=46, ylim=c(-2,2), col='red', add=TRUE)
 legend("bottomright", bty = "n", pch=20, c("Original data values", "Segment medians"), col=c("black", "red"))
 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ manip }

