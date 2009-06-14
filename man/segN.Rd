\name{segN}
\alias{segN}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Segmentation method for noise calculation. }
\description{
Segmentation method for noise calculation.
}
\usage{
segN(r, t)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{r}{		- vector of ratio values. }
  \item{t}{		- value between 0 and 1 (a percentile). }
}
\details{
This method assesses the difference between consecutive points and median normalises segments above the threshold (t).
}
\value{
A ratio which represents the noise of the original ratio(r) as if it were more normal. 
}
\author{ Tomas William Fitzgerald }
\note{ 
This method is useful when looking at highly aneuploid data (e.g. cancer samples). It can be used in an iterative manner if required.\cr
It is important to make sure the data is sorted by genomic position before running this.
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
 seg  = segN(rat,0.50)
 par(mfrow=c(2,1))
 plot(rat, pch=46, ylim=c(-2,2), main="Before_segN")
 plot(seg, pch=46, ylim=c(-2,2), main="After_segN")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ manip }

