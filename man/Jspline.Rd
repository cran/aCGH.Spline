\name{Jspline}
\alias{Jspline}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Spline fitting and interpolation function. }
\description{
Method to carry out robust spline fitting and interpolation.
}
\usage{
Jspline(x, offset=5, knots=1000, ntyp="percentile", p=0.68, fact=4.5, 
	robust=TRUE, segN=FALSE, sn=0.75)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{		- ".temp" formatted data stucture. }
  \item{offset}{	- numeric value 1 or greater dictating how many times to offset the knot points. }
  \item{knots}{	- the number of knot points to use in spline fitting. }
  \item{ntyp}{	- the type of noise calculation to use, c("percentile", "derivative", "combined"), defaults to "percentile". }
  \item{p}{		- numeric value between 0 and 1, the quantile to use. }
  \item{fact}{		- numeric value, the factor by which the noise estimation will be rised. }
  \item{robust}{	- make robust (exclude points from spline) TRUE or FALSE. }
  \item{segN}{		- use segmentation prior to noise estimation TRUE or FALSE (see segN). }
  \item{sn}{	- the segmentation threshold. }
}
\details{
This method carries out natural cubic spline fitting and interpolation on aCGH dual color microarray data. 
}
\value{
Input data structure (x) is returned having had cy5 and cy3 intensity data adjusted. 
}
\author{ Tomas William Fitzgerald }
\note{ 
Jspline performs consistantly over a large range of data qualities and array formats but adds most benefit to noisy, highly rearranged data.\cr\cr 
The number of points included in the spline fit can be adjusted (this ensures that the dye bias is assessed on reliable data points).\cr\cr
This method is written in java and contained within the "Jspline" class.\cr 
}
\examples{

## Set up noisy data with a bias
v = seq(1,100000,0.5)
d = sin(2*pi/500 * v) 
nd = d  +  rnorm(length(d),0,100) + 1000
dd = sin(2*pi/1000 * v) 
ndd = dd  +  rnorm(length(dd),0,120) + 1000

## Create ".temp" data stucture
mat = matrix(ncol=8,nrow=length(dd))
mat[,1] = log2(nd/ndd)
mat[,2] = nd
mat[,3] = ndd
mat[,4] = 0
mat[,5] = 0
mat[,6] = 0
mat[,7] = 0
mat[,8] = 0

## Add a few flags
mat[100:150,8] = 1

## Plot the data before and after Jspline
par(mfrow=c(2,1),mar=c(1,1,1,1))
MAPlot1(mat,ylim=c(-1,1),xlim=c(9.4,10.4))
aa = Jspline(mat)
MAPlot1(aa,ylim=c(-1,1),xlim=c(9.4,10.4))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ manip }
\keyword{ interface }
