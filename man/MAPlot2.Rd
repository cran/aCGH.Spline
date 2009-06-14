\name{MAPlot2}
\alias{MAPlot2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Produce a MAplot showing points fitted vs. points interpolated. }
\description{
Produce a MA-plot showing points fitted and point interpolated.
}
\usage{
MAPlot2(x, raw, ntyp="percentile", p=0.68, fact=4.5, segN=FALSE, sn=0.75, pch=46, Fitcol="green", Intcol="red", ylim=c(-5,5), xlim=c(0,20), xlab="log2(mean(cy5, cy3))", ylab="log2(cy5/cy3)", main="MA-plot")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{	- the ".temp" formatted data structure after normalisation. }
  \item{raw}{		- the".temp" formatted data structure before normalisation.}  \item{ntyp}{	- the type of noise calculation to use, c("percentile", "derivative", "combined"), defaults to "percentile". }
  \item{p}{		- numeric value between 0 and 1, the quantile to use. }
  \item{fact}{		- numeric value, the factor by which the noise estimation will be rised. }
  \item{segN}{		- use segmentation prior to noise estimation TRUE or FALSE (see segN). }
  \item{sn}{	- the segmentation threshold. }
  \item{pch}{		- the point character to use defaults to 46.}
  \item{Fitcol}{		- the point color to use for the points used during spline fitting, defaults to "green".}
  \item{Intcol}{		- the point color to use for the points excluded from spline fitting, defaults to "red".}
  \item{xlim}{		- the x-axis limits defaults to c(0,20). }
  \item{ylim}{		- the y-axis limits defaults to c(-5,5). }
  \item{xlab}{		- the label for the x-axis. }
  \item{ylab}{		- the label for the y-axis. }
  \item{main}{		- the main title for the plot. }
}

\value{
MAplot showing relative number of points considered to be outside the normal distribution.
}
\author{ Tomas William Fitzgerald }
\note{ 
Points outside the interpolation threshold are not used during the spline fitting.
We believe that these points are highly unrealiable for assessing the dye bias.
}
\examples{
# Create some noisy data with a bias
v = seq(1,100000,0.5)
d = sin(2*pi/500 * v) 
nd = d  +  rnorm(length(d),0,100) + 1000
dd = sin(2*pi/1000 * v) 
ndd = dd  +  rnorm(length(dd),0,120) + 1000

## Create ".temp" data stucture
mat = matrix(ncol=7,nrow=length(dd))
mat[,1] = nd
mat[,2] = ndd
mat[,3] = 0
mat[,4] = 0
mat[,5] = 0
mat[,6] = 0
mat[,7] = 0

# Add a few outliers
mat[100:150,1] = mat[100:150,1] * 1.5
mat[100:150,2] = mat[100:150,2] / 1.5

# Run Jspline and plot
aa = Jspline(mat)

MAPlot2(aa, mat, ylim=c(-1.5,1.5),xlim=c(9.4,10.4))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ aplot }
