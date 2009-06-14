\name{MakePDF}
\alias{MakePDF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Write a PDF containing some normalisation plots. }
\description{
This function creates a PDF containing various figures to show the quality of normalisation.
}
\usage{
MakePDF(filename, x, raw, ntyp="percentile", p=0.68, fact=4.5, segN=FALSE, sn=0.75)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filename}{		- character argument, the output file name. }
  \item{x}{	- the ".temp" formatted data structure after normalisation. }
  \item{raw}{		- the".temp" formatted data structure before normalisation.}
  \item{ntyp}{	- the type of noise calculation to use, c("percentile", "derivative", "combined"). }
  \item{p}{		-  numeric value between 0 and 1, the quantile used.}
  \item{fact}{		- numeric value, the factor by which the noise estimation was rised.}
  \item{segN}{		- was segN used TRUE or FASLE. }
  \item{sn}{	- the segmentation threshold. }
}
\details{
Produces a summary ".pdf" files containing various figures.
}
\value{
PDF containing a "ProPlot", "GenomePlot" and two MAplots (MAplot1 and MAplot2).
}
\author{ Tomas William Fitzgerald }
\note{ 
MA-plots a usful way of assessing the quality of the normalisation. Most points should be centered around zero 
and the outliers should have little effect on the overall fit (this is the key point of this normalisation method). 
}
\examples{
# NA
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ IO }
