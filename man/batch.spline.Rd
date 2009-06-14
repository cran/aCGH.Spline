\name{batch.spline}
\alias{batch.spline}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Function for batch processing. }
\description{
Function for batch processing aCGH data though the "Jspline" method.
}
\usage{
batch.spline(dir, format = "FE", robust = TRUE, offset=5, knots=1000,
		ntyp = "percentile", p = 0.68, fact = 4.5, segN = FALSE, 
		sn = 0.75, writeFE = FALSE, QC = FALSE, GFF = FALSE, PDF = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dir}{		- the directory to process. }
  \item{format}{		- the type of file format - c("FE", "NIM"). }
  \item{robust}{		- make the method robust to outliers (TRUE or FALSE). }
  \item{offset}{		- parameters to set the number of time the knot point will be offset. }
  \item{knots}{		- the number of knot points to use for spline fitting. }
  \item{ntyp}{			- type of noise calculation to perform - c("percentile", "derivative", "combined"). }
  \item{p}{			- value between 0 and 1 (the percentile value). }
  \item{fact}{			- numeric value (the factor by which the noise value will be rised). }
  \item{segN}{			- option for robust noise calculation (recommended when data contains many outliers e.g. cancer sample) TRUE or FALSE. }
  \item{sn}{		- value between 0 and 1 - threshold for segmentation (segN). }
  \item{writeFE}{		- Output a Feature Extraction formatted file (only to be used when processing Agilent FE files). }
  \item{QC}{		- produce a QC report TRUE or FALSE. }
  \item{GFF}{		- produce a GFF formatted file TRUE or FALSE.}
  \item{PDF}{		- produce a PDF containing normalisation figures TRUE or FALSE.}
}
\details{
Running this function using defalut parameters (just providing the directory argument) will perform consistantly across a wide range of data qualities.\cr
\cr
If method optimisation is required, the number of knot points and offsets can be adjusted along with the various noise estimation parameters.\cr\cr 
}
\value{
As default this function outputs one ".temp" formatted file per file processed.\cr

		The ".temp" format consists of seven columns:\cr
		
		- Cy5 intensity values.\cr
		- Cy3 intensity values.\cr
		- Chromosome number.\cr
		- Start position.\cr
		- Stop position.\cr
		- An index.\cr
		- A flag.\cr
		
}
\author{ Tomas William Fitzgerald }
\note{

This function assumes that all "*.txt" files within the target directory are in the specified format.\cr
(Agilent "FeatureExtraction.txt" or Nimblegen "SEG_MNT.txt" files)\cr

}
\examples{
# batch.spline(".")
# This would run batch processing on the current directory using default settings.

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ manip }
\keyword{ IO }
\keyword{ interface }
