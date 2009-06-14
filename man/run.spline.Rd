\name{run.spline}
\alias{run.spline}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Function to process a single file. }
\description{
Function to process a single file.
}
\usage{
run.spline(file, format = "FE", robust = TRUE, offset=5, knots=1000, ntyp = "percentile", 
	p = 0.68, fact = 4.5, segN = FALSE, sn = 0.75, writeFE = FALSE, QC = FALSE, GFF = FALSE, PDF = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{file}{	- the file to process. }
  \item{format}{	- the type of file format - c("FE", "NIM"). }
  \item{robust}{	- make the method robust to outliers (TRUE or FALSE). }
  \item{offset}{	- parameters to set the number of time the knot point will be offset. }
  \item{knots}{	- the number of knot points to use for spline fitting. }
  \item{ntyp}{		- type of noise calculation to perform - c("percentile", "derivative", "combined"). }
  \item{p}{		- value between 0 and 1 (the percentile value). }
  \item{fact}{		- numeric value (the factor by which the noise value will be rised). }
  \item{segN}{		- option for robust noise calculation (recommended when data contains many outliers e.g. cancer sample) TRUE or FALSE. }
  \item{sn}{	- value between 0 and 1 - threshold for segmentation (segN). }
  \item{writeFE}{	- Output a Feature Extraction formatted file (only to be used when processing Agilent FE files). }
  \item{QC}{	- produce a QC report TRUE or FALSE. }
  \item{GFF}{	- produce a GFF formatted file TRUE or FALSE.}
  \item{PDF}{	- produce a PDF containing normalisation figures TRUE or FALSE.}

}
\details{
In most cases running this function using defalut parameters (i.e. just providing the file input argument) will perform consistantly across a wide range of data qualities.
If method optimisation is required the number of knot point and offsets can be adjusted along with the various noise parameters. 
}
\value{
Output file containing seven columns (".temp" format):\cr
- Cy5 intensity values.\cr
- Cy3 intensity values.\cr
- Chromosome number.\cr
- Start position.\cr
- Stop position.\cr
- An index.\cr
- A flag.\cr
}
\author{ Tomas William Fitzgerald }

\examples{
# norm <- run.spline("InputFile")
# This would run the method using defalut settings.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ manip }
\keyword{ IO }
\keyword{ interface }
