\name{WriteQC}
\alias{WriteQC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Function to produce a QC report. }
\description{
Write a QC report containing some summery statistics
}
\usage{
WriteQC(filename, batch, x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filename}{		- the name for the output file. }
  \item{batch}{		- TRUE or FALSE (is it being run in batch mode). }
  \item{x}{	- a ".temp" formatted data stucture. }
 }

\value{
 Produces a table containing the filename, 68th percentile, dLRs and chromosme X median.
 }

\author{ Tomas William Fitzgerald }

\examples{

chr = rnorm(1000,1,0)

for (x in 2:24) {
	chr = c(chr, rnorm(1000,x,0)) 
}

mat = matrix(ncol=7, nrow=length(chr))
mat[,1] = rnorm(length(chr),1000,150)
mat[,2] = rnorm(length(chr),800,180)
mat[,3] = chr
mat[,4] = seq(1,length(chr), 1)
mat[,5] = seq(1,length(chr), 1)
mat[,6] = 0
mat[,7] = 0

WriteQC("test", batch=FALSE, mat)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ IO }
