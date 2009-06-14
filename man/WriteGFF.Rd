\name{WriteGFF}
\alias{WriteGFF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Function to write a GFF formatted file. }
\description{
Output a GFF formatted file for use in the Signal-Map software (Nimblegen software).
}
\usage{
WriteGFF(filename, x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filename}{		- the name for the output file. }
  \item{x}{	- a ".temp" formatted data stucture. }
  }
\details{
This will output a GFF formatted file containing the log2 ratio for all probes which were not flagged.
}
\value{
 GFF file for input into the Signal-Map software.
 }
\author{ Tomas William Fitzgerald }

\seealso{ \url{ http://www.nimblegen.com/products/lit/signalmap1.9usersguide.pdf } }
\examples{

chr = rnorm(1000,1,0)

for (x in 2:24) {
	chr = c(chr, rnorm(1000,x,0)) 
}

mat = matrix(ncol=8, nrow=length(chr))
mat[,1] = log2(rnorm(length(chr),1000,150) / rnorm(length(chr),800,180))
mat[,2] = rnorm(length(chr),1000,150)
mat[,3] = rnorm(length(chr),800,180)
mat[,4] = chr
mat[,5] = seq(1,length(chr), 1)
mat[,6] = seq(1,length(chr), 1)
mat[,7] = 0
mat[,8] = 0

WriteGFF("test.gff", mat) 

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ IO }
