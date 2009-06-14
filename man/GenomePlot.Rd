\name{GenomePlot}
\alias{GenomePlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Produces a genome plot. }
\description{
Create a plot of the genome - plots each chromosome separately on the same figure.
}
\usage{
GenomePlot(x, ylim=c(-3,3), axes=FALSE, pch=46, col="black")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{	- a ".temp" formatted data stucture. }
  \item{ylim}{		- the y-axis limits defaults to c(-3,3). }
  \item{axes}{		- should axes be draw around each plot (TRUE or FALSE). }
  \item{pch}{		- the point character to use defaults to 46. }
  \item{col}{		- the point color defaults to "black". }
}
\details{
Plots as many chromosome as are contained within the ".temp" formatted data stucture on the same figure.
}
\value{
A plot showing each chromosome separatly.
}
\author{ Tomas William Fitzgerald }

\examples{

# Create a vector of chromosome numbers
chr = rnorm(100,1,0)
for (x in 2:24) { 
	chr = c(chr, rnorm(100,x,0)) 
}

# Create a ".temp" data structure
mat = matrix(ncol=8, nrow=length(chr))
mat[,1] = log2(rnorm(length(chr),1000,150)/rnorm(length(chr),800,180))
mat[,2] = rnorm(length(chr),1000,150)
mat[,3] = rnorm(length(chr),800,180)
mat[,4] = chr
mat[,5] = seq(1,length(chr), 1)
mat[,6] = seq(1,length(chr), 1)
mat[,7] = 0
mat[,8] = 0

# Plot the data
GenomePlot(mat)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ aplot }
