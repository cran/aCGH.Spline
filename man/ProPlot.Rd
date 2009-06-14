\name{ProPlot}
\alias{ProPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Produce a simple overview plot of the data. }
\description{
Produce a figure showing the log2 ratio sorted by genomic location.
}
\usage{
ProPlot(x, pch=46, col="black", ylim=c(-3,3), xlim=c(0,length(x[,1])), xlab="Index", ylab="log2 ratio", main="Profile Plot")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{x}{	- ".temp" formatted data structure. }
  \item{pch}{		- the point character to use defaults to 46.}
  \item{col}{		- the point color to use, defaults to "black".}
  \item{xlim}{		- the x-axis limits. }
  \item{ylim}{		- the y-axis limits defaults to c(-3,3). }
  \item{xlab}{		- the label for the x-axis. }
  \item{ylab}{		- the label for the y-axis. }
  \item{main}{		- the main title for the plot. }

}
\details{
Provides an overall picture of array data.
}
\value{
Figure showing log2 ratio sorted by genomic location.
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

ProPlot(mat)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ aplot }
