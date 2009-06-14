\name{MAPlot1}
\alias{MAPlot1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Produce a MAplot showing points retained vs. points excluded. }
\description{
 Produce a MA-plot showing points retained and points excluded.
}
\usage{
MAPlot1(x, pch=46, Retcol="black", Excol="red", ylim=c(-5,5), xlim=c(0,20), xlab="log2(mean(cy5, cy3))", ylab="log2(cy5/cy3)", main="MA-plot")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{	- ".temp" formatted data structure. }
  \item{pch}{		- the point character to use defaults to 46.}
  \item{Retcol}{		- the point color to use for the points retained after normalisation, defaults to "black".}
  \item{Excol}{		- the point color to use for the points excluded after normalisation, defaults to "red".}
  \item{xlim}{		- the x-axis limits defaults to c(0,20). }
  \item{ylim}{		- the y-axis limits defaults to c(-5,5). }
  \item{xlab}{		- the label for the x-axis. }
  \item{ylab}{		- the label for the y-axis. }
  \item{main}{		- the main title for the plot. }
}
\value{
Figure showing the log2(mean(red,green)) vs. the log2(red/green). 
}
\author{ Tomas William Fitzgerald }

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
mat[100:150,7] = 1

## Plot the data
par(mfrow=c(2,1),mar=c(1,1,1,1))
MAPlot1(mat,ylim=c(-1,1),xlim=c(9.4,10.4))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ aplot }

