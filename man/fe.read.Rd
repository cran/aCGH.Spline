\name{fe.read}
\alias{fe.read}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Function to read Feature Extraction files. }
\description{
Recommended feature extraction file reading function.
}
\usage{
fe.read(file, Raw=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file}{		- the feature extraction file to read. }
  \item{Raw}{		- the intensity values to extract (TRUE = background subtracted values, FALSE = Agilent normalised values). }
}
\details{
This function is capable of reading a FE file containing 1 million rows of data within 2 minutes.\cr\cr 
Its also flags data points based on Feature extraction exculsion criteria and low intensity probes (cy5 + cy3 < 100).\cr
}
\value{
A dataframe containing eight columns (".temp" format):\cr\cr
- Log2 ratio values.\cr
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
If the flag is equal to 1 the probe should be removed from further analysis.\cr\cr
The index can be used as a unique indentifier, it is the extact row where the data was taken from the FE file.\cr\cr
The function draws on methods contained within the "Jspline" class for increased speed performance.\cr
}
\examples{
# data <- FeRead("FEfile.txt")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ IO }
\keyword{ interface }
