\name{fe.extract}
\alias{fe.extract}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Function to extract Feature Extraction data. }
\description{
 Read Agilent Feature extraction data into R.
}
\usage{
fe.extract(file, list)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file}{		- a feature extraction file to read. }
  \item{list}{	- a list of column names to extract e.g. c("rBGSubSignal", "gBGSubSignal", "SystematicName"). }
}
\details{
Function to extract data from a Feature extraction formatted file.
}
\value{
Returns a dataframe containing all the values speicfied in the list argument. 
}

\author{ Tomas William Fitzgerald }
\note{ 
This function is capable of extracting any data contained within a Feature extraction file. 
}
\seealso{ \url{ http://www.chem.agilent.com/en-us/products/instruments/dnamicroarrays/featureextractionsoftware/pages/default.aspx } }
\examples{
# data <- fe.extract("FEfile.txt", list = c("rBGSubSignal", "gBGSubSignal", "SystematicName")
# This would return and data frame containing background subtracted intensity values and positional information.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ IO }
