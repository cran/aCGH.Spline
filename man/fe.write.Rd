\name{fe.write}
\alias{fe.write}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Function to output a Feature extraction formatted file. }
\description{
  Function used to replace intensity values, log ratio values and flags within a feature extraction input file.
}
\usage{
fe.write(file, x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file}{		- a feature extraction input file to use as a template (ideally this should be the FE file from which the data came from). }
  \item{x}{		- a ".temp" formatted dataframe containing the intensity values and flags to be written. }
}
\details{
  Using a combination of "fe.read" and "fe.write" it is easy to modifying FE data and output results back to the orignal format.  
}
\value{
This outputs a new FE file having exactly the same name as the input FE file except its filename will start with "S".\cr\cr
The new FE file will contain the intensity values, log ratio and flags from the ".temp" format input stucture(x).\cr
NB. If the input FE filename starts with an "S" the file will be overwritten.\cr
}
\author{ Tomas William Fitzgerald }
\note{ 
This function draws on methods contained within the "Jspline" class. 
}
\examples{

# Read "FEfile.txt", flag all data points and write out in FE format.
# data <- fe.read("FEfile.txt", raw=TRUE)
# data[,8] = 1
# fe.write("FEfile.txt", data)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ IO }
\keyword{ interface }
