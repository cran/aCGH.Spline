\name{nim.write}
\alias{nim.write}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Function to output a Nimblegen formatted file. }
\description{
  Function used to replace intensity values and log ratio values within a Nimblegen input file.
}
\usage{
nim.write(file, x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file}{		- a Nimblegen input file to use as a template (ideally this should be the Nimblegen file from which the data came from). }
  \item{x}{		- a ".temp" formatted dataframe containing the intensity values and log2 ratio to be written. }
}
\details{
  Using a combination of "nim.read" and "nim.write" it is easy to modifying Nimblegen data and output results back to the orignal format.  
}
\value{
This outputs a new Nimblegen file having exactly the same filename except it will start with "S".\cr\cr
The new FE file will contain the intensity values, log ratio and flags from the ".temp" format input stucture(x).\cr
NB. If the input filename starts with an "S" the file will be overwritten.\cr
}
\author{ Tomas William Fitzgerald }
\note{ 
This function draws on methods contained within the "Jspline" class. 
}
\examples{

# Read "Nimfile.txt", flag all data points and write out in Nimblgen format.
# data <- nim.read("FEfile.txt", raw=TRUE)
# data[,8] = 1
# nim.write("FEfile.txt", data)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ IO }
\keyword{ interface }
