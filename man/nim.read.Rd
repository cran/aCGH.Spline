\name{nim.read}
\alias{nim.read}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Function to read a Nimblegen "SEG_MNT.txt" file. }
\description{
 Function to read a Nimblgen "SEG_MNT.txt" file.
}
\usage{
nim.read(file)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file}{		- Nimblegen "SEG_MNT.txt" file to read. }
}
\details{
This function reads the contains of a Nimblegen "SEG_MNT.txt" file into a ".temp" formatted dataframe and flags low intensity probes.  
}
\value{
The ".temp" dataframe contains seven columns:\cr\cr
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
This function reads spatially normalised intensity values from the Nimblegen input file.\cr
There is currently no option to output a "SEG_MNT.txt" format, however, signal-map gff format is supported. 
}
\seealso{ \url{ http://www.nimblegen.com/ } }
\examples{
# data <- nim.read("Nimblegen_SEG_MNT.txt")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ IO }
