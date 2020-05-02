\name{as.symMatrix}
\alias{as.symMatrix}
\title{Convert a Character Matrix with Starting Values to a Character Matrix
  without Starting Values
}
\description{It converts a character matrix with starting values to a
  character matrix without the starting values.
}
\usage{
as.symMatrix(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A character or numeric matrix or a list of character or
	numeric matrices.}
}
\details{If there are non-numeric values in \code{x}, they are treated
  as the labels of the free parameters. If a "*" is
  present, the numeric value on the left-hand side will be treated as
  the starting value for a free parameter or a fixed value for a fixed parameter. If it is a matrix of numeric values, there are
  no free parameters in the output matrix. This function removes the
  starting values and "*" in the matrices. 
}
\value{A character matrix.
}

\author{Mike W.-L. Cheung <mikewlcheung@nus.edu.sg>
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link[metaSEM]{as.mxMatrix}}
}
\examples{
## a and b are free parameters with starting values and labels
(a1 <- matrix(c(1:4, "5*a", 6, "7*b", 8, 9), ncol=3, nrow=3))
#      [,1] [,2]  [,3] 
# [1,] "1"  "4"   "7*b"
# [2,] "2"  "5*a" "8"  
# [3,] "3"  "6"   "9"  

(as.symMatrix(a1))
#      [,1] [,2] [,3]
# [1,] "1"  "4"  "b" 
# [2,] "2"  "a"  "8" 
# [3,] "3"  "6"  "9"
}
\keyword{utilities}
