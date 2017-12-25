\name{vec2symMat}
\alias{vec2symMat}
\title{Convert a Vector into a Symmetric Matrix
}
\description{It converts a vector into a symmetric matrix by filling up
  the elements into the lower triangle of the matrix.  
}
\usage{
vec2symMat(x, diag = TRUE, byrow = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A vector of numerics or characters
}
  \item{diag}{Logical. If it is \code{TRUE} (the default), the diagonals
	of the created matrix are replaced by elements of \code{x};
	otherwise, the diagonals of the created matrix are replaced by "1".}
  \item{byrow}{Logical. If it is \code{FALSE} (the default), the created matrix is filled by columns; otherwise, the matrix is filled by rows.
}
}

\value{A symmetric square matrix based on column major
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}

\author{Mike W.-L. Cheung <mikewlcheung@nus.edu.sg>
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link[metaSEM]{matrix2bdiag}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
vec2symMat(1:6)
#      [,1] [,2] [,3]
# [1,]    1    2    3
# [2,]    2    4    5
# [3,]    3    5    6

vec2symMat(1:6, diag=FALSE)
#      [,1] [,2] [,3] [,4]
# [1,]    1    1    2    3
# [2,]    1    1    4    5
# [3,]    2    4    1    6
# [4,]    3    5    6    1

vec2symMat(letters[1:6])
#      [,1] [,2] [,3]
# [1,] "a"  "b"  "c" 
# [2,] "b"  "d"  "e" 
# [3,] "c"  "e"  "f" 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}
