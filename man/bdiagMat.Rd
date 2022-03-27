\name{bdiagMat}
\alias{bdiagMat}
\title{Create a Block Diagonal Matrix
}\description{It creates a block diagonal matrix from a list of numeric
  or character matrices.
}
\usage{
bdiagMat(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A list of numeric or character matrices (or values)
}
}

\value{A numeric or character block diagonal matrix
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{It was based on a function posted by Scott Chasalow at http://www.math.yorku.ca/Who/Faculty/Monette/pub/stmp/0827.html.
%% ~put references to the literature/web site here ~
}
\author{Mike W.-L. Cheung <mikewlcheung@nus.edu.sg>
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link[metaSEM]{bdiagRep}}, \code{\link[metaSEM]{matrix2bdiag}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## Block diagonal matrix of numbers
bdiagMat( list(matrix(1:4,nrow=2,ncol=2),
               matrix(5:6,nrow=1,ncol=2)) )
#      [,1] [,2] [,3] [,4]
# [1,]    1    3    0    0
# [2,]    2    4    0    0
# [3,]    0    0    5    6

## Block diagonal matrix of characters
bdiagMat( list(matrix(letters[1:4],nrow=2,ncol=2),
               matrix(letters[5:6],nrow=1,ncol=2)) )
#      [,1] [,2] [,3] [,4]
# [1,] "a"  "c"  "0"  "0" 
# [2,] "b"  "d"  "0"  "0" 
# [3,] "0"  "0"  "e"  "f" 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ utilities }
