\name{bdiagRep}
\alias{bdiagRep}
\title{Create a Block Diagonal Matrix by Repeating the Input
}
\description{It creates a block diagonal matrix by repeating the input
  matrix several times.
}
\usage{
bdiagRep(x, times)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A numeric or character matrix (or values)
}
  \item{times}{Number of times of \code{x} to be repeated
%%     ~~Describe \code{times} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{A numeric or character block diagonal matrix
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{Mike W.-L. Cheung <mikewlcheung@nus.edu.sg>
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link[metaSEM]{bdiagMat}}, \code{\link[metaSEM]{matrix2bdiag}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## Block diagonal matrix of numerics
bdiagRep( matrix(1:4,nrow=2,ncol=2), 2 )
#      [,1] [,2] [,3] [,4]
# [1,]    1    3    0    0
# [2,]    2    4    0    0
# [3,]    0    0    1    3
# [4,]    0    0    2    4

## Block diagonal matrix of characters
bdiagRep( matrix(letters[1:4],nrow=2,ncol=2), 2 )
#      [,1] [,2] [,3] [,4]
# [1,] "a"  "c"  "0"  "0" 
# [2,] "b"  "d"  "0"  "0" 
# [3,] "0"  "0"  "a"  "c" 
# [4,] "0"  "0"  "b"  "d" 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ utilities }
