\name{vec2symMat}
\alias{vec2symMat}
\title{Convert a Vector into a Symmetric Matrix
}
\description{It converts a vector into a symmetric matrix. 
}
\usage{
vec2symMat(x, diag = TRUE, byrow = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A vector
}
  \item{diag}{Logical. If it is \code{TRUE} (the default), the diagonals of the created matrix are replaced by elements of \code{x}; otherwise, the diagonals of the created matrix are replaced by "1".
}
  \item{byrow}{Logical. If it is \code{FALSE} (the default), the created matrix is filled by columns; otherwise, the matrix is filled by rows.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{A symmetric square matrix
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
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link[metaSEM]{matrix2bdiag}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
vec2symMat(1:6)

vec2symMat(1:6, diag=FALSE)

vec2symMat(1:6, byrow=TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}
