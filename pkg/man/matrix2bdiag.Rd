\name{matrix2bdiag}
\alias{matrix2bdiag}
\title{Convert a Matrix into a Block Diagonal Matrix
}
\description{It converts a matrix into a block diagonal matrix.
}
\usage{
matrix2bdiag(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A \eqn{k}{k} x \eqn{p}{p} matrix of numerics or characters. 
}
  \item{\dots}{Further arguments to be passed to \code{\link[metaSEM]{vec2symMat}}
%%     ~~Describe \code{\dots} here~~
}
}
\details{Each row of \code{x} is converted into a symmetric matrix via
  \code{\link[metaSEM]{vec2symMat}}. Then the list of the symmetric matrices is
  converted into a block diagonal matrix via a function written by Scott
  Chasalow posted at \url{http://www.math.yorku.ca/Who/Faculty/Monette/pub/stmp/0827.html}.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
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

\seealso{ \code{\link[metaSEM]{vec2symMat}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
matrix2bdiag( matrix(1:18, ncol=6, byrow=TRUE) )

matrix2bdiag( matrix(letters[1:18], ncol=6, byrow=TRUE) )
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}
