\name{is.pd}
\alias{is.pd}
\title{Test Positive Definiteness of a List of Square Matrices
}
\description{It tests the positive definiteness of a square matrix or a
  list of square matrices. It returns \code{FALSE} if the matrix is not
  symmetric. Variables with NA in the diagonals will be removed
  before the test.
}
\usage{
is.pd(x, tol = 1e-06)
}
\arguments{
  \item{x}{A square matrix or a list of square matrices}
  \item{tol}{Relative tolerance of positiveness of smallest eigenvalue compared
	to largest eigenvalue. The matrix is considered positive definite if the
	ratio of the smallest eigenvalue to the largest eigenvalue is larger than
	\code{tol}. See \code{\link[Matrix]{nearPD}}
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{\code{TRUE} or \code{FALSE} or a list of it.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{Mike W.-L. Cheung <mikewlcheung@nus.edu.sg>
}
\note{
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
}
\examples{
A <- diag(1,3)
is.pd(A)
# TRUE

B <- matrix(c(1,2,2,1), ncol=2)
is.pd(B)
# FALSE

is.pd(list(A, B))
# TRUE FALSE
}
\keyword{utilities}

