\name{is.pd}
\alias{is.pd}
\title{Test Positive Definiteness of a List of Square Matrices
}
\description{It tests the positive definiteness of a square matrix or a
  list of square matrices. It returns \code{TRUE} if the matrix is
  positive definite. It returns \code{FALSE} if the matrix is either
  non-positive definite or not symmetric. Variables with \code{NA} in the diagonals will be removed
  before testing. It returns \code{NA} when there are missing correlations even after deleting
  the missing variables.
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
\value{If the input is a matrix, it returns \code{TRUE}, \code{FALSE}
  or \code{NA}. If the input is a list of matrices, it returns
  a list of \code{TRUE}, \code{FALSE} or \code{NA}.
}

\author{Mike W.-L. Cheung <mikewlcheung@nus.edu.sg>
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

C <- A
C[2,1] <- C[1,2] <- NA
is.pd(C)
# NA
}
\keyword{utilities}

