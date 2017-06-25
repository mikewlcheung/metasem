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
is.pd(x, check.aCov=FALSE, cor.analysis=TRUE, tol=1e-06)
}
\arguments{
  \item{x}{A square matrix or a list of square matrices}
  \item{check.aCov}{If it is \code{TRUE}, it mirrors the checking in
	\code{\link[metaSEM]{asyCov}}.}
  \item{cor.analysis}{Whether the input matrix is a correlation or a
	covariance matrix. It is ignored when \code{check.aCov=FALSE}.}
  \item{tol}{Tolerance (relative to largest variance) for numerical lack
	of positive-definiteness in \code{x}. It is adopted from \code{\link[MASS]{mvrnorm}}.
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

