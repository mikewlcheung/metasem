\name{create.Fmatrix}
\alias{create.Fmatrix}
\title{Create an F matrix to select observed variables
}
\description{It creates an F matrix to select observed variables for wls
  function.
}
\usage{
create.Fmatrix(x, name, as.mxMatrix=TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A vector of logical type 
}
\item{name}{Name of the matrix. If it is missing, "Fmatrix" will be used.
}
\item{as.mxMatrix}{Logical. If it is \code{TRUE}, the output is a matrix of
  \code{MxMatrix-class}. If it is \code{FALSE}, it is a numeric matrix.}
\item{\dots}{Not used.}
}

\author{Mike W.-L. Cheung <mikewlcheung@nus.edu.sg>
}

\seealso{ \code{\link[metaSEM]{as.mxMatrix}}, \code{\link[metaSEM]{create.mxMatrix}}, \code{\link[metaSEM]{wls}}
}

\examples{
## Select the first 3 variables while the other 2 variables are latent.
create.Fmatrix(c(1,1,1,0,0))
# FullMatrix 'Fmatrix' 
#
# @labels: No labels assigned.
#
# @values
#      [,1] [,2] [,3] [,4] [,5]
# [1,]    1    0    0    0    0
# [2,]    0    1    0    0    0
# [3,]    0    0    1    0    0
#
# @free: No free parameters.
#
# @lbound: No lower bounds assigned.
#
# @ubound: No upper bounds assigned.

create.Fmatrix(c(1,1,1,0,0), as.mxMatrix=FALSE)
#      [,1] [,2] [,3] [,4] [,5]
# [1,]    1    0    0    0    0
# [2,]    0    1    0    0    0
# [3,]    0    0    1    0    0
}
\keyword{utilities}
