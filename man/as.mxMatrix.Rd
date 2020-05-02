\name{as.mxMatrix}
\alias{as.mxMatrix}
\title{Convert a Matrix into MxMatrix-class
}
\description{It converts a matrix into \code{MxMatrix-class} via \code{mxMatrix}.
}
\usage{
as.mxMatrix(x, name, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A character or numeric matrix. If \code{x} is not a matrix,
	\code{as.matrix(x)} is applied first. 
  }
  \item{name}{An optional character string as the name of the
	MxMatrix object created by mxModel function. If the \code{name} is missing,
	the name of \code{x} will be used.}
  \item{\dots}{Further arguments to be passed to
	\code{\link[OpenMx]{mxMatrix}}. It should be noted that \code{type},
	  \code{nrow}, \code{ncol}, \code{values}, \code{free}, \code{name}
	  and \code{labels} will be created automatically. Thus, these
	  arguments except labels should be avoided in \dots 

}
}
\details{If there are non-numeric values in \code{x}, they are treated
  as the labels of the free parameters. If a "*" is
  present, the numeric value on the left-hand side will be treated as
  the starting value for a free parameter or a fixed value for a fixed parameter. If it is a matrix of numeric values, there are
  no free parameters in the output matrix. 
}
\value{A \code{\link[OpenMx]{MxMatrix-class}} object with the same
  dimensions as \code{x}
}

\author{Mike W.-L. Cheung <mikewlcheung@nus.edu.sg>
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link[OpenMx]{mxMatrix}},
	\code{\link[metaSEM]{create.mxMatrix}},
	\code{\link[metaSEM]{create.Fmatrix}},
	\code{\link[metaSEM]{checkRAM}},
    \code{\link[metaSEM]{lavaan2RAM}},
	\code{\link[metaSEM]{as.symMatrix}}
}
\examples{
## a and b are free parameters with starting values and labels
(a1 <- matrix(c(1:4, "5*a", 6, "7*b", 8, 9), ncol=3, nrow=3))
#      [,1] [,2]  [,3] 
# [1,] "1"  "4"   "7*b"
# [2,] "2"  "5*a" "8"  
# [3,] "3"  "6"   "9"  

a1 <- as.mxMatrix(a1)

## a and b are fixed parameters without any labels, name="new2"
(a2 <- matrix(1:9, ncol=3, nrow=3))
#      [,1] [,2] [,3]
# [1,]    1    4    7
# [2,]    2    5    8
# [3,]    3    6    9

new2 <- as.mxMatrix(a2, name="new2")

## Free parameters without starting values
(a3 <- matrix(c(1:4, "*a", 6, "*b", 8, 9), ncol=3, nrow=3))
#      [,1] [,2] [,3]
# [1,] "1"  "4"  "*b"
# [2,] "2"  "*a" "8" 
# [3,] "3"  "6"  "9" 

a3 <- as.mxMatrix(a3, lbound=0)

## A free parameter without label
(a4 <- matrix(c(1:4, "5*", 6, "7*b", 8, 9), ncol=3, nrow=3))
#      [,1] [,2] [,3] 
# [1,] "1"  "4"  "7*b"
# [2,] "2"  "5*" "8"  
# [3,] "3"  "6"  "9"  

a4 <- as.mxMatrix(a4)

## Convert a scalar into mxMatrix object
## "name" is required as "3*a" is not a valid name.
(a5 <- as.mxMatrix("3*a", name="a5"))
}
\keyword{utilities}
