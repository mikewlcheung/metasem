\name{as.mxAlgebra}
\alias{as.mxAlgebra}
\title{Convert a Character Matrix into MxAlgebra-class
}
\description{It converts a character matrix into \code{MxAlgebra} object.
}
\usage{
as.mxAlgebra(x, name="X")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A character or numeric matrix, which consists of valid
	operators in \code{mxAlgebra}. 
  }
  \item{name}{A character string of the names of the objects based
	on.}
}
\details{Suppose the name argument is "X", the output is a list of the
  following elements. 
}
\value{
  \item{X}{The \code{mxAlgebra} object.}
  \item{names}{The names of all the matrices.}
  \item{Avars}{A column vector \code{mxMatrix} of the parameters.}
  \item{Alist}{A list of mxMatrix to form the \code{mxAlgebra} object.}
}

\author{Mike W.-L. Cheung <mikewlcheung@nus.edu.sg>
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link[metaSEM]{as.mxMatrix}}, \code{\link[OpenMx]{mxAlgebra}}
}
\examples{
## a, b, and c are free parameters
(A1 <- matrix(c(1, "a*b", "a^b", "exp(c)"), ncol=2, nrow=2))
##      [,1]  [,2]    
## [1,] "1"   "a^b"   
## [2,] "a*b" "exp(c)"

A <- as.mxAlgebra(A1, name="A")

## Names of all matrices
A$names
## [1] "A"     "Avars" "A1_1"  "A2_1"  "A1_2"  "A2_2" 

## An object of mxAlgebra
A$A
## mxAlgebra 'A' 
## $formula:  rbind(cbind(A1_1, A1_2), cbind(A2_1, A2_2)) 
## $result: (not yet computed) <0 x 0 matrix>
## dimnames: NULL

## A matrix of parameters
A$Avars
## FullMatrix 'Avars' 

## $labels
##      [,1]
## [1,] "a" 
## [2,] "b" 
## [3,] "c" 

## $values
##      [,1]
## [1,]    0
## [2,]    0
## [3,]    0

## $free
##      [,1]
## [1,] TRUE
## [2,] TRUE
## [3,] TRUE

## $lbound: No lower bounds assigned.

## $ubound: No upper bounds assigned.


## A list of matrices of elements for the mxAlgebra
A$Alist
## $A1_1
## mxAlgebra 'A1_1' 
## $formula:  1 
## $result: (not yet computed) <0 x 0 matrix>
## dimnames: NULL

## $A2_1
## mxAlgebra 'A2_1' 
## $formula:  a * b 
## $result: (not yet computed) <0 x 0 matrix>
## dimnames: NULL

## $A1_2
## mxAlgebra 'A1_2' 
## $formula:  a^b 
## $result: (not yet computed) <0 x 0 matrix>
## dimnames: NULL

## $A2_2
## mxAlgebra 'A2_2' 
## $formula:  exp(c) 
## $result: (not yet computed) <0 x 0 matrix>
## dimnames: NULL
}
\keyword{utilities}
