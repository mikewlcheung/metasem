\name{is.mxMatrix}
\alias{is.mxMatrix}
\title{Test whether the object is of MxMatrix-class
}
\description{It tests whether the object is of \code{MxMatrix-class}.
}
\usage{
is.mxMatrix(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{An object 
  }
 \item{\dots}{Not used
}
}

\value{A \code{\link[OpenMx]{MxMatrix-class}} object with the same
  dimensions as \code{x}
}
\references{
}
\author{Mike W.-L. Cheung <mikewlcheung@nus.edu.sg>
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link[OpenMx]{mxMatrix}},
	\code{\link[metaSEM]{create.mxMatrix}},
	\code{\link[metaSEM]{create.Fmatrix}}, \code{\link[metaSEM]{as.mxMatrix}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## a and b are free parameters with starting values and labels
(a1 <- matrix(c(1:4, "5*a", 6, "7*b", 8, 9), ncol=3, nrow=3))
#      [,1] [,2]  [,3] 
# [1,] "1"  "4"   "7*b"
# [2,] "2"  "5*a" "8"  
# [3,] "3"  "6"   "9"  

is.mxMatrix(a1)
# [1] FALSE

a1 <- as.mxMatrix(a1)
is.mxMatrix(a1)
# [1] TRUE
}
\keyword{utilities}
