\name{as.mxMatrix}
\alias{as.mxMatrix}
\title{Convert a Matrix into MxMatrix-class
}
\description{It converts a matrix into  \code{MxMatrix-class} via \code{mxMatrix}.
}
\usage{
as.mxMatrix(x, name, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A character or numeric matrix 
  }
  \item{name}{An optional character string as the name of the
	MxMatrix object created by mxModel function. If \code{name} is missing,
  the name of \code{x} will be used.}
  \item{\dots}{Further arguments to be passed to
	\code{\link[OpenMx]{mxMatrix}}. Please note that \code{type},
	  \code{nrow}, \code{ncol}, \code{values}, \code{free}, \code{name} and possibly
	  \code{labels} will be created automatically. Thus, these arguments
	  excepts labels should be avoided in \dots 
%%     ~~Describe \code{\dots} here~~
}
}
\details{If there are non-numeric values in \code{x}, these values will be treated as free parameters. If an "*" is
  present, the numeric value on the left hand side will be treated as
  the starting value for a free parameter. For example, "1" for a fixed
  parameter with "1" as the value and "5*beta" for a free parameter with
  "5" as the starting value and "beta" as the label. If it is a matrix of numeric values, there is
  no free parameters in the output matrix. 
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

\seealso{ \code{\link[OpenMx]{mxMatrix}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# a and b are free parameters with starting values and labels
a1 <- matrix(c(1:4, "5*a", 6, "7*b", 8, 9), ncol=3, nrow=3)
a1 <- as.mxMatrix(a1)

# Fixed parameters without any labels, name="new2"
a2 <- matrix(1:9, ncol=3, nrow=3)
new2 <- as.mxMatrix(a2, name="new2")

# Free parameters without starting values
a3 <- matrix(c(1:4, "*a", 6, "*b", 8, 9), ncol=3, nrow=3)
a3 <- as.mxMatrix(a3, lbound=0)

# A free parameter without label
a4 <- matrix(c(1:4, "5*", 6, "7*b", 8, 9), ncol=3, nrow=3)
a4 <- as.mxMatrix(a4)
}
\keyword{utilities}
