\name{as.mxMatrix}
\alias{as.mxMatrix}
\title{Convert matrices into MxMatrix-class
}
\description{It converts a matrix into  \code{MxMatrix-class} via \code{mxMatrix}.
}
\usage{
as.mxMatrix(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A matrix of characters or numeric values. If it is a matrix
	of numeric values, there is no free parameters in the output matrix.
}
  \item{\dots}{Further arguments to be passed to
	\code{\link[OpenMx]{mxMatrix}}. Please note that \code{type},
	  \code{nrow}, \code{ncol}, \code{values}, \code{free} and possibly
	  \code{labels} will be created automatically. Thus, these values
	  should not be passed through \dots
%%     ~~Describe \code{\dots} here~~
}
}
\details{The dimensions are the same as \code{x}. If there are
  non-numeric values, it will be treated as a free parameter. If an "*" is
  present, the numeric value on the left hand side will be treated as
  either a fixed value for a fixed parameter or a starting value for a
  free parameter. For example, "1" for a fixed parameter with "1" as the
  value and "5*beta" for a free parameter with "beta" as the label.
%%  ~~ If necessary, more details than the description above ~~
}
\value{A \code{\link[OpenMx]{MxMatrix-class}} object 
}
\references{
}
\author{Mike W.-L. Cheung <mikewlcheung@nus.edu.sg>
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# a and b are free parameters with both starting values and labels
a1 <- matrix(c(1:4, "5*a", 6, "7*b", 8, 9), ncol=3, nrow=3)
a1 <- as.mxMatrix(a1, name="a1")

# fixed parameters without any labels
a2 <- matrix(1:9, ncol=3, nrow=3)
a2 <- as.mxMatrix(a2, name="a2")

# free parameters without starting values
a3 <- matrix(c(1:4, "a", 6, "b", 8, 9), ncol=3, nrow=3)
a3 <- as.mxMatrix(a3, name="a3")

# a free parameter without label
a4 <- matrix(c(1:4, "5*", 6, "7*b", 8, 9), ncol=3, nrow=3)
a4 <- as.mxMatrix(a4, name="a4")
}
\keyword{utilities}
