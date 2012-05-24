\name{create.mxMatrix}
\alias{create.mxMatrix}
\title{Create a Vector into MxMatrix-class
}
\description{It converts a vector into \code{MxMatrix-class} via \code{mxMatrix}.
}
\usage{
create.mxMatrix(x, type=c("Full","Symm","Diag","Stand"), ncol=NA,
nrow=NA, as.mxMatrix=TRUE, ...) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A character or numeric vector 
  }
  \item{type}{Matrix type similar to those listed in
	\code{\link[OpenMx]{mxMatrix}}
}
\item{ncol}{Number of columns. It is necessary when
  \code{type="Full"}. It is ignored and determined by the length of
  \code{x} for the other types of matrices.
}
\item{nrow}{Number of rows. It is necessary when
  \code{type="Full"}. It is ignored and determined by the length of
  \code{x} for the other types of matrices.
}
\item{as.mxMatrix}{Logical. If it is \code{TRUE}, the output is a matrix of
  \code{MxMatrix-class}. If it is \code{FALSE}, it is a numeric matrix.}
  \item{\dots}{Further arguments to be passed to
	\code{\link[OpenMx]{mxMatrix}}. Please note that \code{type},
	  \code{nrow}, \code{ncol}, \code{values}, \code{free} and \code{labels} will be created automatically. Thus, these arguments
	  excepts labels should be avoided in \dots 
%%     ~~Describe \code{\dots} here~~
}
}
\details{If there are non-numeric values in \code{x}, they are treated
  as the labels of the free parameters. If an "*" is
  present, the numeric value on the left hand side will be treated as
  the starting value for a free parameter or a fixed value for a fixed parameter. If it is a matrix of numeric values, there is
  no free parameters in the output matrix. \code{nrow} and \code{ncol}
  will be calculated from the length of \code{x} unless
  \code{type="Full"} is specified.
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
	\code{\link[metaSEM]{create.mxMatrix}},  \code{\link[metaSEM]{create.Fmatrix}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## a and b are free parameters with starting values and labels
(a1 <- c(1:4, "5*a", 6, "7*b", 8, 9))

mat1 <- create.mxMatrix(a1, ncol=3, nrow=3, name="mat1")

(a2 <- c(1:3, "4*f4", "5*f5", "6*f6"))

mat2 <- create.mxMatrix(a2, type="Symm", name="mat2")

a3 <- create.mxMatrix(a2, type="Symm", as.mxMatrix=FALSE)

mat4 <- create.mxMatrix(a2, type="Diag", lbound=1:6, name="mat4")
}
\keyword{utilities}
