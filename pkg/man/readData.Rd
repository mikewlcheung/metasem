\name{readData}
\alias{readFullMat}
\alias{readLowTriMat}
\alias{readStackVec}
\title{Read correlation/covariance data into a list of matrices
}
\description{It reads full/lower triangle/stacked vectors of correlation/covariance data into a list of correlation/covariance matrices.
}
\usage{
readFullMat(file, ...)
readLowTriMat(file, no.var, ...)
readStackVec(file, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file}{File name of the data.
}
  \item{no.var}{Number of variables in the data.
}
  \item{\dots}{Further arguments to be passed to \code{\link[base]{scan}} for \code{readLowTriMat} and to \code{\link[utils]{read.table}} for \code{readFullMat} and \code{readStackVec}.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{A list of correlation/covariance matrices.
}
\references{
%% ~put references to the literature/web site here ~
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
\dontrun{
# Write two full correlation matrices into a file named "fullmat.dat".
# x2 is missing in the second matrix.
cat("1.0 0.3 0.4\n0.3 1.0 0.5\n0.4 0.5 1.0\n1.0 NA 0.4\nNA NA NA\n0.4 NA 1.0", file="fullmat.dat", sep="")

# Read the correlation matrices
my.full <- readFullMat("fullmat.dat")

my.full
#$`1`
#    x1  x2  x3
#x1 1.0 0.3 0.4
#x2 0.3 1.0 0.5
#x3 0.4 0.5 1.0
#
#$`2`
#    x1 x2  x3
#x1 1.0 NA 0.4
#x2  NA NA  NA
#x3 0.4 NA 1.0

# Write two lower triangle correlation matrices into a file named "lowertriangle.dat".
# x2 is missing in the second matrix.
cat("1.0\n0.3 1.0\n0.4 0.5 1.0\n1.0\nNA NA\n0.4 NA 1.0", file="lowertriangle.dat", sep="")

# Read the lower triangle correlation matrices
my.lowertri <- readLowTriMat(file = "lowertriangle.dat", no.var = 3)

my.lowertri
#$`1`
#    x1  x2  x3
#x1 1.0 0.3 0.4
#x2 0.3 1.0 0.5
#x3 0.4 0.5 1.0
#
#$`2`
#    x1 x2  x3
#x1 1.0 NA 0.4
#x2  NA NA  NA
#x3 0.4 NA 1.0

# Write two vectors of correlation coefficients based on column major into a file named "stackvec.dat".
# x2 is missing in the second matrix.
cat("1.0 0.3 0.4 1.0 0.5 1.0\n1.0 NA 0.4 NA NA 1.0", file="stackvec.dat", sep="")

my.vec <- readStackVec("stackvec.dat")

my.vec
#$`1`
#    x1  x2  x3
#x1 1.0 0.3 0.4
#x2 0.3 1.0 0.5
#x3 0.4 0.5 1.0
#
#$`2`
#    x1 x2  x3
#x1 1.0 NA 0.4
#x2  NA NA  NA
#x3 0.4 NA 1.0

}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ utilities }
