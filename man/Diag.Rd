\name{Diag}
\alias{Diag}
\alias{Diag<-}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Matrix Diagonals
}
\description{Extract or replace the diagonal of a matrix, or construct a diagonal matrix with the same behaviors as \code{diag} prior to R-3.0.0.
}
\usage{
Diag(x, ...)
Diag(x) <- value
}

\arguments{
  \item{x}{a matrix, vector or 1D array, or missing.}
  \item{...}{Optional dimensions (\code{nrow} and \code{ncol}) for the result when \code{x} is
    not a matrix.}
  \item{value}{either a single value or a vector of length equal to that
    of the current diagonal.  Should be of a mode which can be coerced
    to that of \code{x}.}
}

\details{Started from R-3.0.0, \code{diag(x)} returns a numeric matrix with NA in the
  diagonals when x is a character vector. Although this follows what the manual says, this
  breaks the metaSEM. The \code{Diag} has the same functions as
  \code{diag} except that \code{Diag(x)} works for a character vector of
  x by returning a square matrix of character "0" with \code{x} as the
  diagonals.
%%  ~~ If necessary, more details than the description above ~~
}

\author{Mike W.-L. Cheung <mikewlcheung@nus.edu.sg>
}
\note{See
  http://r.789695.n4.nabble.com/Behaviors-of-diag-with-character-vector-in-R-3-0-0-td4663735.html
  for the discussion.
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{diag}}
}
\examples{
v <- c("a", "b")
Diag(v)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}

