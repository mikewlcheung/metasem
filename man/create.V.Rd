\name{create.V}
\alias{create.V}
\title{Create a V-known matrix}
\description{It creates a V-known matrix of the sampling covariance
  matrix using definition variables.  
}
\usage{
create.V(x, type = c("Symm", "Diag", "Full"), as.mxMatrix = TRUE)
}
\arguments{
\item{x}{A character vector of variable names of the sampling covariance matrix.}
\item{type}{Either \code{"Symm"}, \code{"Diag"} or
  \code{"Full"}. Suppose the number of variables is \eqn{p}, the
  numbers of variable names for \code{"Symm"}, \code{"Diag"}, and
  \code{"Full"} are \eqn{p(p-1)/2 }{p(p-1)/2}, \eqn{p}{p}, and
  \eqn{p*p}{p*p}, respectively. The elements are arranged in column major.} 
\item{as.mxMatrix}{Logical. Whether to convert the output into \code{MxMatrix-class}.}
}

\value{A list of \code{MxMatrix-class}. The V-known sampling covariance
  matrix is computed in \code{V}.}

\author{Mike W.-L. Cheung <mikewlcheung@nus.edu.sg>
}
\seealso{ \code{\link[metaSEM]{osmasem}},
  \code{\link[metaSEM]{create.Tau2}}, \code{\link[metaSEM]{create.vechsR}}
}

\examples{
\dontrun{
my.df <- Cor2DataFrame(Nohe15A1$data, Nohe15A1$n)

## Create known sampling variance covariance matrix
V0 <- create.V(my.df$vlabels)
V0
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ osmasem }
