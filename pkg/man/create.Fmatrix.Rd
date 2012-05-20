\name{create.Fmatrix}
\alias{create.Fmatrix}
\title{Create an F matrix to select observed variables
}
\description{It creates an F matrix to select observed variables for wls
  function.
}
\usage{
create.Fmatrix(x, name)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A vector of logical type 
}
\item{name}{Name of the matrix. If it is missing, "Fmatrix" will be used.
}
}
\details{
}
\references{
}
\author{Mike W.-L. Cheung <mikewlcheung@nus.edu.sg>
}
\note{
%%  ~~further notes~~
}
\seealso{ \code{\link[metaSEM]{as.mxMatrix}}, \code{\link[metaSEM]{create.mxMatrix}}, \code{\link[metaSEM]{wls}}
}

\examples{
## Select the first 3 variables while the other 2 variables are latent.
create.Fmatrix(c(1,1,1,0,0))
}
\keyword{utilities}
