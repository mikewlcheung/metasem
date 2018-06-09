\name{osmasemR2}
\alias{osmasemR2}
\title{Calculate the R2 in OSMASEM}
\description{It calculate the R2 of the moderators in explaining the
  variances in the heterogeneity variances.  
}
\usage{
osmasemR2(object, R2.truncate=TRUE)
}
\arguments{
\item{object}{An object in class \code{osmasem}.}
\item{R2.truncate}{Whether to truncate the negative R2 to zero.}
}

\value{A list of the diagonals of the heterogeneity variances of the
  models without and with the moderators, and the R2.}

\author{Mike W.-L. Cheung <mikewlcheung@nus.edu.sg>
}
\seealso{ \code{\link[metaSEM]{osmasem}}
}

\keyword{ osmasem }
