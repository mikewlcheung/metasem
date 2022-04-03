\name{osmasemSRMR}
\alias{osmasemSRMR}
\alias{osmasem3LSRMR}
\title{Calculate the SRMR in OSMASEM and OSMASEM3L}
\description{It calculates the standardized root mean squared residuals
  (SRMR) in OSMASEM and OSMASEM3L.  
}
\usage{
osmasemSRMR(x)
osmasem3LSRMR(x)
}
\arguments{
  \item{x}{An OSMASEM object without any moderators.}
}

\value{It calculates the model implied correlation matrix and its
saturated counterpart to calculate the SRMR. It should be noted that
the heterogeneity variances are ignored in the calculations.}

\author{Mike W.-L. Cheung <mikewlcheung@nus.edu.sg>
}
\seealso{ \code{\link[metaSEM]{osmasem}}, \code{\link[metaSEM]{osmasem3L}},
  \code{\link[metaSEM]{Nohe15}}
}

\keyword{osmasem}
\keyword{osmasem3L}
