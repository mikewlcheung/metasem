\name{osmasemSRMR}
\alias{osmasemSRMR}
\title{Calculate the SRMR in OSMASEM}
\description{It calculates the standardized root mean squared residuals (SRMR) in OSMASEM.  
}
\usage{
osmasemSRMR(x)
}
\arguments{
  \item{x}{An OSMASEM without any moderators.}
}

\value{It extracts the model implied correlation matrix and its
saturated counterpart to calculate the SRMR. It should be noted that
the heterogeneity variances are ignored in the calculations.}

\author{Mike W.-L. Cheung <mikewlcheung@nus.edu.sg>
}
\seealso{ \code{\link[metaSEM]{osmasem}}
}

\keyword{ osmasem }
