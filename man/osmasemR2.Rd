\name{osmasemR2}
\alias{osmasemR2}
\alias{osmasem3R2}
\title{Calculate the R2 in OSMASEM and OSMASEM3}
\description{It calculates the R2 of the moderators in explaining the
  variances in the heterogeneity variances.  
}
\usage{
osmasemR2(model1, model0, R2.truncate=TRUE)
}
\arguments{
  \item{model1}{An object in class \code{osmasem}.}
  \item{model0}{An object in class \code{osmasem}.}
  \item{R2.truncate}{Whether to truncate the negative R2 to zero.}
}

\value{\code{model1} and \code{model0} are the models with and without
  the moderators, respectively. The function does not check whether the
  models are nested. It is the users' responsibility to make sure that
  the models with and without the moderators are nested. It returns a
  list of the diagonals of the heterogeneity variances of the
  models without and with the moderators, and the R2.}

\author{Mike W.-L. Cheung <mikewlcheung@nus.edu.sg>
}
\seealso{ \code{\link[metaSEM]{osmasem}, \link[metaSEM]{osmasem3}}
}

\keyword{ osmasem }
