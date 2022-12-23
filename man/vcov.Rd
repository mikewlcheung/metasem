\name{vcov}
\alias{vcov.tssem1FEM}
\alias{vcov.tssem1FEM.cluster}
\alias{vcov.tssem1REM}
\alias{vcov.wls}
\alias{vcov.wls.cluster}
\alias{vcov.meta}
\alias{vcov.meta3LFIML}
\alias{vcov.reml}
\alias{vcov.MxRAMModel}
\alias{vcov.osmasem}
\alias{vcov.osmasem3L}
\alias{vcov.tssemRobust1}
\title{Extract Covariance Matrix Parameter Estimates from Objects of
  Various Classes
}
\description{It extracts the variance-covariance matrix of the
  parameter estimates from objects of various classes.
}
\usage{
\method{vcov}{tssem1FEM}(object, \dots)
\method{vcov}{tssem1FEM.cluster}(object, \dots)
\method{vcov}{tssem1REM}(object, select = c("all", "fixed", "random"), robust=FALSE, \dots)
\method{vcov}{wls}(object, \dots)
\method{vcov}{wls.cluster}(object, \dots)
\method{vcov}{meta}(object, select = c("all", "fixed", "random"), robust=FALSE, \dots)
\method{vcov}{meta3LFIML}(object, select = c("all", "fixed", "random","allX"), robust=FALSE, \dots)
\method{vcov}{reml}(object, \dots)
\method{vcov}{MxRAMModel}(object, \dots)
\method{vcov}{osmasem}(object, select=c("fixed", "all", "random"), robust=FALSE, \dots)
\method{vcov}{osmasem3L}(object, select=c("fixed", "all", "random"), \dots)
\method{vcov}{tssemRobust1}(object, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{An object returned from objects of various classes
  }
  \item{select}{Select \code{all} for both fixed- and random-effects parameters, \code{fixed} for the
  fixed-effects parameters or \code{random} for the random-effects
  parameters. For \code{meta3LFIML} objects, \code{allX} is used to extract
  all parameters including the predictors and auxiliary variables.  
}
  \item{robust}{Logicial. Whether to use robust standard error from \code{\link[OpenMx]{imxRobustSE}}.
}
  \item{\dots}{Further arguments; currently not in use except for
	\code{tssemRobust1}, which to be passed to \code{\link[metafor]{robust}}.}
}

\value{A variance-covariance matrix of the parameter estimates.
}

\author{Mike W.-L. Cheung <mikewlcheung@nus.edu.sg>
}
\note{\code{vcov} returns \code{NA} when the \code{diag.constraints=TRUE}
  argument is used in \code{wls} objects.
}

\seealso{ \code{\link[metaSEM]{tssem1}}, \code{\link[metaSEM]{wls}},
  \code{\link[metaSEM]{meta}}, \code{\link[metaSEM]{reml}}
}
\examples{
## Random-effects meta-analysis
model1 <- meta(y=yi, v=vi, data=Hox02)
vcov(model1)

## Fixed-effects only
vcov(model1, select="fixed")

## Random-effects only
vcov(model1, select="random")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }

