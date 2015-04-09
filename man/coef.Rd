\name{coef}
\alias{coef.tssem1FEM}
\alias{coef.tssem1FEM.cluster}
\alias{coef.tssem1REM}
\alias{coef.wls}
\alias{coef.wls.cluster}
\alias{coef.meta}
\alias{coef.meta3X}
\alias{coef.reml}
\alias{coef.MxRAMModel}

\title{Extract Parameter Estimates from tssem1FEM, tssem1FEM.cluster,
  tssem1REM, wls, wls.cluster, meta, meta3X, reml and MxRAMModel Objects
}
\description{It extracts the parameter estimates from either
  \code{tssem1FEM}, \code{tssem1FEM.cluster}, \code{tssem1REM}, \code{wls}, \code{wls.cluster},
  \code{meta}, \code{meta3X}, \code{reml} or \code{MxRAMModel} objects.
}
\usage{
\method{coef}{tssem1FEM}(object, \dots)
\method{coef}{tssem1FEM.cluster}(object, \dots)
\method{coef}{tssem1REM}(object, select = c("all", "fixed", "random"), \dots)
\method{coef}{wls}(object, \dots)
\method{coef}{wls.cluster}(object, \dots)
\method{coef}{meta}(object, select = c("all", "fixed", "random"), \dots)
\method{coef}{meta3X}(object, select = c("all", "fixed", "random", "allX"), \dots)
\method{coef}{reml}(object, \dots)
\method{coef}{MxRAMModel}(object, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{An object returned from either class
	\code{tssem1FEM}, class \code{tssem1FEM.cluster}, class \code{tssem1REM}, 
	class \code{wls}, class \code{wls.cluster}, class \code{meta}, class
	\code{reml} or class \code{MxRAMModel}
  }
  \item{select}{Select \code{all} for both fixed- and random-effects parameters, \code{fixed} for the
  fixed-effects parameters or \code{random} for the random-effects
  parameters. For \code{meta3X} objects, \code{allX} is used to extract
  all parameters including the predictors and auxiliary variables.  
}  
  \item{\dots}{Further arguments; currently none is used}
}
\note{\code{coef.MxRAMModel} is simply a wraper of
  \code{omxGetParameters}. Extra arguments will be passed to it}

\value{Parameter estimates for both fixed-effects (if any) and random-effects (if any)
}

\author{Mike W.-L. Cheung <mikewlcheung@nus.edu.sg>
}

\seealso{ \code{\link[metaSEM]{tssem1}}, \code{\link[metaSEM]{wls}},
  \code{\link[metaSEM]{meta}}, \code{\link[metaSEM]{reml}}, \code{\link[OpenMx]{omxGetParameters}}
}
\examples{
## Random-effects meta-analysis
model1 <- meta(y=yi, v=vi, data=Hox02)
coef(model1)

## Fixed-effects only
coef(model1, select="fixed")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }

