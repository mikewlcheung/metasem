\name{coef}
\alias{coef.tssem1FEM}
\alias{coef.tssem1FEM.cluster}
\alias{coef.tssem1REM}
\alias{coef.wls}
\alias{coef.wls.cluster}
\alias{coef.meta}
\alias{coef.reml}

\title{Extract Parameter Estimates from tssem1FEM, tssem1FEM.cluster,
  tssem1REM, wls, wls.cluster, meta and reml Objects
}
\description{It extracts the parameter estimates from either
  \code{tssem1FEM}, \code{tssem1FEM.cluster}, \code{tssem1REM}, \code{wls}, \code{wls.cluster},
  \code{meta} or \code{reml} objects.
}
\usage{
\method{coef}{tssem1FEM}(object, \dots)
\method{coef}{tssem1FEM.cluster}(object, \dots)
\method{coef}{tssem1REM}(object, select = c("all", "fixed", "random"), \dots)
\method{coef}{wls}(object, \dots)
\method{coef}{wls.cluster}(object, \dots)
\method{coef}{meta}(object, select = c("all", "fixed", "random"), \dots)
\method{coef}{reml}(object, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{An object returned from either class
	\code{tssem1FEM}, class \code{tssem1FEM.cluster}, class \code{tssem1REM}, 
	class \code{wls}, class \code{wls.cluster}, class \code{meta} or class \code{reml}
  }
  \item{select}{Select "all" for all parameters, "fixed" for the
  fixed-effects parameters or "random" for the random-effects
  parameters. 
}  
  \item{\dots}{Further arguments; currently none is used}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{Parameter estimates for both fixed-effects (if any) and random-effects (if any)
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

\seealso{ \code{\link[metaSEM]{tssem1}}, \code{\link[metaSEM]{wls}},
  \code{\link[metaSEM]{meta}}, \code{\link[metaSEM]{reml}}
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

