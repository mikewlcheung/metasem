\name{coef}
\alias{coef.tssem1}
\alias{coef.tssem1.cluster}
\alias{coef.wls}
\alias{coef.wls.cluster}
\alias{coef.meta}
\alias{coef.reml}

\title{Extract Parameter Estimates from tssem1, tssem1.cluster, wls,
  wls.cluster, meta and reml Objects
}
\description{It extracts the parameter estimates from either
  \code{tssem1}, \code{tssem1.cluster}, \code{wls}, \code{wls.cluster},
  \code{meta} or \code{reml} objects.
}
\usage{
\method{coef}{tssem1}(object, \dots)
\method{coef}{tssem1.cluster}(object, \dots)
\method{coef}{wls}(object, \dots)
\method{coef}{wls.cluster}(object, \dots)
\method{coef}{meta}(object, \dots)
\method{coef}{reml}(object, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{An object returned from either class
	\code{\link[metaSEM]{tssem1}}, class
	\code{tssem1.cluster},
	class \code{\link[metaSEM]{wls}}, class \code{wls.cluster},
	class \code{\link[metaSEM]{meta}} or class \code{\link[metaSEM]{reml}}
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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }

