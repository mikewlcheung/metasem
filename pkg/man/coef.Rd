\name{coef}
\alias{coef.tssem1}
\alias{coef.wls}
\alias{coef.meta}
\alias{coef.reml}
\title{Extract Parameter Estimates from tssem1, wls, meta and reml Objects
}
\description{It extracts the parameter estimates from either
  \code{tssem1}, \code{wls}, \code{meta} or \code{reml} objects.
}
\usage{
\method{coef}{tssem1}(object, \dots)
\method{coef}{wls}(object, \dots)
\method{coef}{meta}(object, \dots)
\method{coef}{reml}(object, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{An object returned from either class
	\code{\link[metaSEM]{tssem1}}, class \code{\link[metaSEM]{wls}}, class \code{\link[metaSEM]{meta}} or class \code{\link[metaSEM]{reml}}
}
  \item{\dots}{Further arguments; currently none is used)}
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
model1 <- with( Hox02, meta(y=yi, v=vi) )
coef(model1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }

