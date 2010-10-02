\name{coef}
\alias{coef.tssem1}
\alias{coef.wls}
\alias{coef.meta}
\title{Extract Parameter Estimates from tssem1, wls and meta Objects
}
\description{It extracts the parameter estimates from either \code{tssem1}, \code{wls} or \code{meta} objects.
}
\usage{
\method{coef}{tssem1}(object, \dots)
\method{coef}{wls}(object, \dots)
\method{coef}{meta}(object, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{An object returned from either class
	\code{\link[metaSEM]{tssem1}}, class \code{\link[metaSEM]{wls}} or
	class \code{\link[metaSEM]{meta}}
}
  \item{\dots}{Further arguments; currently none is used)}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{Parameter estimates for both fixed- and random-effects (if any)
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

\seealso{ \code{\link[metaSEM]{tssem1}}, \code{\link[metaSEM]{wls}}, \code{\link[metaSEM]{meta}}
}
\examples{
## Random-effects meta-analysis
model1 <- with( Hox02, meta(y=yi, v=vi) )
coef(model1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ summary }

