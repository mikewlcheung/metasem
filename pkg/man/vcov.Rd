\name{vcov}
\alias{vcov.tssem1}
\alias{vcov.tssem1.cluster}
\alias{vcov.wls}
\alias{vcov.wls.cluster}
\alias{vcov.meta}
\alias{vcov.reml}
\title{Extract Covariance Matrix Parameter Estimates from Various Objects
}
\description{It extracts the variance-covariance matrix of the
  parameter estimates from either \code{tssem1}, \code{tssem1.cluster},
  \code{wls}, \code{wls.cluster}, \code{meta}
  or \code{reml} objects.
}
\usage{
\method{vcov}{tssem1}(object, \dots)
\method{vcov}{tssem1.cluster}(object, \dots)
\method{vcov}{wls}(object, \dots)
\method{vcov}{wls.cluster}(object, \dots)
\method{vcov}{meta}(object, select = c("all", "fixed", "random"), \dots)
\method{vcov}{reml}(object, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{An object returned from either class
	\code{\link[metaSEM]{tssem1}}, class \code{tssem1.cluster},
	class \code{\link[metaSEM]{wls}}, class \code{wls.cluster}, class \code{\link[metaSEM]{meta}} or class \code{\link[metaSEM]{reml}}
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
\value{A variance-covariance matrix of the parameter estimates.
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
vcov(model1)

## Fixed-effects only
vcov(model1, select="fixed")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }

