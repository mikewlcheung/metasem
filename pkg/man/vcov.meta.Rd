\name{vcov.meta}
\alias{vcov.meta}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Variance-Covariance Matrix of Parameter Estimates for meta Objects
}
\description{It returns the variance-covariance matrix of the
  parameter estimates (both fixed-effects and variance components) for meta objects.
}
\usage{
vcov.meta(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{An object returned from \code{\link[metaSEM]{meta}}
}
  \item{\dots}{Further arguments (not used)}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{A variance-covariance matrix of the fixed- and random-effects
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

\seealso{ \code{\link[metaSEM]{meta}}
}
\examples{
## Random-effects meta-analysis
model1 <- with( Hox02, meta(y=yi, v=vi) )
vcov(model1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ meta-analysis }

