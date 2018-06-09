\name{VarCorr}
\alias{VarCorr}
\title{Extract Variance-Covariance Matrix of the Random Effects
}
\description{It extracts the variance-covariance matrix of the
  random effects (variance component) from either the \code{meta} or
  \code{osmasem} objects.
}
\usage{
VarCorr(x, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{An object returned from either class
	\code{meta} or \code{osmasem}
}
  \item{\dots}{Further arguments; currently none is used}
}

\value{A variance-covariance matrix of the random effects.
}

\author{Mike W.-L. Cheung <mikewlcheung@nus.edu.sg>
}
\note{It is similar to \code{coef(object, select="random")} in tssem. The main
  difference is that \code{coef()} returns a vector while
  \code{VarCorr()} returns its correspondent matrix.
}
\seealso{ \code{\link[metaSEM]{coef}}, \code{\link[metaSEM]{vcov}} }
\examples{
## Multivariate meta-analysis on the log of the odds
## The conditional sampling covariance is 0
bcg <- meta(y=cbind(ln_Odd_V, ln_Odd_NV), data=BCG,
            v=cbind(v_ln_Odd_V, cov_V_NV, v_ln_Odd_NV))
VarCorr(bcg)
}

%% ~Make other sections like Warning with \section{Warning }{....} ~
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }

