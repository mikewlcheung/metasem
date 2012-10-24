\name{indirectEffect}
\alias{indirectEffect}
\title{Estimate the asymptotic covariance matrix of standardized or unstandardized indirect and direct effects
}
\description{It estimates the standardized or unstandardized indirect and direct effects
  and their asymptotic sampling covariance matrix.
}
\usage{
indirectEffect(x, n, standardized = TRUE, direct.effect = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A 3x3 correlation/covariance matrix or a list of
	correlation/covariance matrices. Variables are
	arranged as dependent variable (y), mediator (m) and independent
	variable (x)
}
  \item{n}{Sample size or a vector of sample sizes
  }
  \item{standardized}{Logical. Whether the indirect effect is
  standardized.
}
  \item{direct.effect}{Logical. Whether the direct effect is
  estimated. If it is \code{FALSE}, the direct effect is fixed at zero.
  }
}
\details{Cheung (2009) estimated the standardized indirect effect and
  its standard error with non-linear constraints. Since \code{OpenMx} does not generate standard errors when there
  are non-linear constraints, Kwan and Chan's (2011) approach is used in
  this function. Delta method is used to calculate the asymptotic covariance matrix.
}
\value{A vector (or a matrix if the input is a list of matrices) of
  (standardized) indirect effect, standardized direct effect and their
  asymptotic sampling covariance matrices
}
\references{
  Cheung, M. W. L. (2009). Comparison of methods for constructing confidence intervals of standardized indirect effects. \emph{Behavior
  Research Methods}, \emph{41}, 425-438.
  
  Kwan, J., & Chan, W. (2011). Comparing standardized coefficients in
  structural equation modeling: a model reparameterization
  approach. \emph{Behavior Research Methods}, \emph{43}, 730-745.
}
\author{Mike W.-L. Cheung <mikewlcheung@nus.edu.sg>

}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## A correlation matrix as input
x <- matrix(c(1, 0.4, 0.2, 0.4, 1, 0.3, 0.2, 0.3, 1), ncol=3)
dimnames(x) <- list( c("y", "m", "x"), c("y", "m", "x") )
indirectEffect(x, n=300)

## A list of correlation matrices
indirectEffect( list(x, x), n=c(300,500), standardized=FALSE )
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ compute effect sizes }

