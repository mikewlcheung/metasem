\name{standIndirect}
\alias{standIndirect}
\title{Asymptotic covariance matrix of standardized indirect and direct effects
}
\description{It calculates the asymptotic covariance matrix of standardized indirect
  and direct effects.
}
\usage{
standIndirect(x, n)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A 3x3 correlation/covariance matrix or a list of
	correlation/covariance matrices. Variables are
	arranged as dependent variable, mediator and independent variable
}
  \item{n}{Sample size or a vector of sample sizes
}
}
\details{Cheung (2009) uses non-linear constraints to estimate the
  standardized indirect effect and its standard error. Since \code{OpenMx} does not generate standard errors when there
  are non-linear constraints, Kwan and Chan's (in press) approach is used in this function. The calculated standard errors are the same for both approaches.
}
\value{A vector (or a matrix if the input is a list of matrices) of standardized indirect effect, standardized direct effect and their asymptotic
  sampling covariance matrices
}
\references{
  Cheung, M. W. L. (2009). Comparison of methods for constructing confidence intervals of standardized indirect effects. \emph{Behavior
  Research Methods}, \emph{41}, 425-438.
  
  Kwan, J., & Chan, W. (in press). Comparing standardized coefficients in structural equation modeling: a model reparameterization approach. \emph{Behavior Research Methods}.
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
my.df <- matrix(c(1, 0.4, 0.2, 0.4, 1, 0.3, 0.2, 0.3, 1), ncol=3)
dimnames(my.df) <- list( c("y", "m", "x"), c("y", "m", "x") )
standIndirect(my.df, n=300)

## A list of correlation matrices
standIndirect( list(my.df, my.df), n=c(300,500) )
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ utilities }

