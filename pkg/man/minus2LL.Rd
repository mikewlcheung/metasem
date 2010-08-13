\name{minus2LL}
\alias{minus2LL}
\title{Calculate -2LogLikelihood of an independent or saturated model
}
\description{It calculates the -2LogLikelihood of a saturated or
  independent model on a correlation/covariance matrix. If \code{x} is a
  list of correlation/covariance matrices and \code{n} is a vector of sample sizes, a list of
  -2LogLikehood are calculated. 
}
\usage{
minus2LL(x, n, model = c("saturated", "independent"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A correlation/covariance matrix or a list of correlation/covariance matrices
}
  \item{n}{A sample size or a vector of sample sizes
}
   \item{model}{Either \code{independent} (default if missing) or \code{saturated}
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{-2LogLikelihood of the independent or saturated model or a list
  of it.
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

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# One study
minus2LL(Digman97[[1]], Digman97.n[1], model="independent")

# 14 independent studies
sum(minus2LL(Digman97, Digman97.n, model="saturated"))
}
\keyword{utilities}
