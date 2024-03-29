\name{bootuniR1}
\alias{bootuniR1}
\title{Parametric bootstrap on the univariate R (uniR) object
}
\description{It generates correlation matrices with the parametric bootstrap on the univariate R
  (uniR) object. 
}
\usage{
bootuniR1(x, Rep, nonPD.pop=c("replace", "nearPD", "accept"))
}
\arguments{
  \item{x}{An object of class 'uniR1'}
  \item{Rep}{Number of replications of the parametric bootstrap}
  \item{nonPD.pop}{If it is \code{replace}, generated non-positive
	definite matrices are replaced by generated new ones which are
	positive definite. If it is \code{nearPD}, they are replaced by
	nearly positive definite matrices by calling
	\code{Matrix::nearPD()}. If it is \code{accept}, they are accepted.}
}
\value{An object of the generated correlation matrices.
}
\details{
  This function implements the parametric bootstrap approach suggested
  by Yu et al. (2016). It is included in this package for research
  interests. Please refer to Cheung (2018) for the issues associated
  with this parametric bootstrap approach.
}
\references{
  Cheung, M. W.-L. (2018). Issues in solving the problem of effect size
  heterogeneity in meta-analytic structural equation modeling: A
commentary and simulation study on Yu, Downes, Carter, and O'Boyle
(2016). \emph{Journal of Applied Psychology}, \bold{103}, 787-803.

  Yu, J. (Joya), Downes, P. E., Carter, K. M., & O'Boyle,
  E. H. (2016). The problem of effect size heterogeneity 
  in meta-analytic structural equation modeling.
  \emph{Journal of Applied Psychology}, \bold{101}, 1457-1473. 
}
\author{Mike W.-L. Cheung <mikewlcheung@nus.edu.sg>
}
\seealso{ \code{\link[metaSEM]{rCor}}, \code{\link[metaSEM]{bootuniR2}},
    \code{\link[metaSEM]{Nohe15}}
}
\keyword{bootuniR}
