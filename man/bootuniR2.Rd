\name{bootuniR2}
\alias{bootuniR2}
\title{Fit Models on the bootstrapped correlation matrices
}
\description{It fits structural equation models on the bootstrapped
  correlation matrices. 
}
\usage{
bootuniR2(model, data, n, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{A model in \code{\link[lavaan]{sem}} syntax.}  
  \item{data}{A list of correlation matrices.}
  \item{n}{Sample size in fitting the structural equation models}  
  \item{\dots}{Further arguments to be passed to \code{\link[lavaan]{sem}}.
  }
}

\value{A list of the fitted object from \code{\link[lavaan]{sem}}.
}
\details{This function fits the lavaan model with the bootstrapped correlation
  matrices. It implements the parametric bootstrap approach suggested
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
\seealso{ \code{\link[metaSEM]{bootuniR2}},
  \code{\link[metaSEM]{tssemParaVar}}, \code{\link[metaSEM]{Nohe15}}
}
\keyword{bootuniR}
