\name{bootMASEM2}
\alias{bootMASEM2}
\title{Fit Models on the bootstrapped objects of the univariate R (uniR) approach
}
\description{It fits structural equation models on the boostrapped
  objects on the univariate R (uniR). 
}
\usage{
bootMASEM2(model, data, n, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{Futher arguments to be passed to either
	\code{\link[OpenMx]{mxRun}} or \code{\link[lavaan]{sem}}}  
  \item{data}{A list of correlation matrices.}
  \item{n}{Sample size in fitting the structural equation models}  
  \item{\dots}{Futher arguments to be passed to \code{\link[lavaan]{sem}}.
  }
}

\value{A list of fitted object from \code{\link[lavaan]{sem}}.
}
\details{It fitts the lavaan model with the bootstrapped correlation matrices.
}
\references{
  %% Cheung, M. W.-L. (2017). Issues in solving the problem of effect size
  %% heterogeneity in meta-analytic structural equation modeling:
  %% A commentary and simulation study on Yu, Downes, Carter, and O'Boyle (2016). 
  %% \emph{Manuscript under review}.

  Yu, J. (Joya), Downes, P. E., Carter, K. M., & O'Boyle,
  E. H. (2016). The problem of effect size heterogeneity 
  in meta-analytic structural equation modeling.
  \emph{Journal of Applied Psychology}, \emph{101}, 1457-1473.   
}
\author{Mike W.-L. Cheung <mikewlcheung@nus.edu.sg>
}
\seealso{ \code{\link[metaSEM]{bootMASEM2}}, \code{\link[metaSEM]{Nohe15}}
}
\keyword{bootMASEM}
