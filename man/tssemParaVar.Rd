\name{tssemParaVar}
\alias{tssemParaVar}
\title{Estimate the heterogeneity (SD) of the parameter estimates of the
  TSSEM object
}
\description{It estimates the heterogeneity of the parameter estimates
  of the TSSEM objects using either the bootstrap or the delta methods.
}
\usage{
tssemParaVar(tssem1.obj, tssem2.obj, method=c("bootstrap", "delta"),
             interval=0.8, Rep=50, output=c("data.frame", "matrices"),
             nonPD.pop=c("replace", "nearPD", "accept"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tssem1.obj}{An object of class \code{tssem1REM} returned from \code{tssem1()}
}
  \item{tssem2.obj}{An object of class \code{wls} returned from
\code{tssem2()} or \code{wls()}
}
  \item{method}{If it is \code{bootstrap}, random correlation matrices
are sampled from the \code{tssem1.obj} by the parametric bootstrap. If
it is \code{delta}, the delta method is used to estimate the
heterogeneity of the parameter estimates.
}
  \item{interval}{The desired interval, e.g., .8 or .95.
}
  \item{Rep}{The number of parametric bootstrap. It is ignored when the
method is \code{delta}.
}
  \item{output}{Either a \code{data.frame} or \code{matrices} of the output.
}
  \item{nonPD.pop}{If it is \code{replace}, generated non-positive
	definite matrices are replaced by generated new ones which are
	positive definite. If it is \code{nearPD}, they are replaced by
	nearly positive definite matrices by calling
	\code{Matrix::nearPD()}. If it is \code{accept}, they are accepted.}
}

\value{Either a \code{data.frame} or \code{matrices} of the output.
}
\details{
  The bootstrap method is based on the discussion in Cheung (in press)
and Yu et al. (2016). The delta method is an alternative method to
obtain the heterogeneity.
}
\references{
   Cheung, M. W.-L. (in press). Issues in solving the problem of effect size
  heterogeneity in meta-analytic structural equation modeling:
  A commentary and simulation study on Yu, Downes, Carter, and O'Boyle (2016). 
  \emph{Journal of Applied Psychology}.

  Yu, J. (Joya), Downes, P. E., Carter, K. M., & O'Boyle,
  E. H. (2016). The problem of effect size heterogeneity 
  in meta-analytic structural equation modeling.
  \emph{Journal of Applied Psychology}, \emph{101}, 1457-1473.   
}
\author{Mike W.-L. Cheung <mikewlcheung@nus.edu.sg>
}
\seealso{ \code{\link[metaSEM]{bootuniR1}}, \code{\link[metaSEM]{bootuniR2}}, \code{\link[metaSEM]{Nohe15}}
}
\keyword{tssem}
