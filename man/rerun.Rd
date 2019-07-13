\name{rerun}
\alias{rerun}
\title{Rerun models via mxTryHard()
}
\description{It reruns models via mxTryHard().
}
\usage{
rerun(object, autofixtau2=FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{An object of either class \code{tssem1FEM},
  class \code{tssem1REM}, class \code{wls}, class \code{meta}, class
  \code{reml}, class \code{osmasem}, and class \code{MxModel}.
}
  \item{autofixtau2}{Logical. Whether automatically fixes elements of tau2
  with NA of standard errors. It only works for objects of
  class \code{tssem1REM}, class \code{meta}, and class \code{osmasem}.}
  \item{\dots}{{Further arguments to be passed to \code{\link[OpenMx]{mxTryHard}}}}
}

\author{Mike W.-L. Cheung <mikewlcheung@nus.edu.sg>
}

\examples{
\dontrun{
random1 <- tssem1(Digman97$data, Digman97$n, method="REM", RE.type="Diag")
random1_rerun <- rerun(random1)
summary(random1_rerun)
}
}
\keyword{tssem}
