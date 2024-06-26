\name{rerun}
\alias{rerun}
\title{Rerun models via mxTryHard()
}
\description{It reruns models via mxTryHard().
}
\usage{
rerun(object, autofixtau2=FALSE, extraTries=10, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{An object of either class \code{tssem1FEM},
  class \code{tssem1REM}, class \code{wls}, class \code{meta}, class
  \code{reml}, class \code{osmasem}, class \code{osmasem3L}, and class \code{MxModel}.
}
  \item{autofixtau2}{Logical. Whether automatically fixes elements of tau2
  with NA of standard errors. It only works for objects of
  class \code{tssem1REM}, class \code{meta}, and class \code{osmasem}.}
\item{extraTries}{The number of attempts to run the model in addition to the first.}  
  \item{\dots}{{Further arguments to be passed to \code{\link[OpenMx]{mxTryHard}}}}
}

\author{Mike W.-L. Cheung <mikewlcheung@nus.edu.sg>
}

\examples{
\donttest{
random1 <- tssem1(Digman97$data, Digman97$n, method="REM", RE.type="Diag")
random1_rerun <- rerun(random1)
summary(random1_rerun)
}
}
\keyword{tssem}
\keyword{meta}
\keyword{osmasem}
\keyword{osmasem3L}
\keyword{wls}
