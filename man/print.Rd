\name{print}
\alias{print.tssem1FEM}
\alias{print.tssem1FEM.cluster}
\alias{print.tssem1REM}
\alias{print.wls}
\alias{print.meta}
\alias{print.meta3LX}
\alias{print.reml}
\alias{print.uniR1}
\alias{print.impliedR}
\title{Print Methods for various Objects 
}
\description{Print methods for the \code{tssem1FEM},
  \code{tssem1FEM.cluster}, \code{tssem1REM}, \code{wls},
  \code{meta}, \code{meta3LX}, \code{reml}, \code{uniR1} and
  \code{impliedR} objects.
}
\usage{
\method{print}{tssem1FEM}(x, \dots)
\method{print}{tssem1FEM.cluster}(x, \dots)
\method{print}{tssem1REM}(x, \dots)
\method{print}{wls}(x, \dots)
\method{print}{meta}(x, \dots)
\method{print}{meta3LX}(x, \dots)
\method{print}{reml}(x, \dots)
\method{print}{uniR1}(x, \dots)
\method{print}{impliedR}(x, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{An object returned from either class
	\code{tssem1FEM}, class \code{tssem1FEM.cluster}, class \code{tssem1REM},
	class \code{wls}, class \code{meta}, class \code{meta3LX},
	class \code{reml}, class \code{uniR1} or class \code{impliedR}
}
  \item{\dots}{Further arguments to be passed to \code{summary.default}
	or unused.
}
}

\author{Mike W.-L. Cheung <mikewlcheung@nus.edu.sg>
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link[metaSEM]{tssem1}}, \code{\link[metaSEM]{wls}},
  \code{\link[metaSEM]{meta}},  \code{\link[metaSEM]{reml}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}

\keyword{ methods }
