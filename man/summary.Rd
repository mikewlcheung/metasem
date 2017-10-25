\name{summary}
\alias{summary.tssem1FEM}
\alias{summary.tssem1FEM.cluster}
\alias{summary.tssem1REM}
\alias{summary.wls}
\alias{summary.wls.cluster}
\alias{summary.meta}
\alias{summary.meta3X}
\alias{summary.reml}
\alias{summary.CorPop}
\alias{summary.bootuniR2}
\alias{print.summary.tssem1FEM}
\alias{print.summary.tssem1FEM.cluster}
\alias{print.summary.wls}
\alias{print.summary.meta}
\alias{print.summary.meta3X}
\alias{print.summary.reml}
\alias{print.summary.CorPop}
\alias{print.summary.bootuniR2}
\title{Summary Method for tssem1, wls, meta and meta3X Objects
}
\description{It summaries results for either class \code{tssem1FEM},
  class \code{tssem1FEM.cluster}, class \code{tssem1REM}, class
  \code{wls}, class \code{wls.cluster}, class \code{meta}, class
  \code{meta3X}, class \code{reml}, class \code{CorPop}, or class \code{bootuniR2}.
}
\usage{
\method{summary}{tssem1FEM}(object, \dots)
\method{summary}{tssem1FEM.cluster}(object, \dots)
\method{summary}{tssem1REM}(object, \dots)
\method{summary}{wls}(object, df.adjustment=0, R=50, \dots)
\method{summary}{wls.cluster}(object, df.adjustment=0, R=50, \dots)
\method{summary}{meta}(object, homoStat=TRUE, \dots)
\method{summary}{meta3X}(object, allX=FALSE, \dots)
\method{summary}{reml}(object, \dots)
\method{summary}{CorPop}(object, \dots)
\method{summary}{bootuniR2}(object, probs=c(0, 0.1, 0.5, 0.9, 1),
cutoff.chisq.pvalue=0.05, cutoff.CFI=0.9, cutoff.SRMR=0.1, cutoff.RMSEA=0.05, \dots)
\method{print.summary}{tssem1FEM}(x, \dots)
\method{print.summary}{wls}(x, \dots)
\method{print.summary}{meta}(x, \dots)
\method{print.summary}{meta3X}(x, \dots)
\method{print.summary}{reml}(x, \dots)
\method{print.summary}{CorPop}(x, \dots)
\method{print.summary}{bootuniR2}(x, \dots)
}
\arguments{
  \item{object}{An object returned from either class
	\code{tssem1FEM}, class \code{tssem1FEM.cluster},
	class \code{tssem1REM}, class \code{wls}, class \code{wls.cluster},
	class \code{meta}, class \code{meta3X}, class \code{reml} or class \code{CorPop}.
}
  \item{x}{An object returned from either class \code{summary.tssem1FEM},
	class \code{tssem1FEM.cluster}, class \code{summary.wls}, class
	\code{summary.meta}, class \code{summary.meta3X}, class
	\code{summary.reml} or class \code{summary.CorPop}.
}
  \item{homoStat}{Logical. Whether to conduct a homogeneity test on the
	effect sizes.}
  \item{allX}{Logical. Whether to report the predictors and the
	auxiliary variables.}
  \item{df.adjustment}{Numeric. Adjust the degrees of freedom
	manually. It may be necessary if the df calculated is incorrect when
	\code{diag.constraints=TRUE}.
  }
  \item{R}{Positive integer. The number of parameter bootstrap
	replicates when \code{diag.constraints=TRUE}.}
  \item{probs}{Quantiles for the parameter estimates.}
  \item{cutoff.chisq.pvalue}{Cutoff of the p-value for the chi-square
	statistic.}
  \item{cutoff.CFI}{Cutoff of the CFI.}
  \item{cutoff.SRMR}{Cutoff of the SRMR.}
  \item{cutoff.RMSEA}{Cutoff of the RMSEA.}
  \item{\dots}{Further arguments to be passed to \code{\link[stats]{printCoefmat}}
}
}
\author{Mike W.-L. Cheung <mikewlcheung@nus.edu.sg>
%%  ~~who you are~~
}
\note{If the OpenMx status1 is either 0 or 1, the estimation is
  considered fine. If the OpenMx status1 is other values, it indicates
  estimation problems. Users should refer to \url{http://openmx.psyc.virginia.edu/wiki/errors} for more details.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link[metaSEM]{tssem1}}, \code{\link[metaSEM]{wls}},
  \code{\link[metaSEM]{meta}}, \code{\link[metaSEM]{reml}}, \code{\link[metaSEM]{rCor}}, \code{\link[metaSEM]{bootuniR2}}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }

