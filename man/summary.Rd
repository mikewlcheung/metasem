\name{summary}
\alias{summary.tssem1FEM}
\alias{summary.tssem1FEM.cluster}
\alias{summary.tssem1REM}
\alias{summary.wls}
\alias{summary.wls.cluster}
\alias{summary.meta}
\alias{summary.meta3LFIML}
\alias{summary.reml}
\alias{summary.CorPop}
\alias{summary.Cor3L}
\alias{summary.bootuniR2}
\alias{summary.osmasem}
\alias{print.summary.tssem1FEM}
\alias{print.summary.tssem1FEM.cluster}
\alias{print.summary.wls}
\alias{print.summary.meta}
\alias{print.summary.meta3LFIML}
\alias{print.summary.reml}
\alias{print.summary.CorPop}
\alias{print.summary.Cor3L}
\alias{print.summary.bootuniR2}
\title{Summary Method for tssem1, wls, meta, and meta3LFIML Objects
}
\description{It summaries results for various class.
}
\usage{
\method{summary}{tssem1FEM}(object, \dots)
\method{summary}{tssem1FEM.cluster}(object, \dots)
\method{summary}{tssem1REM}(object, robust=FALSE, \dots)
\method{summary}{wls}(object, df.adjustment=0, \dots)
\method{summary}{wls.cluster}(object, df.adjustment=0, \dots)
\method{summary}{meta}(object, homoStat=TRUE, robust=FALSE, \dots)
\method{summary}{meta3LFIML}(object, allX=FALSE, robust=FALSE, \dots)
\method{summary}{reml}(object, \dots)
\method{summary}{CorPop}(object, \dots)
\method{summary}{Cor3L}(object, \dots)
\method{summary}{bootuniR2}(object, probs=c(0, 0.1, 0.5, 0.9, 1),
        cutoff.chisq.pvalue=0.05, cutoff.CFI=0.9, cutoff.SRMR=0.1,
        cutoff.RMSEA=0.05, \dots)
\method{summary}{osmasem}(object, fitIndices=FALSE, numObs, robust=FALSE, \dots)
\method{print.summary}{tssem1FEM}(x, \dots)
\method{print.summary}{wls}(x, \dots)
\method{print.summary}{meta}(x, \dots)
\method{print.summary}{meta3LFIML}(x, \dots)
\method{print.summary}{reml}(x, \dots)
\method{print.summary}{CorPop}(x, \dots)
\method{print.summary}{Cor3L}(x, \dots)
\method{print.summary}{bootuniR2}(x, \dots)
}
\arguments{
  \item{object}{An object returned from either class
	\code{tssem1FEM}, class \code{tssem1FEM.cluster},
	class \code{tssem1REM}, class \code{wls}, class \code{wls.cluster},
	class \code{meta}, class \code{meta3LFIML}, class \code{reml} or class \code{CorPop}.
}
  \item{x}{An object returned from either class \code{summary.tssem1FEM},
	class \code{tssem1FEM.cluster}, class \code{summary.wls}, class
	\code{summary.meta}, class \code{summary.meta3LFIML}, class
	\code{summary.reml} or class \code{summary.CorPop}.
}
  \item{homoStat}{Logical. Whether to conduct a homogeneity test on the
	effect sizes.}
  \item{allX}{Logical. Whether to report the predictors and the
	auxiliary variables.}
  \item{robust}{Logicial. Whether to use robust standard error from \code{\link[OpenMx]{imxRobustSE}}.
}
  \item{df.adjustment}{Numeric. Adjust the degrees of freedom
	manually. It may be necessary if the df calculated is incorrect when
	\code{diag.constraints=TRUE}.
  }
  \item{probs}{Quantiles for the parameter estimates.}
  \item{cutoff.chisq.pvalue}{Cutoff of the p-value for the chi-square
	statistic.}
  \item{cutoff.CFI}{The cutoff of the CFI.}
  \item{cutoff.SRMR}{The cutoff of the SRMR.}
  \item{cutoff.RMSEA}{The cutoff of the RMSEA.}
  \item{fitIndices}{Whether to calculate the chi-square statistic and
	various goodness-of-fit indices in osmasem. Note. It may take a
	while since statistics of the saturated and independence models are required.}
  \item{numObs}{The number of observations in calculating the fit statistics
	in osmasem. If it is missing, the total number of observations is used.}
  \item{\dots}{Further arguments to be passed to \code{\link[stats]{printCoefmat}}
}
}
\author{Mike W.-L. Cheung <mikewlcheung@nus.edu.sg>
%%  ~~who you are~~
}
\note{If the OpenMx status1 is either 0 or 1, the estimation is
  considered fine. If the OpenMx status1 is other values, it indicates
  estimation problems. Users should refer to \url{https://openmx.ssri.psu.edu//wiki/errors} for more details.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link[metaSEM]{tssem1}}, \code{\link[metaSEM]{wls}},
  \code{\link[metaSEM]{meta}}, \code{\link[metaSEM]{reml}},
  \code{\link[metaSEM]{rCor}}, \code{\link[metaSEM]{bootuniR2}},
  \code{\link[metaSEM]{osmasem}}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }

