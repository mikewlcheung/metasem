\name{tssem1}
\alias{tssem1}
\alias{tssem1FE}
\alias{tssem1RE}
\title{First Stage of the Two-Stage Structural Equation Modeling (TSSEM)
}
\description{It conducts the first stage analysis of TSSEM by pooling
  correlation/covariance matrices. \code{tssem1FE()} and
  \code{tssem1RE()} use fixed- and random-effects models,
  respectively. \code{tssem1()} is a wrapper of these functions. 
}
\usage{
tssem1(my.df, n, method=c("FE","RE"), cor.analysis = TRUE, cluster=NULL,
       RE_diag=FALSE, RE.startvalues=0.1, RE.lbound = 1e-10,
       model.name=NULL, suppressWarnings=TRUE, ...)
tssem1FE(my.df, n, cor.analysis = TRUE, model.name = NULL,
       cluster = NULL, suppressWarnings = TRUE, ...)
tssem1RE(my.df, n, cor.analysis = TRUE, RE_diag = FALSE,
         RE.startvalues = 0.1, RE.lbound = 1e-10, model.name = NULL,
         suppressWarnings = TRUE, ...)
}
\arguments{
  \item{my.df}{A list of correlation/covariance matrices
}
  \item{n}{A vector of sample sizes
}
\item{method}{Either \code{"FE"} (default if missing) or \code{"RE"}.
  If it is "FE", fixed-effects meta-analysis will be applied. If it is "FE",
   random-effects meta-analysis will be applied.
}
  \item{cor.analysis}{Logical. The output is either a pooled correlation or a covariance matrix.
  }
  \item{cluster}{A vector of characters or numbers indicating the
	clusters. Analyses will be conducted for each cluster. It will be
	ignored when \code{method="RE"}.
  }
  \item{RE_diag}{Logical. If it is \code{TRUE} (default if missing),
	there are only random effects on the specific correlation/covariance elements
	while the correlation/covariance elments are assumed independent. If
	it is \code{FALSE}, there are random effects on the covariances
	among the correlation/covariance elements. It will be ignored when \code{method="FE"}.}
  \item{RE.startvalues}{Starting values on the
	diagonals of the variance component of the random effects. It will be ignored when \code{method="FE"}.}
  \item{RE.lbound}{Lower bounds on the diagonals of the variance
	component of the random effects. It will be ignored when \code{method="FE"}.}
  \item{model.name}{A string for the model name in \code{\link[OpenMx]{mxModel}}. 
  }  
  \item{suppressWarnings}{Logical. If \code{TRUE}, warnings are
	suppressed. Argument to be passed to \code{\link[OpenMx]{mxRun}}.}
  \item{\dots}{Further arguments to be passed to \code{\link[OpenMx]{mxRun}}
}
}
\details{
}
\value{Either an object of class \code{tssem1FE} for fixed-effects TSSEM,
  an object of class \code{tssem1FE.cluster} for fixed-effects TSSEM
  with \code{cluster} argument, or an object of class \code{tssem1RE}
  for random-effects TSSEM. 
}
\references{
  Cheung, M. W. L. (2011). Applications of meta-analytic structural
  equation modeling. \emph{Manuscript in preparation.}

  Cheung, M. W. L. (2011). Multivariate meta-analaysis as structural equation models. \emph{Manuscript submitted for publication.}

  Cheung, M. W. L., & Chan, W. (2005). Meta-analytic structural equation modeling: A two-stage approach. \emph{Psychological Methods}, \bold{10}, 40-64.

Cheung, M. W. L., & Chan, W. (2009). A two-stage approach to synthesizing covariance matrices in meta-analytic structural equation modeling. \emph{Structural Equation Modeling}, \bold{16}, 28-53.
}
\author{Mike W.-L. Cheung <mikewlcheung@nus.edu.sg>
}
\note{It assumes that there is no missing values in the first
  study. If the \code{cluster} argument is used, it returns a list of
  results based for each cluster.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link[metaSEM]{wls}}
}
\examples{
## Fixed-effects TSSEM
digmanFixed1 <- tssem1(Digman97$data, Digman97$n)
summary(digmanFixed1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ tssem }
