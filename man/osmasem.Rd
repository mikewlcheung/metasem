\name{osmasem}
\alias{osmasem}
\title{One-stage meta-analytic structural equation modeling
}
\description{It fits MASEM with the one-stage MASEM (OSMASEM) approach.  
}
\usage{
osmasem(model.name="osmasem", Mmatrix, Tmatrix, data,
        intervals.type=c("z", "LB"), mxModel.Args=NULL,
        mxRun.Args=NULL, suppressWarnings=TRUE,
        silent=FALSE, run=TRUE, ...)
}
\arguments{
  \item{model.name}{A string for the model name in \code{\link[OpenMx]{mxModel}}. 
}
\item{Mmatrix}{A list of matrices of the model implied correlation
  matrix created by the \code{create.vechsR}.}
\item{Tmatrix}{A list of matrices of the heterogeneity
  variance-covariance matrix created by the \code{create.Tau2}.}
\item{data}{A list of data created by the \code{Cor2DataFrame}.}
\item{intervals.type}{Either \code{z} (default if missing) or
  \code{LB}. If it is \code{z}, it calculates the 95\% confidence
  intervals (CIs) based on the estimated standard error. If it
  is \code{LB}, it calculates the 95\% likelihood-based CIs on the parameter estimates.}
\item{mxModel.Args}{A list of arguments passed to \code{\link[OpenMx]{mxModel}}.}
\item{mxRun.Args}{A list of arguments passed to \code{\link[OpenMx]{mxRun}}.}
  \item{suppressWarnings}{Logical. If it is \code{TRUE}, warnings are
	suppressed. This argument is passed to \code{\link[OpenMx]{mxRun}}.}
  \item{silent}{Logical. The argument is passed to \code{\link[OpenMx]{mxRun}}}
  \item{run}{Logical. If \code{FALSE}, only return the mx model without running the analysis.}
  \item{\dots}{Not used yet.}
  
}

\value{An object of class \code{osmasem}}

\author{Mike W.-L. Cheung <mikewlcheung@nus.edu.sg>
}
\seealso{ \code{\link[metaSEM]{Cor2DataFrame}}, \code{\link[metaSEM]{create.vechsR}},
  \code{\link[metaSEM]{create.Tau2}}, \code{\link[metaSEM]{create.V}}, \code{\link[metaSEM]{Nohe15}}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ osmasem }
