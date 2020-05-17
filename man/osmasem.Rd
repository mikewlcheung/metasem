\name{osmasem}
\alias{osmasem}
\title{One-stage meta-analytic structural equation modeling
}
\description{It fits MASEM with the one-stage MASEM (OSMASEM) approach.  
}
\usage{
osmasem(model.name="osmasem", RAM=NULL, Mmatrix=NULL,
        Tmatrix=NULL, Jmatrix=NULL, Ax=NULL, Sx=NULL,
        RE.type=c("Diag", "Symm"), data,
        subset=NULL, intervals.type = c("z", "LB"),
        mxModel.Args=NULL, mxRun.Args=NULL,
        suppressWarnings=TRUE, silent=TRUE, run=TRUE, ...) 
}
\arguments{
\item{model.name}{A string for the model name in \code{\link[OpenMx]{mxModel}}. 
  }
\item{RAM}{A RAM object including a list of matrices of the model
  returned from \code{\link[metaSEM]{lavaan2RAM}}. If it is given,
  \code{Mmatrix} and \code{Tmatrix} arguments will be ignored. 
}  
\item{Mmatrix}{A list of matrices of the model implied correlation
  matrix created by the \code{create.vechsR}. It is only required when
  \code{RAM} is null.}
\item{Tmatrix}{A list of matrices of the heterogeneity
  variance-covariance matrix created by the \code{create.Tau2}. It is only required when \code{RAM} is null.}
\item{Jmatrix}{The Jacobian matrix of the mean structure in
  mxMatrix. The covariance structure is Jmatrix \%&\% Tau2 + Vi. If it is
not givin, an identity matrix will be used. }
\item{Ax}{A Amatrix of a list of Amatrix with definition variables as
  the moderators of the Amatrix. It is used to create the \code{Mmatrix}.}
\item{Sx}{A Smatrix of a list of Smatrix with definition variables as
  the moderators of the Smatrix. It is used to create the \code{Mmatrix}.}
\item{RE.type}{Type of the random effects.}
\item{data}{A list of data created by the \code{Cor2DataFrame}.}
\item{subset}{A character vector of the observed variables selected for
  the analysis.}
\item{intervals.type}{Either \code{z} (default if missing) or
  \code{LB}. If it is \code{z}, it calculates the 95\% confidence
  intervals (CIs) based on the estimated standard error. If it
  is \code{LB}, it calculates the 95\% likelihood-based CIs on the parameter estimates.}
\item{mxModel.Args}{A list of arguments passed to \code{\link[OpenMx]{mxModel}}.}
\item{mxRun.Args}{A list of arguments passed to \code{\link[OpenMx]{mxRun}}.}
  \item{suppressWarnings}{Logical. If it is \code{TRUE}, warnings are
	suppressed. This argument is passed to \code{\link[OpenMx]{mxRun}}.}
  \item{silent}{Logical. An argument is passed to \code{\link[OpenMx]{mxRun}}}
  \item{run}{Logical. If \code{FALSE}, only return the mx model without running the analysis.}
  \item{\dots}{Not used yet.}
  
}

\value{An object of class \code{osmasem}}

\author{Mike W.-L. Cheung <mikewlcheung@nus.edu.sg>
}

\references{
Jak, S., & Cheung, M. W.-L. (2019). Meta-analytic structural equation
modeling with moderating effects on SEM parameters. \emph{Psychological Methods}.
}  
\seealso{ \code{\link[metaSEM]{Cor2DataFrame}}, \code{\link[metaSEM]{create.vechsR}},
  \code{\link[metaSEM]{create.Tau2}}, \code{\link[metaSEM]{create.V}}, \code{\link[metaSEM]{osmasem}}, \code{\link[metaSEM]{Nohe15}}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ osmasem }
