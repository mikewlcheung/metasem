\name{osmasem3}
\alias{osmasem3}
\title{Three-level one-stage meta-analytic structural equation modeling
}
\description{It fits MASEM with the one-stage MASEM with a three-level (OSMASEM3) approach.  
}
\usage{
osmasem3(model.name="osmasem3", RAM=NULL, cluster=NULL,
                     Mmatrix=NULL, TmatrixB=NULL, TmatrixW=NULL,
                     Ax=NULL, Sx=NULL, A.lbound=NULL, A.ubound=NULL,
                     RE.type=c("Diag", "Symm"), data,
                     subset.variables=NULL, subset.rows=NULL, 
                     intervals.type = c("z", "LB"),
                     mxModel.Args=NULL, mxRun.Args=NULL,
                     suppressWarnings=TRUE, silent=TRUE, run=TRUE, ...) 
}
\arguments{
\item{model.name}{A string for the model name in \code{\link[OpenMx]{mxModel}}. 
  }
\item{RAM}{A RAM object including a list of matrices of the model
  returned from \code{\link[metaSEM]{lavaan2RAM}}. If it is given,
  \code{Mmatrix}, \code{TmatrixB} and \code{TmatrixW} arguments will be ignored. 
}
\item{cluster}{A vector of string or number indicating the clusters.
}
\item{Mmatrix}{A list of matrices of the model implied correlation
  matrix created by the \code{create.vechsR}. It is only required when
  \code{RAM} is null.}
\item{TmatrixB}{A list of matrices of the between-level heterogeneity
  variance-covariance matrix created by the \code{create.Tau2}. It is
  only required when \code{RAM} is null.}
\item{TmatrixW}{A list of matrices of the within-level heterogeneity
  variance-covariance matrix created by the \code{create.Tau2}. It is
  only required when \code{RAM} is null.}
\item{Ax}{A Amatrix of a list of Amatrix with definition variables as
  the moderators of the Amatrix. It is used to create the \code{Mmatrix}.}
\item{Sx}{A Smatrix of a list of Smatrix with definition variables as
  the moderators of the Smatrix. It is used to create the
  \code{Mmatrix}.}
\item{A.lbound}{A matrix of lower bound of the Amatrix. If a scalar is
  given, the lbound matrix will be filled with this scalar.}
\item{A.ubound}{A matrix of upper bound of the Amatrix. If a scalar is
  given, the ubound matrix will be filled with this scalar.}
\item{RE.type}{Type of the random effects.}
\item{data}{A list of data created by the \code{Cor2DataFrame}.}
\item{subset.variables}{A character vector of the observed variables selected for
  the analysis.}
\item{subset.rows}{A logical vector of the same length as the number of rows in the data to select which rows are used in the analysis.}  
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
  Cheung, M. W.-L. (2014). Modeling dependent effect sizes with three-level meta-analyses: A structural equation modeling approach. \emph{Psychological Methods}, \bold{19}, 211-229.

  Jak, S., & Cheung, M. W.-L. (2020). Meta-analytic structural equation
modeling with moderating effects on SEM parameters. \emph{Psychological
  Methods}, \bold{25} (4), 430-455. https://doi.org/10.1037/met0000245
}  
\seealso{ \code{\link[metaSEM]{Cor2DataFrame}}, \code{\link[metaSEM]{create.vechsR}},
  \code{\link[metaSEM]{create.Tau2}}, \code{\link[metaSEM]{create.V}}, \code{\link[metaSEM]{osmasem}}, \code{\link[metaSEM]{Nohe15}}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ osmasem }
