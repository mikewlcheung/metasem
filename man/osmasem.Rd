\name{osmasem}
\alias{osmasem}
\alias{osmasem2}
\title{One-stage meta-analytic structural equation modeling
}
\description{It fits MASEM with the one-stage MASEM (OSMASEM) approach.  
}
\usage{
osmasem(model.name="osmasem", RAM=NULL, Mmatrix=NULL,
        Tmatrix=NULL, Jmatrix=NULL, Ax=NULL, Sx=NULL,
        A.lbound=NULL, A.ubound=NULL,
        RE.type=c("Diag", "Symm", "Zero"), data,
        subset.variables=NULL, subset.rows=NULL, 
        intervals.type = c("z", "LB"),
        mxModel.Args=NULL, mxRun.Args=NULL,
        suppressWarnings=TRUE, silent=TRUE, run=TRUE, ...)
osmasem2(model.name="osmasem2", RAM, data, cor.analysis=TRUE,
         RE.type.Sigma=c("Diag", "Symm", "Zero"),
         RE.type.Mu=c("Symm", "Diag", "Zero"),
         RE.type.SigmaMu=c("Zero", "Full"),
         mean.analysis=FALSE, startvalues=NULL,
         replace.constraints=FALSE,
         mxModel.Args=NULL, run=TRUE, ...)
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
\item{cor.analysis}{Whether to analyze correlation or covariance structure analysis.}
\item{RE.type.Sigma}{Type of the random effects of the correlation or covariance vectors.}
\item{RE.type.Mu}{Type of the random effects of the mean vectors.}
\item{RE.type.SigmaMu}{Type of the random effects between the correlation/covariance vectors and the mean vectors.}
\item{mean.analysis}{Whether to include the analysis of the mean structure.}
\item{startvalues}{An optional list of starting values. It is useful when there are new parameters in RAM.}
\item{replace.constraints}{It is relevant only when there are constraints in RAM. If it is \code{FALSE}, these constraints will be impose. If it is \code{FALSE}, the parameters on the left-hand side will be replaced by the algebras on the right-hand side.}
}
\details{\code{osmasem()} was implemented based on Jak and Cheung (2020) for meta-analyzing correlation matrices. \code{osmasem2()} was a rewritten version to handle correlation or covariance matrices with the means. There are several major differences between them:
1. \code{osmasem()} allows the use of \code{RAM} or (\code{Mmatrix} and \code{Tmatrix}), while \code{osmasem2()} only accepts the \code{RAM} input.
2. \code{RE.type} is used to specify the type of random effects on the correlations in \code{osmasem()}. On the contrary, \code{osmasem2()} has three types of random effects: correlations/covariances, means, and covariance between correlations/covariance and means.
3. \code{osmasem()} reports the transformed random effects in the parameter table. Users have to use \code{VarCorr()} to obtain the heterogeneity matrix of the random effects. In contrast, \code{osmasem2()} reports the heterogeneity matrix in the parameter table.
4. \code{osmasem()} accepts either \code{intervals.type="z"} or \code{intervals.type="LB"}, whereas \code{osmasem2()} only uses \code{intervals.type="z"}. Thus, this argument is removed in \code{osmasem2()}.
5. \code{osmasem2()} allows the imposing linear and nonlinear constraints and the creation of parameter functions in \code{RAM}, which \code{osmasem} cannot.
}
\value{An object of class \code{osmasem}}

\author{Mike W.-L. Cheung <mikewlcheung@nus.edu.sg>
}

\references{
  Jak, S., & Cheung, M. W.-L. (2020). Meta-analytic structural equation modeling with moderating effects on SEM parameters. \emph{Psychological
  Methods}, \bold{25} (4), 430-455. https://doi.org/10.1037/met0000245

  Jak, S., Cheung, M. W.-L., Acar, S., & Kindred, R. (December 18, 2024). Evaluating differences in latent means across studies: Extending meta-analytic confirmatory factor analysis with the analysis of means. OSF. https://doi.org/10.31234/osf.io/35gtz.

}  
\seealso{ \code{\link[metaSEM]{Cor2DataFrame}}, \code{\link[metaSEM]{create.vechsR}},
  \code{\link[metaSEM]{create.Tau2}}, \code{\link[metaSEM]{create.V}}, \code{\link[metaSEM]{osmasem}}, \code{\link[metaSEM]{Nohe15}}, \code{\link[metaSEM]{issp05}}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{meta-analysis}
