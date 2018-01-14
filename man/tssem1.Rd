\name{tssem1}
\alias{tssem1}
\alias{tssem1FEM}
\alias{tssem1REM}
\title{First Stage of the Two-Stage Structural Equation Modeling (TSSEM)
}
\description{It conducts the first stage analysis of TSSEM by pooling
  correlation/covariance matrices. \code{tssem1FEM()} and
  \code{tssem1REM()} use fixed- and random-effects models,
  respectively. \code{tssem1()} is a wrapper of these functions. 
}
\usage{
tssem1(Cov, n, method=c("REM","FEM"), cor.analysis = TRUE, cluster=NULL,
       RE.type=c("Diag", "Symm", "Zero", "User"), RE.startvalues=0.1,
       RE.lbound=1e-10, RE.constraints=NULL, I2="I2q",
       acov=c("weighted", "individual", "unweighted"),
       model.name=NULL, suppressWarnings=TRUE, silent=TRUE, run=TRUE, ...)
tssem1FEM(Cov, n, cor.analysis=TRUE, model.name=NULL,
          cluster=NULL, suppressWarnings=TRUE, silent=TRUE, run=TRUE, ...)
tssem1REM(Cov, n, cor.analysis=TRUE, RE.type=c("Diag", "Symm", "Zero","User"),
          RE.startvalues=0.1, RE.lbound=1e-10, RE.constraints=NULL,
          I2="I2q", acov=c("weighted", "individual", "unweighted"),
          model.name=NULL, suppressWarnings=TRUE,
          silent=TRUE, run=TRUE, ...)
}
\arguments{
  \item{Cov}{A list of correlation/covariance matrices
}
  \item{n}{A vector of sample sizes
}
\item{method}{Either \code{"REM"} (default if missing) or \code{"FEM"}.
  If it is "REM",a random-effects meta-analysis will be applied. If it
  is "FEM", a fixed-effects meta-analysis will be applied. 
}
  \item{cor.analysis}{Logical. The output is either a pooled correlation or a covariance matrix.
  }
  \item{cluster}{A vector of characters or numbers indicating the
	clusters. Analyses will be conducted for each cluster. It will be
	ignored when \code{method="REM"}.
  }
  \item{RE.type}{Either \code{"Diag"}, \code{"Symm"},  
	\code{"Zero"} or \code{"User"}. If it is
	\code{"Diag"}, a diagonal matrix is used for the random effects
	meaning that the random effects are independent. If it is \code{"Symm"} (default if missing), a
	symmetric matrix is used for the random effects on the covariances
	among the correlation (or covariance) vectors. If it is
	\code{"Zero"}, there is no random effects which is similar to the
	conventional Generalized Least Squares (GLS) approach to
	fixed-effects analysis.
    \code{"User"}, user has to specify the variance component via the
	\code{RE.constraints} argument. This argument will be ignored when \code{method="FEM"}.}
  \item{RE.startvalues}{Starting values on the
	diagonals of the variance component of the random effects. It will be ignored when \code{method="FEM"}.}
  \item{RE.lbound}{Lower bounds on the diagonals of the variance
	component of the random effects. It will be ignored when
	\code{method="FEM"}.}
  \item{RE.constraints}{A \eqn{p*}{p*} x \eqn{p*}{p*} matrix
	specifying the variance components of the random effects, where
	\eqn{p*}{p*} is the number of effect sizes. If the input
	is not a matrix, it is converted into a matrix by
	\code{as.matrix()}. The default is that all
	covariance/variance components are free. The format of this matrix
	follows \code{\link[metaSEM]{as.mxMatrix}}. Elements of the variance
	components can be constrained equally by using the same labels. If a zero matrix is
	specified, it becomes a fixed-effects meta-analysis.
  }  
  \item{I2}{Possible options are \code{"I2q"}, \code{"I2hm"} and
	\code{"I2am"}. They represent the \code{I2} calculated by using a
	typical within-study sampling variance from the Q statistic, the
	harmonic mean and the arithmetic mean of the within-study sampling
	variances (Xiong, Miller, & Morris, 2010). More than one options are possible. If \code{intervals.type="LB"}, 95\% confidence intervals on the heterogeneity indices will be constructed. 
  }
  \item{acov}{If it is \code{individual}, the sampling variance-covariance
	matrices are calculated based on individual correlation/covariance
	matrix. If it is either \code{unweighted} or \code{weighted} (the default), the average
	correlation/covariance matrix is calculated based on the unweighted
	or weighted mean with the sample sizes. The average
	correlation/covariance matrix is used to calculate the sampling
	variance-covariance matrices. This argument is ignored with the
	\code{method="FEM"} argument.}
  \item{model.name}{A string for the model name in \code{\link[OpenMx]{mxModel}}. 
  }  
  \item{suppressWarnings}{Logical. If \code{TRUE}, warnings are
	suppressed. It is passed to \code{\link[OpenMx]{mxRun}}.}
  \item{silent}{Logical. Argument to be passed to \code{\link[OpenMx]{mxRun}}}
  \item{run}{Logical. If \code{FALSE}, only return the mx model without running the analysis.}  
  \item{\dots}{Further arguments to be passed to \code{\link[OpenMx]{mxRun}}
}
}
\value{Either an object of class \code{tssem1FEM} for fixed-effects TSSEM,
  an object of class \code{tssem1FEM.cluster} for fixed-effects TSSEM
  with \code{cluster} argument, or an object of class \code{tssem1REM}
  for random-effects TSSEM. 
}
\references{
  Cheung, M. W.-L. (2014). Fixed- and random-effects meta-analytic structural equation modeling: Examples and analyses in R. \emph{Behavior Research Methods}, \bold{46}, 29-40.

  Cheung, M. W.-L. (2013). Multivariate meta-analysis as structural equation models. \emph{Structural Equation Modeling}, \bold{20},
429-454.

 Cheung, M. W.-L., & Chan, W. (2005). Meta-analytic structural equation modeling: A two-stage approach. \emph{Psychological Methods}, \bold{10}, 40-64.

 Cheung, M. W.-L., & Chan, W. (2009). A two-stage approach to synthesizing covariance matrices in meta-analytic structural equation modeling. \emph{Structural Equation Modeling}, \bold{16}, 28-53.
}
\author{Mike W.-L. Cheung <mikewlcheung@nus.edu.sg>
}
\note{If the \code{cluster} argument is used, it returns a list of
  results on each cluster.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link[metaSEM]{wls}}, \code{\link[metaSEM]{Cheung09}},
    \code{\link[metaSEM]{Becker92}}, \code{\link[metaSEM]{Digman97}}, \code{\link[metaSEM]{issp89}}, \code{\link[metaSEM]{issp05}}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ tssem }
