\name{wls}
\alias{wls}
\alias{tssem2}
\title{Conduct a Correlation/Covariance Structure Analysis with WLS
}
\description{It fits a correlation or covariance structure with
  weighted least squares (WLS) estimation method where the inverse of the asymptotic covariance matrix is
  used as the weight matrix. \code{tssem2} conducts the second stage
  analysis of the two-stage structural equation modeling (TSSEM). \code{tssem2} is a wrapper of \code{wls}.
}
\usage{
wls(Cov, aCov, n, RAM=NULL, Amatrix=NULL, Smatrix=NULL, Fmatrix=NULL,
    diag.constraints=FALSE, cor.analysis=TRUE, intervals.type=c("z","LB"),
    mx.algebras=NULL, model.name=NULL, suppressWarnings=TRUE,
    silent=TRUE, run=TRUE, ...)
tssem2(tssem1.obj, RAM=NULL, Amatrix=NULL, Smatrix=NULL, Fmatrix=NULL,
    diag.constraints=FALSE, intervals.type=c("z", "LB"), mx.algebras=NULL,
    model.name=NULL, suppressWarnings=TRUE, silent=TRUE, run=TRUE, ...) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tssem1.obj}{An object of either class
	\code{tssem1FEM}, class \code{tssem1FEM.cluster} or class
	\code{tssem1REM} returned from \code{tssem1()}}
  \item{Cov}{A \eqn{p}{p} x \eqn{p}{p} sample correlation/covariance matrix
	where \eqn{p}{p} is the number of variables.
}
  \item{aCov}{A \eqn{p*}{p*} x \eqn{p*}{p*} asymptotic sampling covariance
	matrix of either \code{\link[OpenMx]{vechs}} \code{(Cov)} or
	\code{\link[OpenMx]{vech}} \code{(Cov)} where \eqn{p* = p(p-1)/2 }{p*
	  = p(p-1)/2} for correlation matrix and \eqn{p* = p(p+1)/2 }{p* = p(p+1)/2} for covariance matrix.
}
  \item{n}{Sample size.
}
  \item{RAM}{A RAM object including a list of matrices of the model
  returned from \code{\link[metaSEM]{lavaan2RAM}}.  
}  
  \item{Amatrix}{If \code{RAM} is not specified, an \code{Amatrix} is
  required. An asymmetric matrix in the RAM specification with
	\code{\link[OpenMx]{MxMatrix-class}}. If it is \code{NULL}, a matrix
	of zero will be created. If it is a matrix, it will be converted into \code{\link[OpenMx]{MxMatrix-class}} by the \code{as.mxMatrix} function.
}
  \item{Smatrix}{If \code{RAM} is not specified, an \code{Smatrix} is
  required. A symmetric matrix in the RAM specification with
	\code{\link[OpenMx]{MxMatrix-class}}. If it is a matrix, it will be converted into \code{\link[OpenMx]{MxMatrix-class}} by the \code{as.mxMatrix} function.
}
  \item{Fmatrix}{A filter matrix in the RAM specification with
	\code{\link[OpenMx]{MxMatrix-class}}. If it is \code{NULL} (the default), an
	identity matrix with the same dimensions of \code{Cov} will be
	created. If it is a matrix, it will be converted into \code{\link[OpenMx]{MxMatrix-class}} by the \code{as.mxMatrix} function. It is not required when there is no latent variable.
  }
  \item{diag.constraints}{Logical. This argument is ignored when
	\code{cor.analysis=FALSE}. If
	\code{diag.constraints=TRUE}, the diagonals of the model implied matrix would be constrained at 1 by nonlinear constraints. The drawback is that standard
	error will not be generated. Parametric bootstrap is used to
	estimate the standard error by drawing samples from \eqn{\mathcal{N}(vech(Cov), asyCov)}{N(vech(Cov), asyCov)} for covariance analysis and \eqn{\mathcal{N}(vechs(Cov), asyCov)}{N(vechs(Cov), asyCov)} for
	correlation analysis while asyCov is treated as fixed. This process
	is computationally intensive. A better approach is to request likelihood-based
	confidence intervals (CIs) by specifying \code{intervals.type="LB"}. 
	If \code{diag.constraints=FALSE} and \code{cor.analysis=TRUE}, the diagonals are automatically constrained as ones by
  treating the error variances as computed values rather than as
  parameters. Since the error variances are not parameters, they are not reported. 
}  
  \item{cor.analysis}{Logical. Analysis of correlation or covariance structure. If \code{cor.analysis=TRUE}, \code{\link[OpenMx]{vechs}} is used to vectorize \code{S}; otherwise, \code{\link[OpenMx]{vech}} is used to vectorize \code{S}.
}
  \item{intervals.type}{Either \code{z} (default if missing) or
	\code{LB}. If it is \code{z}, it calculates the 95\% Wald CIs based on the z statistic. If it is \code{LB}, it
	calculates the 95\% likelihood-based CIs on the
	parameter estimates. Please note that the z values and their
	associated p values are based on the z statistic. They are not
	related to the likelihood-based CIs.
  }
  \item{mx.algebras}{A list of \code{\link[OpenMx]{mxMatrix}} or  \code{\link[OpenMx]{mxAlgebra}} objects on the \code{Amatrix},
	\code{Smatrix}, and \code{Fmatrx}. It can be used to define new functions
	of parameters and their LBCIs. For example, if the regression
	coefficients to calculate an indirect effect are stored in A[1,2] and
	A[1,3], we may define
	\code{list(ind=mxAlgebra(Amatrix[1,2]*Amatrix[1,3], name="ind"))}
	See the examples in \code{\link[metaSEM]{Becker92}} and
	\code{\link[metaSEM]{Hunter83}}. It should be noted that Fmatrix,
	Amatrix, Smatrix, Iden (a \eqn{p}{p} x \eqn{p}{p} identity matrix), sampleS (sample correlation
	or covariance matrix), impliedS1, impliedS (model implied
	correlation or covariance matrix), vecS, invAcov, obj, One, select
	and constraint and Ematrix (computed error variances when
	\code{diag.constraints=FALSE}) have been defined internally. You
	should not create new matrices using these names.
  }
  \item{model.name}{A string for the model name in
	\code{\link[OpenMx]{mxModel}}. If it is missing, the default is
	"TSSEM2 (or WLS) Analysis of Correlation Structure" for \code{cor.analysis=TRUE} and
	"TSSEM2 (or WLS) Analysis of Covariance Structure" for \code{cor.analysis=FALSE}.
	}
  \item{suppressWarnings}{Logical. If \code{TRUE}, warnings are
	suppressed. The argument to be passed to \code{\link[OpenMx]{mxRun}}.}
  \item{silent}{Logical. An argument to be passed to \code{\link[OpenMx]{mxRun}}}
  \item{run}{Logical. If \code{FALSE}, only return the mx model without running the analysis.}  
  \item{\dots}{Further arguments to be passed to \code{\link[OpenMx]{mxRun}}.
}
}

\value{An object of class \code{wls} with a list of
\item{call}{The matched call}
\item{Cov}{Input data of either a covariance or correlation matrix}
\item{asyCov}{The asymptotic covariance matrix of the input data}
\item{noObservedStat}{Number of observed statistics}
\item{n}{Sample size}
\item{cor.analysis}{logical}
\item{noConstraints}{Number of constraints imposed on S}
\item{indepModelChisq}{Chi-square statistic of the independent model
  returned by \code{.indepwlsChisq} }
\item{indepModelDf}{Degrees of freedom of the independent model returned
by \code{.indepwlsChisq}}
\item{mx.fit}{A fitted object returned from
  \code{\link[OpenMx]{mxRun}}}
}
\references{
Bentler, P.M., & Savalei, V. (2010). Analysis of correlation structures: current status and open problems. In Kolenikov, S., Thombs, L., & Steinley, D. (Eds.). \emph{Recent Methodological Developments in Social Science Statistics} (pp. 1-36). Hoboken, NJ: Wiley. 

Cheung, M. W.-L. (2010). Fixed-effects meta-analyses as multiple-group structural equation models. \emph{Structural Equation Modeling}, \bold{17}, 481-509.

Cheung, M. W.-L. (2014). Fixed- and random-effects meta-analytic structural equation modeling: Examples and analyses in R. \emph{Behavior Research Methods}, \bold{46}, 29-40.

Cheung, M. W.-L., & Chan, W. (2005). Meta-analytic structural equation modeling: A two-stage approach. \emph{Psychological Methods}, \bold{10}, 40-64.

Cheung, M. W.-L., & Chan, W. (2009). A two-stage approach to synthesizing covariance matrices in meta-analytic structural equation modeling. \emph{Structural Equation Modeling}, \bold{16}, 28-53.

Joreskog, K. G., Sorbom, D., Du Toit, S., & Du Toit,
M. (1999). \emph{LISREL 8: New Statistical Features.} Chicago: Scientific Software International.

McArdle, J. J., & MacDonald, R. P. (1984). Some algebraic properties of the Reticular Action Model for moment structures. \emph{British Journal of Mathematical and Statistical Psychology}, \bold{37}, 234-251. 
}
\author{Mike W.-L. Cheung <mikewlcheung@nus.edu.sg>
}
\note{If the input is a list of \code{tssem1.obj}, it returns a list of
  results for each cluster.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link[metaSEM]{tssem1}},
  \code{\link[metaSEM]{Becker92}}, \code{\link[metaSEM]{Digman97}},
  \code{\link[metaSEM]{Hunter83}}, \code{\link[metaSEM]{issp89}}, \code{\link[metaSEM]{issp05}}
}
\examples{
#### Analysis of correlation structure
R1.labels <- c("a1", "a2", "a3", "a4")

R1 <- matrix(c(1.00, 0.22, 0.24, 0.18,
               0.22, 1.00, 0.30, 0.22,
               0.24, 0.30, 1.00, 0.24,
               0.18, 0.22, 0.24, 1.00), ncol=4, nrow=4,
               dimnames=list(R1.labels, R1.labels))
n <- 1000
acovR1 <- asyCov(R1, n)

#### One-factor CFA model using lavaan specification
model1 <- "f =~ a1 + a2 + a3 + a4"

RAM1 <- lavaan2RAM(model1, obs.variables=R1.labels)

wls.fit1a <- wls(Cov=R1, aCov=acovR1, n=n, RAM=RAM1,
                 cor.analysis=TRUE, intervals="LB")
summary(wls.fit1a)

## One-factor CFA model using RAM specification
(A1 <- cbind(matrix(0, nrow=5, ncol=4),
             matrix(c("0.2*a1","0.2*a2","0.2*a3","0.2*a4",0),
             ncol=1)))

(S1 <- Diag(c("0.2*e1","0.2*e2","0.2*e3","0.2*e4",1)))

## The first 4 variables are observed while the last one is latent.
(F1 <- create.Fmatrix(c(1,1,1,1,0), name="F1"))

wls.fit1b <- wls(Cov=R1, aCov=acovR1, n=n, Fmatrix=F1, Smatrix=S1, Amatrix=A1,
                 cor.analysis=TRUE, intervals="LB")
summary(wls.fit1b)

#### Multiple regression analysis using lavaan specification
R2.labels <- c("y", "x1", "x2")

R2 <- matrix(c(1.00, 0.22, 0.24, 
               0.22, 1.00, 0.30, 
               0.24, 0.30, 1.00, 
               0.18, 0.22, 0.24), ncol=3, nrow=3,
               dimnames=list(R2.labels, R2.labels))
acovR2 <- asyCov(R2, n)

model2 <- "y ~ x1 + x2
           ## Variances of x1 and x2 are 1
           x1 ~~ 1*x1
           x2 ~~ 1*x2
           ## x1 and x2 are correlated
           x1 ~~ x2"

RAM2 <- lavaan2RAM(model2, obs.variables=R2.labels)

wls.fit2a <- wls(Cov=R2, aCov=acovR2, n=n, RAM=RAM2,
                 cor.analysis=TRUE, intervals="LB")
summary(wls.fit2a)


#### Multiple regression analysis using RAM specification

## A2: Regression coefficents
#    y x1 x2
# y  F T  T 
# x1 F F  F 
# x2 F F  F 
(A2 <- mxMatrix("Full", ncol=3, nrow=3, byrow=TRUE,
               free=c(FALSE, rep(TRUE, 2), rep(FALSE, 6)), name="A2"))

## S2: Covariance matrix of free parameters
#    y x1 x2
# y  T F  F 
# x1 F F  F 
# x2 F T  F
(S2 <- mxMatrix("Symm", ncol=3, nrow=3, values=c(0.2,0,0,1,0.2,1),
                labels=c("Var_y", NA, NA, NA, "Cov_x1_x2", NA),
                free=c(TRUE,FALSE,FALSE,FALSE,TRUE,FALSE), name="S2"))

## F may be ignored as there is no latent variable.
wls.fit2b <- wls(Cov=R2, aCov=acovR2, n=n, Amatrix=A2, Smatrix=S2,
                 cor.analysis=TRUE, intervals="LB")
summary(wls.fit2b)


#### Analysis of covariance structure using lavaan specification
R3.labels=c("a1", "a2", "a3", "a4")

R3 <- matrix(c(1.50, 0.22, 0.24, 0.18,
               0.22, 1.60, 0.30, 0.22,
               0.24, 0.30, 1.80, 0.24,
               0.18, 0.22, 0.24, 1.30), ncol=4, nrow=4,
               dimnames=list(R3.labels, R3.labels))
n <- 1000
acovS3 <- asyCov(R3, n, cor.analysis=FALSE)

model3 <- "f =~ a1 + a2 + a3 + a4"

RAM3 <- lavaan2RAM(model3, obs.variables=R3.labels)

wls.fit3a <- wls(Cov=R3, aCov=acovS3, n=n, RAM=RAM3,
                 cor.analysis=FALSE)
summary(wls.fit3a)

#### Analysis of covariance structure using RAM specification
(A3 <- cbind(matrix(0, nrow=5, ncol=4),
             matrix(c("0.2*a1","0.2*a2","0.2*a3","0.2*a4",0),ncol=1)))

(S3 <- Diag(c("0.2*e1","0.2*e2","0.2*e3","0.2*e4",1)))

F3 <- c(TRUE,TRUE,TRUE,TRUE,FALSE)
(F3 <- create.Fmatrix(F3, name="F3", as.mxMatrix=FALSE))

wls.fit3b <- wls(Cov=R3, aCov=acovS3, n=n, Amatrix=A3, Smatrix=S3,
                Fmatrix=F3, cor.analysis=FALSE)
summary(wls.fit3b)
}
\keyword{ tssem }
