\name{wls}
\alias{wls}
\alias{tssem2}
\title{Conduct a Correlation/Covariance Structure Analysis with WLS
}
\description{It fits a correlation or covariance structure with
  weighted least squares (WLS) estimation method where the inverse of the asymptotic covariance matrix is
  used as the weight matrix. \code{tssem2} conducts the second stage
  analysis of the two-stage strutural equation modeling (TSSEM). \code{tssem2} is a wrapper of \code{wls}.
}
\usage{
wls(Cov, asyCov, n, Amatrix=NULL, Smatrix=NULL, Fmatrix=NULL, diag.constraints=FALSE,
    cor.analysis=TRUE, intervals.type=c("z", "LB"), model.name=NULL,
    suppressWarnings=TRUE, ...)
tssem2(tssem1.obj, Amatrix=NULL, Smatrix=NULL, Fmatrix=NULL, diag.constraints=FALSE,
       intervals.type=c("z", "LB"), model.name=NULL, suppressWarnings=TRUE, ...) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tssem1.obj}{An object of either class
	\code{tssem1FEM}, class \code{tssem1FEM.cluster} or class \code{tssem1REM}}
  \item{Cov}{A \eqn{p}{p} x \eqn{p}{p} sample correlation/covariance matrix
	where \eqn{p}{p} is the number of variables.
}
  \item{asyCov}{A \eqn{p*}{p*} x \eqn{p*}{p*} asymptotic sampling covariance
	matrix of either \code{\link[OpenMx]{vechs}} \code{(S)} or \code{\link[OpenMx]{vech}} \code{(S)} where \eqn{p* = p(p-1)/2 }{p* = p(p-1)/2} for correlation matrix and \eqn{p* = p(p+1)/2 }{p* = p(p+1)/2} for covariance matrix.
}
  \item{n}{Sample size.
}
  \item{Amatrix}{An asymmetric matrix in the RAM approach with
	\code{\link[OpenMx]{MxMatrix-class}}. If it is \code{NULL}, a matrix
	with zeros will be created.
}
  \item{Smatrix}{A symmetric matrix in the RAM approach with \code{\link[OpenMx]{MxMatrix-class}}.
}
  \item{Fmatrix}{A filter matrix in the RAM approach with
	\code{\link[OpenMx]{MxMatrix-class}}. If it is \code{NULL} (the default), an
	identity matrix with the same dimensions of \code{S} will be created.
}  
\item{diag.constraints}{Logical. If \code{cor.analysis=TRUE}, this
  argument constrains the diagonals at 1. The drawback is that standard
  error will not be generated. Users may use \code{intervals.type="LB"}
  to obtain the likelihood-based confidence intervals of the parameter estimates.
}
  \item{cor.analysis}{Logical. Analysis of correlation or covariance structure. If \code{cor.analysis=TRUE}, \code{\link[OpenMx]{vechs}} is used to vectorize \code{S}; otherwise, \code{\link[OpenMx]{vech}} is used to vectorize \code{S}.
}
  \item{intervals.type}{Either \code{z} (default if missing) or
	\code{LB}. If it is \code{z}, it calculates the 95\% Wald confidence
	intervals (CIs) based on the z statistic. If it is \code{LB}, it
	calculates the 95\% likelihood-based CIs on the
	parameter estimates. Please note that the z values and their
	associated p values are based on the z statistic. They are not
	related to the likelihood-based CIs.
  }
  \item{model.name}{A string for the model name in
	\code{\link[OpenMx]{mxModel}}. If it is missing, the default is
	"TSSEM2 (or WLS) Analysis of Correlation Structure" for \code{cor.analysis=TRUE} and
	"TSSEM2 (or WLS) Analysis of Covariance Structure" for \code{cor.analysis=FALSE}.
	}
  \item{suppressWarnings}{Logical. If \code{TRUE}, warnings are
	suppressed. Argument to be passed to \code{\link[OpenMx]{mxRun}}.}
  \item{\dots}{Futher arguments to be passed to \code{\link[OpenMx]{mxRun}}.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{An object of class \code{wls} with a list of
\item{call}{The matched call}
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

Cheung, M. W. L. (2010). Fixed-effects meta-analyses as multiple-group structural equation models. \emph{Structural Equation Modeling}, \bold{17}, 481-509.

Cheung, M. W. L. (2011). Applications of meta-analytic structural equation modeling. \emph{Manuscript in preparation.}

Cheung, M. W. L., & Chan, W. (2005). Meta-analytic structural equation modeling: A two-stage approach. \emph{Psychological Methods}, \bold{10}, 40-64.

Cheung, M. W. L., & Chan, W. (2009). A two-stage approach to synthesizing covariance matrices in meta-analytic structural equation modeling. \emph{Structural Equation Modeling}, \bold{16}, 28-53.

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
  \code{\link[metaSEM]{Becker92}}, \code{\link[metaSEM]{Digman97}}, \code{\link[metaSEM]{issp89}}, \code{\link[metaSEM]{issp05}}
}
\examples{
#### Analysis of correlation structure
R1 <- matrix(c(1.00, 0.22, 0.24, 0.18,
               0.22, 1.00, 0.30, 0.22,
               0.24, 0.30, 1.00, 0.24,
               0.18, 0.22, 0.24, 1.00), ncol=4, nrow=4)
n <- 1000
acovR1 <- asyCov(R1, n)

## One-factor CFA model
(A1 <- cbind(matrix(0, nrow=5, ncol=4),
             matrix(c("0.2*a1","0.2*a2","0.2*a3","0.2*a4",0),
             ncol=1)))
A1 <- as.mxMatrix(A1)

(S1 <- diag(c("0.2*e1","0.2*e2","0.2*e3","0.2*e4",1)))
S1 <- as.mxMatrix(S1)

## The first 4 variables are observed while the last one is latent.
(F1 <- create.Fmatrix(c(1,1,1,1,0), name="F1"))
wls.fit1 <- wls(Cov=R1, asyCov=acovR1, n=n, Fmatrix=F1, Smatrix=S1, Amatrix=A1,
                 cor.analysis=TRUE, diag.constraints=TRUE, intervals="LB")
summary(wls.fit1)


#### Multiple regression analysis
## Variables in R2: y, x1, x2
R2 <- matrix(c(1.00, 0.22, 0.24, 
               0.22, 1.00, 0.30, 
               0.24, 0.30, 1.00, 
               0.18, 0.22, 0.24), ncol=3, nrow=3)
acovR2 <- asyCov(R2, n)

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
               free=c(TRUE,FALSE,FALSE,FALSE,TRUE,FALSE), name="S2"))

## F may be ignored as there is no latent variable.
wls.fit2 <- wls(Cov=R2, asyCov=acovR2, n=n, Amatrix=A2, Smatrix=S2,
                cor.analysis=TRUE, diag.constraints=TRUE, intervals="LB")
summary(wls.fit2)


#### Analysis of covariance structure
R3 <- matrix(c(1.50, 0.22, 0.24, 0.18,
               0.22, 1.60, 0.30, 0.22,
               0.24, 0.30, 1.80, 0.24,
               0.18, 0.22, 0.24, 1.30), ncol=4, nrow=4)
n <- 1000
acovS3 <- asyCov(R3, n, cor.analysis=FALSE)

(A3 <- cbind(matrix(0, nrow=5, ncol=4), matrix(c("0.2*a1","0.2*a2","0.2*a3","0.2*a4",0),ncol=1)))
A3 <- as.mxMatrix(A3)

(S3 <- diag(c("0.2*e1","0.2*e2","0.2*e3","0.2*e4",1)))
S3 <- as.mxMatrix(S3)

F3 <- c(TRUE,TRUE,TRUE,TRUE,FALSE)
(F3 <- create.Fmatrix(F3, name="F3"))

wls.fit3 <- wls(Cov=R3, asyCov=acovS3, n=n, Amatrix=A3, Smatrix=S3,
                Fmatrix=F3, cor.analysis=FALSE)
summary(wls.fit3)
}
\keyword{ tssem }
