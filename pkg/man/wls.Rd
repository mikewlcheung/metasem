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
wls(S, acovS, n, impliedS, matrices, cor.analysis = TRUE,
    intervals.type = c("z", "LB"), model.name, suppressWarnings = TRUE, ...)
tssem2(tssem1.obj, impliedS, matrices, intervals.type = c("z", "LB"),
       model.name, suppressWarnings = TRUE, ...) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tssem1.obj}{An object of either class
	\code{tssemFE1}, class \code{tssemFE1.cluster} or class \code{tssemRE1}}
  \item{S}{A \eqn{p}{p} x \eqn{p}{p} sample correlation/covariance matrix
	where \eqn{p}{p} is the number of variables.
}
  \item{acovS}{A \eqn{p*}{p*} x \eqn{p*}{p*} asymptotic sampling covariance
	matrix of either \code{\link[OpenMx]{vechs}} \code{(S)} or \code{\link[OpenMx]{vech}} \code{(S)} where \eqn{p* = p(p-1)/2 }{p* = p(p-1)/2} for correlation matrix and \eqn{p* = p(p+1)/2 }{p* = p(p+1)/2} for covariance matrix.
}
  \item{n}{Sample size.
}
  \item{impliedS}{Model implied correlation/covariance matrix of an object of either \code{\link[OpenMx]{MxMatrix-class}} or \code{\link[OpenMx]{MxAlgebra-class}}. 
}
  \item{matrices}{A list of matrices used to calculate
	\code{impliedS}. They are objects of either \code{\link[OpenMx]{MxMatrix-class}} or \code{\link[OpenMx]{MxAlgebra-class}}. 
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
\item{indepModelChisq}{Chi-square statistic of the independent model
  returned by \code{.indepwlsChisq} }
\item{indepModelDf}{Degrees of freedom of the independent model returned
by \code{.indepwlsChisq}}
\item{wls.fit}{A fitted object returned from
  \code{\link[OpenMx]{mxRun}}}
}
\references{
Bentler, P.M., & Savalei, V. (2010). Analysis of correlation structures: current status and open problems. In Kolenikov, S., Thombs, L., & Steinley, D. (Eds.). \emph{Recent Methodological Developments in Social Science Statistics} (pp. 1-36). Hoboken, NJ: Wiley. 

  Cheung, M. W. L. (2011). Applications of meta-analytic structural
  equation modeling. \emph{Manuscript in preparation.}

Cheung, M. W. L., & Chan, W. (2005). Meta-analytic structural equation modeling: A two-stage approach. \emph{Psychological Methods}, \bold{10}, 40-64.

Cheung, M. W. L., & Chan, W. (2009). A two-stage approach to synthesizing covariance matrices in meta-analytic structural equation modeling. \emph{Structural Equation Modeling}, \bold{16}, 28-53.

Joreskog, K. G., Sorbom, D., Du Toit, S., & Du Toit, M. (1999). \emph{LISREL 8: New Statistical Features.} Chicago: Scientific Software International. 
}
\author{Mike W.-L. Cheung <mikewlcheung@nus.edu.sg>
}
\note{If the input is a list of \code{tssem1.obj}, it returns a list of
  results for each cluster.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link[metaSEM]{tssem1}}
}
\examples{
#### Analysis of correlation structure
R1 <- matrix(c(1.00, 0.22, 0.24, 0.18,
               0.22, 1.00, 0.30, 0.22,
               0.24, 0.30, 1.00, 0.24,
               0.18, 0.22, 0.24, 1.00), ncol=4, nrow=4)
n <- 1000
acovR1 <- asyCov(R1, n)

## One-factor CFA model- P1: Factor variance; L1: Factor loadings 
P1 <- mxMatrix("Full", ncol=1, nrow=1, value=1, free=FALSE, name="P1")
L1 <- mxMatrix("Full", ncol=1, nrow=4, free=TRUE, name="L1")

## Model to be fitted
impliedR1 <- mxAlgebra(L1 \%&\% P1, name="impliedR1")
wls.fit1 <- wls(S=R1, acovS=acovR1, n=n, impliedS=impliedR1,
                matrices=c(P1, L1), cor.analysis=TRUE)
summary(wls.fit1)

#### Multiple regression analysis with RAM specification
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
A2 <- mxMatrix("Full", ncol=3, nrow=3, byrow=TRUE,
               free=c(FALSE, rep(TRUE, 2), rep(FALSE, 6)), name="A2")

## S2: Covariance matrix of free parameters
## Note that the diagonal elements are not involved in
## the analysis of correlation structure
#    y x1 x2
# y  F F  F 
# x1 F F  F 
# x2 F T  F
S2 <- mxMatrix("Stand", ncol=3, nrow=3, free=c(FALSE, FALSE, TRUE),name="S2")

## Identity matrix
Id <- mxMatrix("Iden", ncol=3, nrow=3, name="Id")

## Model implied correlation matrix: (Id-A2)^-1 \%*\% S2 \%*\% ((Id-A2)^-1)'
impliedR2 <- mxAlgebra( solve(Id-A2) \%&\% S2, name="impliedR2")

wls.fit2 <- wls(S=R2, acovS=acovR2, n=n, impliedS=impliedR2,
                matrices=c(A2, S2, Id), cor.analysis=TRUE,
                model.name="Regression analysis")
summary(wls.fit2)


#### Analysis of covariance structure
S3 <- matrix(c(1.50, 0.22, 0.24, 0.18,
               0.22, 1.60, 0.30, 0.22,
               0.24, 0.30, 1.80, 0.24,
               0.18, 0.22, 0.24, 1.30), ncol=4, nrow=4)
n <- 1000
acovS3 <- asyCov(S3, n, cor.analysis=FALSE)
P3 <- mxMatrix("Full", ncol=1, nrow=1, value=1, free=FALSE, name="P3")
L3 <- mxMatrix("Full", ncol=1, nrow=4, value=c(0.3, 0.4, 0.5, 0.4),
               free=TRUE, name="L3")
E3 <- mxMatrix("Diag", ncol=4, nrow=4, value=0.2, free=TRUE, name="E3") 
impliedS3 <- mxAlgebra(L3 \%&\% P3 + E3, name="impliedS3")

## Use likelihood-based CI
wls.fit3 <- wls(S=S3, acovS=acovS3, n=n, impliedS=impliedS3,
                matrices=c(P3, L3, E3), cor.analysis=FALSE,
                intervals.type="LB",
                model.name="Covariance structure with LB CI")
summary(wls.fit3)


##### Example of Fixed-effects TSSEM
fixed1 <- tssem1(Digman97$data, Digman97$n)
summary(fixed1)

P4 <- mxMatrix("Stand", ncol=2, nrow=2, value=.2, free=TRUE, name="P4")
L4 <- mxMatrix("Full", ncol=2, nrow=5, value=c(0,.3,.3,.3,0,.3,0,0,0,.3),
free=c(FALSE,TRUE,TRUE,TRUE,FALSE,TRUE,FALSE,FALSE,FALSE,TRUE),
name="L4")

## Model to be fitted
impliedR4 <- mxAlgebra(L4 \%&\% P4, name="impliedR4")

fixed2 <- tssem2(fixed1, impliedS=impliedR4, matrices=c(P4, L4),
                  model.name="TSSEM2 Digman97")
summary(fixed2)

##### Example of Fixed-effects TSSEM with several clusters
data.cluster <- tapply(Digman97$data, Digman97$cluster, function(x) {x})
n.cluster <- tapply(Digman97$n, Digman97$cluster, function(x) {x})

fixed1.cluster <- tssem1(Digman97$data, Digman97$n, cluster=Digman97$cluster)
summary(fixed1.cluster)

fixed2.cluster <- tssem2(fixed1.cluster, impliedS=impliedR4, matrices=c(P4, L4))
summary(fixed2.cluster)

##### Example of Random-effects TSSEM
random1 <- tssem1(Digman97$data, Digman97$n, method="RE")
summary(random1)

random2 <- tssem2(random1, impliedS=impliedR4, matrices=c(P4, L4))
summary(random2)
}
\keyword{ tssem }
