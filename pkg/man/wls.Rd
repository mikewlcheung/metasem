\name{wls}
\alias{wls}
\alias{tssem2}
\title{Weighted least squares on correlation or covariance structure
}
\description{Analysis of correlation or covariance structure with
  weighted least squares (WLS) where the asymptotic covariance matrix is
  used as the weight matrix. \code{tssem2} conducts the second stage of
  two-stage strutural equation modeling (TSSEM). \code{tssem2} is simply a wrapper of \code{wls}.
}
\usage{
wls(S, acovS, n, impliedS, matrices, cor.analysis = TRUE, intervals = FALSE, ...)
tssem2(tssem1.obj, impliedS, matrices, intervals = FALSE, ...) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tssem1.obj}{An object returned from \code{\link[metaSEM]{tssem1}}}
  \item{S}{A \eqn{p}{p} x \eqn{p}{p} sample correlation/covariance matrix
	where \eqn{p}{p} is the number of variables
}
  \item{acovS}{A \eqn{p*}{p*} x \eqn{p*}{p*} asymptotic sampling covariance
	matrix of either \code{\link[OpenMx]{vechs}} \code{(S)} or \code{\link[OpenMx]{vech}} \code{(S)} where \eqn{p* = p(p-1)/2 }{p* = p(p-1)/2} for correlation matrix and \eqn{p* = p(p+1)/2 }{p* = p(p+1)/2} for covariance matrix
}
  \item{n}{Sample size
}
  \item{impliedS}{Model implied correlation/covariance matrix of an object of either \code{\link[OpenMx]{MxMatrix-class}} or \code{\link[OpenMx]{MxAlgebra-class}} 
}
  \item{matrices}{A vector of matrices used to calculate
	\code{impliedS}. They are objects of either \code{\link[OpenMx]{MxMatrix-class}} or \code{\link[OpenMx]{MxAlgebra-class}} 
}
  \item{cor.analysis}{Logical. Analysis of correlation or covariance structure. If \code{cor.analysis=TRUE}, \code{\link[OpenMx]{vechs}} is used to vectorize \code{S}; otherwise, \code{\link[OpenMx]{vech}} is used to vectorize \code{S}.
}
  \item{intervals}{logical. If \code{TRUE}, it requests likelihood-based confidence intervals on the parameter estimates; otherwise, Wald-based confidence intervals are calculated.
}
  \item{\dots}{Futher arguments to be passed to \code{\link[OpenMx]{mxRun}}
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{An object of class \code{wls} with a list of
\item{call}{The matched call}
\item{noObservedStat}{Number of observed statistics}
\item{n}{Sample size}
\item{indepModelChisq}{Chi-square statistic of the independent model}
\item{indepModelDf}{Degrees of freedom of the independent model}
\item{wls.fit}{A fitted object returned from \code{\link[OpenMx]{mxRun}}}
}
\references{
%Bentler, P.M., & Savalei, V. (2010). Analysis of correlation structures: current status and open problems. In Kolenikov, S., Thombs, L., & Steinley, D. (Eds.). \emph{Recent Methodological Developments in Social Science Statistics} (pp. 1-36). Hoboken, NJ: Wiley. 
  
Cheung, M. W. L., & Chan, W. (2005). Meta-analytic structural equation modeling: A two-stage approach. \emph{Psychological Methods}, \bold{10}, 40-64.

Cheung, M. W. L., & Chan, W. (2009). A two-stage approach to synthesizing covariance matrices in meta-analytic structural equation modeling. \emph{Structural Equation Modeling}, \bold{16}, 28-53.

Joreskog, K. G., Sorbom, D., Du Toit, S., & Du Toit, M. (1999). \emph{LISREL 8: New Statistical Features.} Chicago: Scientific Software International. 
}
\author{Mike W.-L. Cheung <mikewlcheung@nus.edu.sg>
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link[metaSEM]{indepwlsChisq}}, \code{\link[metaSEM]{tssem1}}
}
\examples{
\dontrun{
# Analysis of correlation structure
R1 <- matrix(c(1.00, 0.22, 0.24, 0.18,
               0.22, 1.00, 0.30, 0.22,
               0.24, 0.30, 1.00, 0.24,
               0.18, 0.22, 0.24, 1.00), ncol=4, nrow=4)
n <- 1000
acovR <- asyCov(R1, n)
P1 <- mxMatrix("Full", ncol=1, nrow=1, value=1, free=FALSE, name="P1")
L1 <- mxMatrix("Full", ncol=1, nrow=4, value=c(0.3, 0.4, 0.5, 0.4), free=TRUE, name="L1")
impliedR1 <- mxAlgebra(L1 \%&\% P1, name="impliedR1")
wls.fit1 <- wls(S=R1, acovS=acovR, n=n, impliedS=impliedR1, matrices=c(P1, L1), cor.analysis=TRUE)
summary(wls.fit1)

# Analysis of covariance structure
S2 <- matrix(c(1.50, 0.22, 0.24, 0.18,
               0.22, 1.60, 0.30, 0.22,
               0.24, 0.30, 1.80, 0.24,
               0.18, 0.22, 0.24, 1.30), ncol=4, nrow=4)
n <- 1000
acovS <- asyCov(S2, n, cor.analysis=FALSE)
P2 <- mxMatrix("Full", ncol=1, nrow=1, value=1, free=FALSE, name="P2")
L2 <- mxMatrix("Full", ncol=1, nrow=4, value=c(0.3, 0.4, 0.5, 0.4), free=TRUE, name="L2")
E2 <- mxMatrix("Diag", ncol=4, nrow=4, value=0.2, free=TRUE, name="E2") 
impliedS2 <- mxAlgebra(L2 \%&\% P2 + E2, name="impliedS2")
wls.fit2 <- wls(S=S2, acovS=acovS, n=n, impliedS=impliedS2, matrices=c(P2, L2, E2), cor.analysis=FALSE, intervals=TRUE)
summary(wls.fit2)

# Example of tssem2
digman1 <- tssem1(Digman97, Digman97.n)
P3 <- mxMatrix("Stand", ncol=2, nrow=2, value=.2, free=TRUE, name="P3")
L3 <- mxMatrix("Full", ncol=2, nrow=5, value=c(0,.3,.3,.3,0,.3,0,0,0,.3),
free=c(FALSE,TRUE,TRUE,TRUE,FALSE,TRUE,FALSE,FALSE,FALSE,TRUE), name="L3")
impliedR3 <- mxAlgebra(L3 \%&\% P3, name="impliedR3")
digman2 <- tssem2(digman1, impliedS=impliedR3, matrices=c(P3, L3))
summary(digman2)
}
}
\keyword{ tssem }
