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
  \item{intervals}{logical. Request likelihood-based confidence intervals on the parameter estimates.
}
  \item{\dots}{Futher arguments to pass to \code{\link[OpenMx]{mxRun}}
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{A fitted object returned from \code{\link[OpenMx]{mxRun}}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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

\seealso{ \code{\link[metaSEM]{tssem1}}
}
\examples{
\dontrun{
# Example of wls
my.cor <- matrix(c(1.00, 0.22, 0.24, 0.18,
                   0.22, 1.00, 0.30, 0.22,
                   0.24, 0.30, 1.00, 0.24,
                   0.18, 0.22, 0.24, 1.00), ncol=4, nrow=4)
n <- 1000
acovS <- asyCov(my.cor, n)
P1 <- mxMatrix("Full", ncol=1, nrow=1, value=1, free=FALSE, name="P1")
L1 <- mxMatrix("Full", ncol=1, nrow=4, value=c(0.3, 0.4, 0.5, 0.4), free=TRUE, name="L1")
impliedR1 <- mxAlgebra(L1 \%&\% P1, name="impliedR1")
wls.fit <- wls(S=my.cor, acovS=acovS, n=n, impliedS=impliedR1,matrices=c(P1, L1), cor.analysis=TRUE)

# Example of tssem2
digman1 <- tssem1(Digman97, Digman97.n)
P2 <- mxMatrix("Stand", ncol=2, nrow=2, value=.2, free=TRUE, name="P2")
L2 <- mxMatrix("Full", ncol=2, nrow=5,
value=c(0,.3,.3,.3,0,.3,0,0,0,.3),
free=c(FALSE,TRUE,TRUE,TRUE,FALSE,TRUE,FALSE,FALSE,FALSE,TRUE), name="L2")
impliedR2 <- mxAlgebra(L2 \%&\% P2, name="impliedR2")
digman2 <- tssem2(digman1, impliedS=impliedR2, matrices=c(P2, L2))
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ tssem }
