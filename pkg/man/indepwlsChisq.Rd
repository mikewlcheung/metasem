\name{indepwlsChisq}
\alias{indepwlsChisq}
\title{Calculate chi-square statistic of an independent model based on WLS
}
\description{It calculates the chi-square statistic of an independent model based on WLS. The analysis can be either a correlation or a covariance structure. It is usually called internally by \code{wls}.
}
\usage{
indepwlsChisq(S, acovS, cor.analysis = TRUE)
}
\arguments{
  \item{S}{A \eqn{p}{p} x \eqn{p}{p} sample correlation/covariance matrix
	where \eqn{p}{p} is the number of variables
}
  \item{acovS}{A \eqn{p*}{p*} x \eqn{p*}{p*} asymptotic sampling covariance
	matrix of either \code{\link[OpenMx]{vechs}} \code{(S)} or \code{\link[OpenMx]{vech}} \code{(S)} where \eqn{p* = p(p-1)/2 }{p* = p(p-1)/2} for correlation matrix and \eqn{p* = p(p+1)/2 }{p* = p(p+1)/2} for covariance matrix
}
  \item{cor.analysis}{Logical. Analysis of correlation or covariance structure. If \code{cor.analysis=TRUE}, \code{\link[OpenMx]{vechs}} is used to vectorize \code{S}; otherwise, \code{\link[OpenMx]{vech}} is used to vectorize \code{S}.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{Chi-square statistic of an independent model on either a correlation or covariance structure with WLS.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{Mike W.-L. Cheung <mikewlcheung@nus.edu.sg>
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link[metaSEM]{wls}}, \code{\link[metaSEM]{minus2LL}}
}
\examples{
S <- matrix(c(1.00, 0.22, 0.24, 0.18,
              0.22, 1.00, 0.30, 0.22,
              0.24, 0.30, 1.00, 0.24,
              0.18, 0.22, 0.24, 1.00), ncol=4, nrow=4)
n <- 1000
acovS <- asyCov(S, n)
indepwlsChisq(S=S, acovS=acovS, cor.analysis=TRUE)
}
\keyword{ utilities }
