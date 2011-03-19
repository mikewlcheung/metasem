\name{tssem1}
\alias{tssem1}
\title{First Stage of the Two-Stage Structural Equation Modeling (TSSEM)
}
\description{It conducts the first stage analysis of TSSEM by pooling
  correlation/covariance matrices with a fixed-effets model. The function expects that there is no missing data in the first group. 
}
\usage{
tssem1(my.df, n, start.values, cor.analysis = TRUE, model.name,
       suppressWarnings = TRUE, ...)
}
\arguments{
  \item{my.df}{A list of correlation/covariance matrices
}
  \item{n}{A vector of sample sizes
}
 \item{start.values}{A symmetric matrix of starting values for the pooled correlation/covariance matrix. If it is missing, \code{.startValues} will be used to generate the starting values.
}
  \item{cor.analysis}{Logical. The output is either a pooled correlation or a covariance matrix.
  }
  \item{model.name}{A string for the model name in \code{\link[OpenMx]{mxModel}}. If it is missing, the default is "TSSEM1
	Analysis of Correlation Matrix" for \code{cor.analysis=TRUE} and
	"TSSEM1 Analysis of Covariance Matrix" for \code{cor.analysis=FALSE}
	}
  \item{suppressWarnings}{Logical. If \code{TRUE}, warnings are
	suppressed. Argument to be passed to \code{\link[OpenMx]{mxRun}}.}
  \item{\dots}{Further arguments to be passed to \code{\link[OpenMx]{mxRun}}
}
}
\details{
}
\value{
  \item{call}{The matched call}
  \item{data}{A list of correlation/covariance matrices from input}
  \item{pooledS}{The pooled correlation/covariance matrix}
  \item{acovS }{The asymptotic sampling covariance matrix of the pooled
	correlation/covariance matrix}
  \item{total.n}{Total sample size of all studies}
  \item{modelMinus2LL}{-2LogLikelihood of the model}
  \item{independentMinus2LL}{-2LogLikelihood of the independent model
	returned by \code{.minus2LL}}
  \item{saturatedMinus2LL}{-2LogLikelihood of the saturated model
	returned by \code{.minus2LL}}
  \item{tssem1.fit}{A fitted object returned from \code{\link[OpenMx]{mxRun}}}
}
\references{
Cheung, M. W. L., & Chan, W. (2005). Meta-analytic structural equation modeling: A two-stage approach. \emph{Psychological Methods}, \bold{10}, 40-64.

Cheung, M. W. L., & Chan, W. (2009). A two-stage approach to synthesizing covariance matrices in meta-analytic structural equation modeling. \emph{Structural Equation Modeling}, \bold{16}, 28-53.
}
\author{Mike W.-L. Cheung <mikewlcheung@nus.edu.sg>
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link[metaSEM]{wls}}
}
\examples{
digman1 <- tssem1(Digman97$data, Digman97$n)
summary(digman1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ tssem }
