\name{meta}
\alias{meta}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Univariate and multivariate meta-analysis
}
\description{
It conducts univariate and multivariate meta-analysis based on maximum
likelihood estimation. Mixed-effects meta-analysis is allowed by
including study characteristics as predictors. Equality constraints on
intercepts, regression coefficients and variance components are possible.
}
\usage{
meta(y, v, x, intercept.constraints, coeff.constraints, RE.constraints, RE.lbound = 1e-10, intervals = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{An effect size for univariate meta-analysis or a \eqn{k}{k} x
	\eqn{p}{p} matrix of effect sizes for multivariate meta-analysis
	where \eqn{k}{k} is the number of studies and \eqn{p}{p} is the
	number of effect sizes. For multiple effect sizes, they can be combined via \code{\link[base]{cbind}}.
}
  \item{v}{The sampling variance of the effect size for univarite
	meta-analysis or a \eqn{k}{k} x \eqn{p*}{p*} matrix of the sampling
	covariance matrix of the effect sizes for multivariate meta-analysis
	where \eqn{p* = p(p+1)/2 }{p* = p(p+1)/2}. It is arranged by column
	major as used by \code{\link[OpenMx]{vech}}.
}
  \item{x}{A predictor or a \eqn{k}{k} x \eqn{m}{m} matrix of predictors
	where \eqn{m}{m} is the number of predictors. They can be combined via \code{\link[base]{cbind}}.
}
  \item{intercept.constraints}{A \eqn{1}{1} x \eqn{p}{p} matrix
	specifying whether the intercepts of the effect sizes are fixed or
	free. The default is that all intercepts are free. When there is no predictor, these intercepts are the same as
	the pooled effect sizes. The format of this matrix follows
	\code{\link[metaSEM]{as.mxMatrix}}. The parameter estimates will be
	constrained equally if the labels are the same.
}
  \item{coeff.constraints}{A \eqn{p}{p} x \eqn{m}{m} matrix
	specifying how the predictors predict the effect sizes. The default
	is that all \eqn{m}{m} predictors predict all \eqn{p}{p} effect sizes. The
	format of this matrix follows \code{\link[metaSEM]{as.mxMatrix}}. The parameter estimates will be
	constrained equally if the labels are the same.
}
  \item{RE.constraints}{A \eqn{p}{p} x \eqn{p}{p} matrix
	specifying the random effects. The default is that all
	covariance/variance components are free. The format of this matrix follows \code{\link[metaSEM]{as.mxMatrix}}. The parameter estimates will be
	constrained equally if the labels are the same. If a zero matrix is
	specified, it becomes a fixed-effects meta-analysis.
}
  \item{RE.lbound}{It sets the lower bound for the variance components. It
	should be noted that no lower bound is imposed on the covariance
	components for the multivariate meta-analysis
}
  \item{intervals}{logical. If \code{TRUE}, it calculates the likelihood-based confidence intervals on the parameter estimates; otherwise, Wald-based confidence intervals are calculated.
}
  \item{\dots}{Futher arguments to be passed to \code{\link[OpenMx]{mxRun}}
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
  Cheung, M. W. L. (2008). A model for integrating fixed-, random-, and mixed-effects meta-analyses into structural equation modeling. \emph{Psychological Methods}, \bold{13}, 182-202. 
}
\author{Mike W.-L. Cheung <mikewlcheung@nus.edu.sg>
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ meta-analysis }
