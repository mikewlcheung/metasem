\name{meta}
\alias{meta}
\title{Fixed-, Random-, and Mixed-Effects Univariate and Multivariate Meta-Analysis
}
\description{
It conducts univariate and multivariate meta-analysis based on maximum
likelihood estimation method. Mixed-effects meta-analysis is possible by
adding study characteristics as predictors. Equality constraints on
intercepts, regression coefficients and variance components are possible.
}
\usage{
meta(y, v, x, intercept.constraints, coeff.constraints, RE.constraints,
RE.startvalues=0.1, RE.lbound = 1e-10, intervals.type = c("z", "LB"), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{An effect size for univariate meta-analysis or a \eqn{k}{k} x
	\eqn{p}{p} matrix of effect sizes for multivariate meta-analysis
	where \eqn{k}{k} is the number of studies and \eqn{p}{p} is the
	number of effect sizes. 
}
  \item{v}{The sampling variance of the effect size for univarite
	meta-analysis or a \eqn{k}{k} x \eqn{p*}{p*} matrix of the sampling
	covariance matrix of the effect sizes for multivariate meta-analysis
	where \eqn{p* = p(p+1)/2 }{p* = p(p+1)/2}. It is arranged by column
	major as used by \code{\link[OpenMx]{vech}}.
}
  \item{x}{A predictor or a \eqn{k}{k} x \eqn{m}{m} matrix of predictors
	where \eqn{m}{m} is the number of predictors. 
}
  \item{intercept.constraints}{A \eqn{1}{1} x \eqn{p}{p} matrix
	specifying whether the intercepts of the effect sizes are fixed or
	free. The default is that the intercepts are free. When there is no predictor, these intercepts are the same as
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
	specifying the variance componets of the random effects. The default is that all
	covariance/variance components are free. The format of this matrix follows \code{\link[metaSEM]{as.mxMatrix}}. The parameter estimates will be
	constrained equally if the labels are the same. If a zero matrix is
	specified, it becomes a fixed-effects meta-analysis.
  }
  \item{RE.startvalues}{A vector of \eqn{p}{p} starting values on the
	diagonals of the variance component of the random effects. If only one
	scalar is given, it will be duplicated across the diagonals. Starting values for the
	off-diagonals of the variance component are all 0. A \eqn{p}{p} x
	\eqn{p}{p} symmetric matrix of starting values is also acceptable. 
}
  \item{RE.lbound}{A vector of \eqn{p}{p} lower bounds on the
	diagonals of the variance component of the random effects. If only one
	scalar is given, it will be duplicated across the diagonals. Lower bounds for the off-diagonals of the variance component are set
	at \code{NA}. A \eqn{p}{p} x
	\eqn{p}{p} symmetric matrix of the lower bounds is also acceptable
}
  \item{intervals.type}{Either \code{z} (default if missing) or
	\code{LB}. If it is \code{z}, it calculates 95\% Wald confidence
	intervals (CIs) based on the z statistic. If it is \code{LB}, it
	calculates 95\% likelihood-based CIs on the
	parameter estimates. Note that the z values and their
	associated p values are based on the z statistic. They are not
	related to the likelihood-based CIs.
}
  \item{\dots}{Futher arguments to be passed to \code{\link[OpenMx]{mxRun}}
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{An object of class \code{meta} with a list of
  \item{call}{Object returned by \code{\link[base]{match.call}}}
  \item{data}{A data matrix of y, v and x  }
  \item{no.y}{No. of effect sizes}
  \item{no.x}{No. of predictors}
  \item{miss.x}{A vector indicating whether the predictors are
	missing. Studies will be removed before the analysis if they are \code{TRUE}}
  \item{meta.fit}{A fitted object returned from
	\code{\link[OpenMx]{mxRun}}}
}
\references{
  Cheung, M. W. L. (2008). A model for integrating fixed-, random-, and
  mixed-effects meta-analyses into structural equation
  modeling. \emph{Psychological Methods}, \bold{13}, 182-202.

  Cheung, M. W. L. (2009). Constructing approximate confidence intervals
  for parameters with structural equation models. \emph{Structural Equation Modeling}, \bold{16}, 267-294.

  Cheung, M. W. L. (2010). Modeling multivariate effect sizes with structural equation models. \emph{Manuscript submitted for publication.}
  
  Hardy, R. J., & Thompson, S. G. (1996). A likelihood approach to
  meta-analysis with random effects. \emph{Statistics in Medicine},
  \bold{15}, 619-629.
  
  Neale, M. C., & Miller, M. B. (1997). The use of likelihood-based
  confidence intervals in genetic models. \emph{Behavior
Genetics}, \bold{27}, 113-120.
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

## Random-effects meta-analysis
attach(Hox02)
summary( meta(y=yi, v=vi) )

## Fixed-effects meta-analysis
summary( meta(y=yi, v=vi, RE.constraints=matrix(0, ncol=1, nrow=1)) )

## Mixed-effects meta-analysis with "weeks" as a predictor
## Use likelihood-based CI
summary( meta(y=yi, v=vi, x=weeks, intervals.type="LB") )
detach(Hox02)


## Multivariate meta-analysis
attach(Berkey98)
summary( meta(y=cbind(PD, AL), v=cbind(var_PD, cov_PD_AL, var_AL)) )

## Multivariate meta-analysis with "publication year-1979" as the predictor
summary( meta(y=cbind(PD, AL), v=cbind(var_PD, cov_PD_AL, var_AL),
         x=scale(pub_year, center=1979)) )

## Multivariate meta-analysis with an equality constraint on regression coefficients
summary( meta(y=cbind(PD, AL), v=cbind(var_PD, cov_PD_AL, var_AL),
         x=scale(pub_year, center=1979),
         coeff.constraints=matrix(c("0.3*Eq_slope", "0.3*Eq_slope"), nrow=2)) )
detach(Berkey98)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ meta-analysis }
