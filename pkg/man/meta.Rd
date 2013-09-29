\name{meta}
\alias{meta}
\title{Univariate and Multivariate Meta-Analysis with Maximum
  Likelihood Estimation
}
\description{
It conducts univariate and multivariate meta-analysis with maximum
likelihood estimation method. Mixed-effects meta-analysis can be
conducted by including study characteristics as predictors. Equality constraints on
intercepts, regression coefficients and variance components can be
easily imposed by setting the same labels on the parameter estimates.
}
\usage{
meta(y, v, x, data, intercept.constraints = NULL, coef.constraints = NULL,
     RE.constraints = NULL, RE.startvalues=0.1, RE.lbound = 1e-10,
     intervals.type = c("z", "LB"), I2="I2q", R2=TRUE,
     model.name="Meta analysis with ML", suppressWarnings = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{A vector of effect size for univariate meta-analysis or a \eqn{k}{k} x
	\eqn{p}{p} matrix of effect sizes for multivariate meta-analysis
	where \eqn{k}{k} is the number of studies and \eqn{p}{p} is the
	number of effect sizes. 
}
  \item{v}{A vector of the sampling variance of the effect size for univarite
	meta-analysis or a \eqn{k}{k} x \eqn{p*}{p*} matrix of the sampling
	covariance matrix of the effect sizes for multivariate meta-analysis
	where \eqn{p* = p(p+1)/2 }{p* = p(p+1)/2}. It is arranged by column
	major as used by \code{\link[OpenMx]{vech}}.
}
  \item{x}{A predictor or a \eqn{k}{k} x \eqn{m}{m} matrix of predictors
	where \eqn{m}{m} is the number of predictors. 
}
  \item{data}{An optional data frame containing the variables in the model.
}  
  \item{intercept.constraints}{A \eqn{1}{1} x \eqn{p}{p} matrix
	specifying whether the intercepts of the effect sizes are fixed or
	free. If the input is not a matrix, the input is converted into a
	\eqn{1}{1} x \eqn{p}{p} matrix with
	\code{t(as.matrix(intercept.constraints))}. The default is that the intercepts are free. When there is no predictor, these intercepts are the same as
	the pooled effect sizes. The format of this matrix follows
	\code{\link[metaSEM]{as.mxMatrix}}. The intercepts can be
	constrained equally by using the same labels.
}
  \item{coef.constraints}{A \eqn{p}{p} x \eqn{m}{m} matrix
	specifying how the predictors predict the effect sizes. If the input
	is not a matrix, it is converted into a matrix by
	\code{as.matrix()}. The default is that all \eqn{m}{m} predictors predict all \eqn{p}{p} effect sizes. The
	format of this matrix follows
	\code{\link[metaSEM]{as.mxMatrix}}. The regression coefficients can be
	constrained equally by using the same labels.
}
  \item{RE.constraints}{A \eqn{p}{p} x \eqn{p}{p} matrix
	specifying the variance componets of the random effects. If the input
	is not a matrix, it is converted into a matrix by
	\code{as.matrix()}. The default is that all
	covariance/variance components are free. The format of this matrix
	follows \code{\link[metaSEM]{as.mxMatrix}}. Elements of the variance
	components can be constrained equally by using the same labels. If a zero matrix is
	specified, it becomes a fixed-effects meta-analysis.
  }
  \item{RE.startvalues}{A vector of \eqn{p}{p} starting values on the
	diagonals of the variance component of the random effects. If only one
	scalar is given, it will be duplicated across the diagonals. Starting values for the
	off-diagonals of the variance component are all 0. A \eqn{p}{p} x
	\eqn{p}{p} symmetric matrix of starting values is also accepted. 
}
  \item{RE.lbound}{A vector of \eqn{p}{p} lower bounds on the
	diagonals of the variance component of the random effects. If only one
	scalar is given, it will be duplicated across the diagonals. Lower
	bounds for the off-diagonals of the variance component are set at \code{NA}. A \eqn{p}{p} x
	\eqn{p}{p} symmetric matrix of the lower bounds is also accepted.
}
  \item{intervals.type}{Either \code{z} (default if missing) or
	\code{LB}. If it is \code{z}, it calculates the 95\% Wald confidence
	intervals (CIs) based on the z statistic. If it is \code{LB}, it
	calculates the 95\% likelihood-based CIs on the
	parameter estimates. Note that the z values and their
	associated p values are based on the z statistic. They are not
	related to the likelihood-based CIs.
  }
  \item{I2}{Possible options are \code{"I2q"}, \code{"I2hm"} and
	\code{"I2am"}. They represent the \code{I2} calculated by using a
	typical within-study sampling variance from the Q statistic, the
	harmonic mean and the arithmatic mean of the within-study sampling
	variances (Xiong, Miller, & Morris, 2010). More than one options are possible. If \code{intervals.type="LB"}, 95\% confidence intervals on the heterogeneity indices will be constructed. 
  }
  \item{R2}{Logical. If \code{TRUE} and there are predictors, R2 is
	calculated (Raudenbush, 2009). 
  }
  \item{model.name}{A string for the model name in \code{\link[OpenMx]{mxModel}}. 
	}
  \item{suppressWarnings}{Logical. If \code{TRUE}, warnings are
	suppressed. Argument to be passed to \code{\link[OpenMx]{mxRun}}.}
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
	missing. Studies will be removed before the analysis if they are
	\code{TRUE}}
  \item{I2}{Types of I2 calculated}
  \item{R2}{Logical}
  \item{mx.fit}{A fitted object returned from
	\code{\link[OpenMx]{mxRun}}}
  \item{mx0.fit}{A fitted object without any predictor returned from
	\code{\link[OpenMx]{mxRun}}}
}
\references{
  Cheung, M. W. L. (2008). A model for integrating fixed-, random-, and
  mixed-effects meta-analyses into structural equation
  modeling. \emph{Psychological Methods}, \bold{13}, 182-202.

  Cheung, M. W. L. (2009). Constructing approximate confidence intervals
  for parameters with structural equation models. \emph{Structural
	Equation Modeling}, \bold{16}, 267-294.

  Cheung, M. W. L. (2013). Multivariate meta-analysis as structural equation models. \emph{Structural Equation Modeling}, \bold{20}, 429-454.
  
  Hardy, R. J., & Thompson, S. G. (1996). A likelihood approach to
  meta-analysis with random effects. \emph{Statistics in Medicine},
  \bold{15}, 619-629.
  
  Neale, M. C., & Miller, M. B. (1997). The use of likelihood-based
  confidence intervals in genetic models. \emph{Behavior Genetics},
  \bold{27}, 113-120.

  Raudenbush, S. W. (2009). Analyzing effect sizes: random effects models. In H. M. Cooper, L. V. Hedges, & J. C. Valentine (Eds.), \emph{The handbook of research synthesis and meta-analysis} (2nd ed., pp. 295-315). New York: Russell Sage Foundation.

  Xiong, C., Miller, J. P., & Morris, J. C. (2010). Measuring study-specific heterogeneity in meta-analysis: application to an antecedent biomarker study of alzheimer's disease. \emph{Statistics in Biopharmaceutical Research}, \bold{2(3)}, 300-309. doi:10.1198/sbr.2009.0067
}
\author{Mike W.-L. Cheung <mikewlcheung@nus.edu.sg>
}
\note{Missing values (NA) in y and their related elements in v
  will be removed automatically. When there are missing values in v but
  not in y, missing values will be replaced by 1e5. Effectively, these
  effect sizes will have little impact on the analysis.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link[metaSEM]{reml}}, \code{\link[metaSEM]{Hox02}},
  \code{\link[metaSEM]{Berkey98}}, \code{\link[metaSEM]{wvs94}}
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ meta-analysis }
