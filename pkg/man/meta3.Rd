\name{meta3}
\alias{meta3}
\title{Three-Level Univariate Meta-Analysis with Maximum Likelihood Estimation
}
\description{
It conducts three-level univariate meta-analysis with maximum
likelihood estimation method. Mixed-effects meta-analysis can be conducted by
including study characteristics as predictors. Equality constraints on
the intercepts, regression coefficients and variance components on the level-2
and on the level-3 can be easily imposed by setting the same labels on the parameter estimates.
}
\usage{
meta3(y, v, cluster, x, data, intercept.constraints,
      coef.constraints, RE2.constraints, RE2.lbound = 1e-10,
      RE3.constraints, RE3.lbound = 1e-10,
      intervals.type = c("z", "LB"), I2="I2q",
      R2=TRUE, model.name = "Meta analysis with ML",
      suppressWarnings = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{A vector of \eqn{k}{k} studies of effect size. 
}
  \item{v}{A vector of \eqn{k}{k} studies of sampling variance.
}
  \item{cluster}{A vector of \eqn{k}{k} characters or numbers indicating the
	clusters. 
}
  \item{x}{A predictor or a \eqn{k}{k} x \eqn{m}{m} matrix of level-2
	and level-3 predictors where \eqn{m}{m} is the number of predictors. 
}
  \item{data}{An optional data frame containing the variables in the model.
}
  \item{intercept.constraints}{A \eqn{1}{1} x \eqn{p}{p} matrix
	specifying whether the intercepts of the effect sizes are fixed or
	free. The default is that the intercepts are free. When there is no predictor, these intercepts are the same as
	the pooled effect sizes. The format of this matrix follows
	\code{\link[metaSEM]{as.mxMatrix}}. The intercepts can be
	constrained equally by using the same labels.
}
  \item{coef.constraints}{A \eqn{p}{p} x \eqn{m}{m} matrix
	specifying how the level-2 and level-3 predictors predict the effect sizes. The default
	is that all \eqn{m}{m} predictors predict all \eqn{p}{p} effect sizes. The
	format of this matrix follows
	\code{\link[metaSEM]{as.mxMatrix}}. The regression coefficients can be
	constrained equally by using the same labels.
}
  \item{RE2.constraints}{A \eqn{p}{p} x \eqn{p}{p} matrix
	specifying the variance componets of the random effects. The default is that all
	covariance/variance components are free. The format of this matrix
	follows \code{\link[metaSEM]{as.mxMatrix}}. Elements of the variance
	components can be constrained equally by using the same labels. If a zero matrix is
	specified, it becomes a fixed-effects meta-analysis.
}
  \item{RE2.lbound}{A vector of \eqn{p}{p} lower bounds on the
	diagonals of the variance component of the random effects at level-2. If only one
	scalar is given, it will be duplicated across the diagonals. Lower
	bounds for the off-diagonals of the variance component are set at \code{NA}. A \eqn{p}{p} x
	\eqn{p}{p} symmetric matrix of the lower bounds is also accepted.
}
  \item{RE3.constraints}{A \eqn{p}{p} x \eqn{p}{p} matrix
	specifying the variance componets of the random effects at level-3. The default is that all
	covariance/variance components are free. The format of this matrix
	follows \code{\link[metaSEM]{as.mxMatrix}}. Elements of the variance
	components can be constrained equally by using the same labels. If a zero matrix is
	specified, it becomes a fixed-effects meta-analysis.
}
  \item{RE3.lbound}{A vector of \eqn{p}{p} lower bounds on the
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
  \item{I2}{Possible options are \code{"I2q"}, \code{"I2hm"}, \code{"I2am"}
	and \code{"ICC"}. More than one options are possible. If \code{intervals.type="LB"}, 95\% confidence intervals on the heterogeneity indices will be constructed. (See the details.)
  }
  \item{R2}{Logical. If \code{TRUE} and there are predictors, R2 are
	calculated.
  }
  \item{model.name}{A string for the model name in \code{\link[OpenMx]{mxModel}}. 
}
  \item{suppressWarnings}{Logical. If \code{TRUE}, warnings are
	suppressed. Argument to be passed to \code{\link[OpenMx]{mxRun}}.
}
  \item{\dots}{Futher arguments to be passed to \code{\link[OpenMx]{mxRun}}}
}
\details{
  \deqn{y_{ij} = \beta_0 + \beta_1*x_{ij} + \dots + \beta_m*x_{ij} + u_{ij(2)}
	+ u_{ij(3)} + e_{ij}
	} where \eqn{u_{ij(2)} \sim N(0, \tau^2_2)}{u_{ij(2)}~ N(0,
	  tau^2_2)} and \eqn{u_{ij(3)} \sim N(0, \tau^2_3)}{u_{ij(3)}~ N(0,
	  tau^2_3)} are the level-2 and level-3 heterogeneity variances,
	respectively, and \eqn{e_{ij} \sim N(0, v_{ij})}{e_{ij}~ N(0,
	  v_{ij})} is the conditional known sampling variance.
	
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
  Cheung, M. W. L. (2012). Three-level meta-analyses as structural
  equation models. \emph{Manuscript in preparation}.

  Konstantopoulos, S. (2011). Fixed effects and variance components
  estimation in three-level meta-analysis. \emph{Research Synthesis
  Methods}, \bold{2}, 61-76.
}
\author{Mike W.-L. Cheung <mikewlcheung@nus.edu.sg>
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link[metaSEM]{Cooper03}}
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ meta-analysis }