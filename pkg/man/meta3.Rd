\name{meta3}
\alias{meta3}
\alias{meta3X}
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
meta3(y, v, cluster, x, data, intercept.constraints = NULL,
      coef.constraints = NULL , RE2.constraints = NULL,
      RE2.lbound = 1e-10, RE3.constraints = NULL, RE3.lbound = 1e-10,
      intervals.type = c("z", "LB"), I2="I2q",
      R2=TRUE, model.name = "Meta analysis with ML",
      suppressWarnings = TRUE, ...)
meta3X(y, v, cluster, x2, x3, av2, av3, data, intercept.constraints=NULL,
       coef.constraints=NULL, RE2.constraints=NULL, RE2.lbound=1e-10,
       RE3.constraints=NULL, RE3.lbound=1e-10, intervals.type=c("z", "LB"),
       R2=TRUE, model.name="Meta analysis with ML", suppressWarnings=TRUE, ...) 
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
  \item{x2}{A predictor or a \eqn{k}{k} x \eqn{m}{m} matrix of level-2
	predictors where \eqn{m}{m} is the number of predictors. 
}
  \item{x3}{A predictor or a \eqn{k}{k} x \eqn{m}{m} matrix of level-3
	predictors where \eqn{m}{m} is the number of predictors. 
}
  \item{av2}{A predictor or a \eqn{k}{k} x \eqn{m}{m} matrix of level-2
	auxiliary variables where \eqn{m}{m} is the number of variables. 
}
  \item{av3}{A predictor or a \eqn{k}{k} x \eqn{m}{m} matrix of level-3
	auxiliary variables where \eqn{m}{m} is the number of variables. 
}    
  \item{data}{An optional data frame containing the variables in the model.
}
  \item{intercept.constraints}{A \eqn{1}{1} x \eqn{1}{1} matrix
	specifying whether the intercept of the effect size is fixed or
	constrained. The format of this matrix follows
	\code{\link[metaSEM]{as.mxMatrix}}. The intercept can be
	constrained with other parameters by using the same label.
}
  \item{coef.constraints}{A \eqn{1}{1} x \eqn{m}{m} matrix
	specifying how the level-2 and level-3 predictors predict the effect sizes. If the input
	is not a matrix, it is converted into a matrix by
	\code{as.matrix()}. The default
	is that all \eqn{m}{m} predictors predict the effect size. The
	format of this matrix follows
	\code{\link[metaSEM]{as.mxMatrix}}. The regression coefficients can be
	constrained equally by using the same labels.
}
  \item{RE2.constraints}{A scalar or a \eqn{1}{1} x \eqn{1}{1} matrix
	specifying the variance componets of the random effects. The default
	is that the variance components is free. The format of this matrix
	follows \code{\link[metaSEM]{as.mxMatrix}}. Elements of the variance
	components can be constrained equally by using the same label. 
}
  \item{RE2.lbound}{A scalar or a \eqn{1}{1} x \eqn{1}{1} matrix of lower bound on the level-2 variance component of the random effects. 
}
  \item{RE3.constraints}{A scalar of a \eqn{1}{1} x \eqn{1}{1} matrix
	specifying the variance componets of the random effects at
	level-3. The default is that the variance components is free. The format of this matrix
	follows \code{\link[metaSEM]{as.mxMatrix}}. Elements of the variance
	components can be constrained equally by using the same label. 
}
  \item{RE3.lbound}{A scalar or a \eqn{1}{1} x \eqn{1}{1} matrix of lower bound on the level-3 variance component of the random effects. 
}
  \item{intervals.type}{Either \code{z} (default if missing) or
	\code{LB}. If it is \code{z}, it calculates the 95\% Wald confidence
	intervals (CIs) based on the z statistic. If it is \code{LB}, it
	calculates the 95\% likelihood-based CIs on the
	parameter estimates. Note that the z values and their
	associated p values are based on the z statistic. They are not
	related to the likelihood-based CIs.
}
  \item{I2}{Possible options are \code{"I2q"}, \code{"I2hm"},
	\code{"I2am"} and \code{"ICC"}. They represent the \code{I2} calculated by using a
	typical within-study sampling variance from the Q statistic, the
	harmonic mean, the arithmatic mean of the within-study sampling
	variances, and the intra-class correlation. More than one options are possible. If \code{intervals.type="LB"}, 95\% confidence intervals on the heterogeneity indices will be constructed. 
  }
  \item{R2}{Logical. If \code{TRUE} and there are predictors, R2 is
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
  \deqn{y_{ij} = \beta_0 + \mathbf{\beta'}*\mathbf{x}_{ij} + u_{(2)ij} + u_{(3)j} + e_{ij}
  } where \eqn{y_{ij}} is the effect size for the ith study in the jth cluster,
  \eqn{\beta_0} is the intercept, \eqn{\mathbf{\beta}} is the
  regression coefficients, \eqn{\mathbf{x}_{ij}} is a vector of predictors, \eqn{u_{(2)ij} \sim N(0, \tau^2_2)}{u_{(2)ij}~ N(0,
	  tau^2_2)} and \eqn{u_{(3)j} \sim N(0, \tau^2_3)}{u_{(3)j}~ N(0,
	  tau^2_3)} are the level-2 and level-3 heterogeneity variances,
	respectively, and \eqn{e_{ij} \sim N(0, v_{ij})}{e_{ij}~ N(0,
	  v_{ij})} is the conditional known sampling variance.

 \code{meta3()} does not differentiate between level-2 or level-3
 variables in \code{x} since both variables are treated as a design
 matrix. When there are missing values in \code{x}, the data will be
 deleted. \code{meta3X()} treats the predictors \code{x2} and \code{x3}
 as level-2 and level-3 variables. Thus, their means and covariance
 matrix will be estimated. Missing values in \code{x2} and \code{x3}
 will be handled by (full information) maximum likelihood (FIML) in \code{meta3X()}. Moreover,
 auxiliary variables \code{av2} at level-2 and \code{av3} at level-3 may
 be included to improve the estimation. Although \code{meta3X()} is more
 flexible in handling missing covariates, it is more likely to encounter
 estimation problems. 	
}

\references{
  Cheung, M. W. L. (in press). Modeling dependent effect sizes with three-level meta-analyses: A structural equation modeling approach. \emph{Psychological Methods}.

  Enders, C. K. (2010). \emph{Applied missing data analysis}. New York:
  Guilford Press.

  Graham, J. (2003). Adding missing-data-relevant variables to
 FIML-based structural equation models. \emph{Structural Equation
   Modeling: A Multidisciplinary Journal}, \bold{10(1)}, 80-100.

 Konstantopoulos, S. (2011). Fixed effects and variance components
  estimation in three-level meta-analysis. \emph{Research Synthesis
  Methods}, \bold{2}, 61-76.
}
\author{Mike W.-L. Cheung <mikewlcheung@nus.edu.sg>
%%  ~~who you are~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link[metaSEM]{reml3}}, \code{\link[metaSEM]{Cooper03}}, \code{\link[metaSEM]{Bornmann07}}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ meta-analysis }
