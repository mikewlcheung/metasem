\name{reml3}
\alias{reml3}
\title{Extract Variance Components in Three-Level Univariate
  Meta-Analysis with Restricted (Residual) Maximum
  Likelihood Estimation
}
\description{
It estimates the variance components of random-effects in three-level
univariate meta-analysis with restricted (residual) maximum likelihood (REML) estimation method.
}
\usage{
reml3(y, v, cluster, x, data, RE2.startvalue=0.1, RE2.lbound=1e-10,
      RE3.startvalue=RE2.startvalue, RE3.lbound=RE2.lbound, RE.equal=FALSE,
      intervals.type=c("z", "LB"), model.name="Variance component with REML",
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
  \item{RE2.startvalue}{Starting value for the level-2 variance.
}
  \item{RE2.lbound}{Lower bound for the level-2 variance.
}
  \item{RE3.startvalue}{Starting value for the level-3 variance.
}
  \item{RE3.lbound}{Lower bound for the level-3 variance.
}
  \item{RE.equal}{Logical. Whether the variance components at level-2
	and level-3 are constrained equally.
}
  \item{intervals.type}{Either \code{z} (default if missing) or
	\code{LB}. If it is \code{z}, it calculates the 95\% Wald confidence
	intervals (CIs) based on the z statistic. If it is \code{LB}, it
	calculates the 95\% likelihood-based CIs on the
	parameter estimates. Note that the z values and their
	associated p values are based on the z statistic. They are not
	related to the likelihood-based CIs.
}
  \item{model.name}{A string for the model name in \code{\link[OpenMx]{mxModel}}. 
}
  \item{suppressWarnings}{Logical. If \code{TRUE}, warnings are
	suppressed. Argument to be passed to \code{\link[OpenMx]{mxRun}}.
}
  \item{\dots}{Futher arguments to be passed to \code{\link[OpenMx]{mxRun}}}
}
\details{
Restricted (residual) maximum likelihood obtains the parameter estimates on the transformed data that do not include the fixed-effects parameters. A transformation matrix \eqn{M=I-X(X'X)^{-1}X}{M=I-X(X'X)^{-1}X'} is created based on the design matrix \eqn{X}{X} which is just a column vector when there is no predictor in \code{x}. The last \eqn{N}{N} redundant rows of \eqn{M}{M} is removed where \eqn{N}{N} is the rank of \eqn{X}{X}. After pre-multiplying by \eqn{M} on \code{y}, the parameters of fixed-effects are removed from the model. Thus, only the parameters of random-effects are estimated.

An alternative but equivalent approach is to minimize the
  -2*log-likelihood function: \deqn{
	\log(\det|V+T^2|)+\log(\det|X'(V+T^2)^{-1}X|)+(y-X\hat{\alpha})'(V+T^2)^{-1}(y-X\hat{\alpha})}{
	log(det|V+T^2|)+log(det|X'(V+T^2)^{-1}X|)+(y-X\alpha)'(V+T^2)^{-1}(y-X*\alpha)}
  where \eqn{V}{V} is the known conditional sampling covariance matrix
  of \eqn{y}{y}, \eqn{T^2}{T^2} is the variance component combining
  level-2 and level-3 random effects, and \eqn{\hat{\alpha}=(X'(V+T^2)^{-1}X)^{-1}
  X'(V+T^2)^{-1}y}{\hat{\alpha}=(t(X)(V+T^2)^{-1}X)^{-1}t(X)(V+T^2)^{-1}y}. \code{reml()}
  minimizes the above likelihood function to obtain the parameter estimates.

}
\value{An object of class \code{reml} with a list of
  \item{call}{Object returned by \code{\link[base]{match.call}}}
  \item{data}{A data matrix of y, v and x  }
  % \item{no.y}{No. of effect sizes}
  % \item{no.x}{No. of predictors}
  % \item{miss.vec}{A vector indicating missing data. Studies will be removed before the analysis if they are \code{TRUE}}
  \item{mx.fit}{A fitted object returned from	\code{\link[OpenMx]{mxRun}}}
}
\references{
  Cheung, M. W. L. (in press). Implementing restricted maximum likelihood (REML) estimation in structural equation models. \emph{Structural Equation Modeling.}

  Cheung, M. W. L. (2012). Modeling three-level meta-analyses as
  structural equation models. \emph{Manuscript in preparation}.
	
  Mehta, P. D., & Neale, M. C. (2005). People Are Variables Too: Multilevel Structural Equations Modeling. \emph{Psychological
  	Methods}, \bold{10(3)}, 259-284.

  Searle, S. R., Casella, G., & McCulloch, C. E. (1992). \emph{Variance components}. New York: Wiley.
  
  Viechtbauer, W. (2005). Bias and efficiency of meta-analytic variance estimators in the random-effects model. \emph{Journal of Educational and Behavioral Statistics}, \bold{30(3)}, 261-293.
  
}
\author{Mike W.-L. Cheung <mikewlcheung@nus.edu.sg>
%%  ~~who you are~~
}
\note{ \code{reml} is more computional intensive than \code{meta}. Moreover, \code{reml} is more
  likely than \code{meta} to encounter errors during optimization. Since
  a likelihood function is directly employed to obtain the parameter
  estimates, there is no number of studies and number of observed statistics
  returned by \code{\link[OpenMx]{mxRun}}. Ad-hoc calculations are used
  to modify \code{mx.fit@runstate$objectives[[1]]@numObs} and \code{mx.fit@runstate$objectives[[1]]@numStats}.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link[metaSEM]{meta3}}, \code{\link[metaSEM]{reml}}, \code{\link[metaSEM]{Cooper03}}, \code{\link[metaSEM]{Bornmann07}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ meta-analysis }
