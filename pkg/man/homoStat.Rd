\name{homoStat}
\alias{homoStat}
\title{Test the Homogeneity of Effect Sizes
}
\description{It tests the homogeneity of univariate and multivariate effect sizes.
}
\usage{
homoStat(y, v)
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
	major as used by \code{\link[OpenMx]{vech}}. It is assumed that
	there is no missing value in \code{v} if \code{y} is complete. If there are missing values in \code{v}
	due to the missingness on \code{y}, the missing values in
	\code{v} will be removed automatically.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{A list of
  \item{Q}{Q statistic on the null hypothesis of homogeneity of effect
	sizes. It has an approximate chi-square distribution under the null
	hypothesis.}
  \item{Q.df}{Degrees of freedom of the Q statistic}
  \item{pval}{p value on the test of homogeneity of effect sizes}
}
\references{
Becker, B. J. (1992). Using results from replicated studies to
  estimate linear models. \emph{Journal of Educational Statistics},
  \bold{17}, 341-362. 
  
Cheung, M. W. L. (2010). Fixed-effects meta-analyses as multiple-group
structural equation models. \emph{Structural Equation Modeling},
\bold{17}, 481-509.

Cochran, W. G. (1954). The combination of estimates from different experiments. \emph{Biometrics}, \bold{10}, 101-129.
}
\author{Mike W.-L. Cheung <mikewlcheung@nus.edu.sg>
}
\note{It is usually called internally by \code{\link[metaSEM]{meta}}.
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{  \code{\link[metaSEM]{meta}}
}
\examples{
with( Hox02, homoStat(yi, vi) )

with( HedgesOlkin85, homoStat(y=cbind(d_att, d_ach),
      v=cbind(var_att, cov_att_ach, var_ach)) )
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ meta-analysis }
