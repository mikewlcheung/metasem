\name{smdMTS}
\alias{smdMTS}
\title{Compute Effect Sizes for Multiple Treatment Studies
}
\description{It computes the standardized mean differences and their asymptotic sampling
  covariance matrix for \emph{k} multiple treatment studies. The first group is
  assumed as the control group. 
}
\usage{
smdMTS(m, v, n, homogeneity=c("variance", "none"), bias.adjust=TRUE, 
       all.comparisons=FALSE, list.output=TRUE, lavaan.output=FALSE) 
}
\arguments{
  \item{m}{A vector of \emph{k} sample means.}
  \item{v}{A vector of \emph{k} sample variances.}
  \item{n}{A vector of \emph{k} sample sizes.}
  \item{homogeneity}{If it is \code{variance} (the default), homogeneity of variances
	is assumed. The common standard deviation is used as the
	standardizer in calculating the effect sizes. If it is \code{none}, homogeneity of variances is not
	assumed. The standard deviation of the first group is used as the
	standardizer in calculating the effect sizes.}
  \item{bias.adjust}{If it is \code{TRUE} (the default), the effect sizes are
	adjusted for small bias by multiplying \eqn{1-3/(4*(n1+n2)-9)}.}  
  \item{all.comparisons}{If it is \code{FALSE} (the default), all groups
	(except the first group) are compared against the first group. If it
	is \code{TRUE}, all pairwise comparisons are calculated. This may be
	useful in network meta-analysis.}
  \item{list.output}{If it is \code{TRUE} (the default), the effect
	sizes and their sampling covariance matrix are outputed as a
	list. If it is \code{FALSE}, they will be stacked into a vector.}  
  \item{lavaan.output}{If it is \code{FALSE} (the default), the effect
	sizes and its sampling covariance matrix are reported. If it is
	\code{TRUE}, it outputs the fitted
	\code{\link[lavaan]{lavaan-class}} object.}
}
\details{
  Gleser and Olkin (2009) introduce formulas to calculate the
  standardized mean differences and their sampling covariance matrix for
  multiple treatment studies under the assumption of homogeneity of the
  covariance matrix. This function uses a structural equation modeling (SEM)
  approach introduced in Chapter 3 of Cheung (2015) to calculate the
  same estimates. The SEM approach is more flexible in three ways: (1)
  it allows homogeneity of variances or not; (2) it allows users to test
  the assumption of homogeneity of variances by checking the fitted
	\code{\link[lavaan]{lavaan-class}} object; and (3) it may calculate all
  pairwise comparisons.
}
\references{
  Cheung, M. W.-L. (2015). \emph{Meta-analysis: A structural equation modeling approach}. Chichester, West Sussex: John Wiley & Sons, Inc.

  Cheung, M. W.-L. (2017). Computing multivariate effect sizes and their
  sampling covariance matrices with structural equation modeling:
  Theory, examples, and computer simulations. \emph{Manuscript under
	review}.
  
  Gleser, L. J., & Olkin, I. (2009). Stochastically dependent effect sizes. In H. Cooper, L. V. Hedges, & J. C. Valentine (Eds.), \emph{The handbook of research synthesis and meta-analysis}. (2nd ed., pp. 357-376). New York: Russell Sage Foundation.
}
\author{Mike W.-L. Cheung <mikewlcheung@nus.edu.sg>
}
\seealso{\code{\link[metaSEM]{Gleser94}}
}

\examples{
\dontrun{  
## Sample means for groups 1 to 3
m <- c(5,7,9)

## Sample variances
v <- c(10,11,12)

## Sample sizes
n <- c(50,52,53)

## Assuming homogeneity of variances
smdMTS(m, v, n, homogeneity = "var", bias.adjust=TRUE, all.comparisons=FALSE,
       lavaan.output=FALSE)

## Not assuming homogeneity of variances and comparing all pairwise groups
## Please note that the SD of the first group is used as the standardizer    
smdMTS(m, v, n, homogeneity = "none", bias.adjust=TRUE, all.comparisons=TRUE,
       lavaan.output=FALSE)

## Output the fitted lavaan model
## It provides a likelihood ratio test to test the null hypothesis of
## homogeneity of variances.    
fit <- smdMTS(m, v, n, homogeneity = "var", bias.adjust=FALSE, all.comparisons=FALSE,
              lavaan.output=TRUE)

lavaan::summary(fit)
    
lavaan::parameterestimates(fit)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ meta-analysis }
