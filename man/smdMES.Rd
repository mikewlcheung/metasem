\name{smdMES}
\alias{smdMES}
\title{Compute Effect Sizes for Multiple End-point Studies
}
\description{It computes the standardized mean differences and their asymptotic sampling
  covariance matrix for two multiple end-point studies with \emph{p}
  effect sizes. 
}
\usage{
smdMES(m1, m2, V1, V2, n1, n2,
       homogeneity=c("covariance", "correlation", "none"), 
       bias.adjust=TRUE, list.output=TRUE, lavaan.output=FALSE) 
}
\arguments{
  \item{m1}{A vector of \emph{p} sample means of the first group.}
  \item{m2}{A vector of \emph{p} sample means of the second group.}
  \item{V1}{A \emph{p} by \emph{p} sample covariance matrix of the first
	group.}
  \item{V2}{A \emph{p} by \emph{p} sample covariance matrix of the second
	group.}  
  \item{n1}{The sample size of the first group.}
  \item{n2}{The sample size of the second group.}
  \item{homogeneity}{If it is \code{covariance} (the default),
	homogeneity of covariance matrices is assumed. The common standard
	deviations are used as the standardizers in calculating the effect
	sizes. If it is \code{correlation}, homogeneity of correlation is
	not assumed. The standard deviations of the first group are used as the
	standardizer in calculating the effect sizes. If it is \code{none},
	no homogeneity assumption is made. The standard deviations of the first group are used as the
	standardizer in calculating the effect sizes.}
  \item{bias.adjust}{If it is \code{TRUE} (the default), the effect sizes are
	adjusted for small bias by multiplying \eqn{1-3/(4*(n1+n2)-9)}.}
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
  multiple end-point studies under the assumption of homogeneity of the
  covariance matrix. This function uses a structural equation modeling (SEM)
  approach introduced in Chapter 3 of Cheung (2015) to calculate the
  same estimates. The SEM approach is more flexible in two ways: (1)
  it allows homogeneity of covariance or correlation matrices or not; and
  (2) it allows users to test this assumption by checking the fitted
	\code{\link[lavaan]{lavaan-class}} object.
}
\references{
  Cheung, M. W.-L. (2015). \emph{Meta-analysis: A structural equation
	modeling approach}. Chichester, West Sussex: John Wiley & Sons, Inc.

  Cheung, M. W.-L. (in press). Computing multivariate effect sizes and
  their sampling covariance matrices with structural equation modeling:
  Theory, examples, and computer simulations. \emph{Frontiers in
	Psychology}.
  
  Gleser, L. J., & Olkin, I. (2009). Stochastically dependent effect
  sizes. In H. Cooper, L. V. Hedges, & J. C. Valentine (Eds.), \emph{The handbook of research synthesis and meta-analysis}. (2nd ed., pp. 357-376). New York: Russell Sage Foundation.
}
\author{Mike W.-L. Cheung <mikewlcheung@nus.edu.sg>
}
\seealso{\code{\link[metaSEM]{Gleser94}}
}

\examples{
\dontrun{    
## Sample means for the two constructs in Group 1  
m1 <- c(2.5, 4.5)

## Sample means for the two constructs in Group 2     
m2 <- c(3, 5)

## Sample covariance matrix in Group 1    
V1 <- matrix(c(3,2,2,3), ncol=2)

## Sample covariance matrix in Group 2
V2 <- matrix(c(3.5,2.1,2.1,3.5), ncol=2)

## Sample size in Group 1
n1 <- 20

## Sample size in Group 2    
n2 <- 25

## SMD with the assumption of homogeneity of covariance matrix    
smdMES(m1, m2, V1, V2, n1, n2, homogeneity="cov", bias.adjust=TRUE,
       lavaan.output=FALSE)

## SMD with the assumption of homogeneity of correlation matrix 
smdMES(m1, m2, V1, V2, n1, n2, homogeneity="cor", bias.adjust=TRUE,
       lavaan.output=FALSE)

## SMD without any assumption of homogeneity
smdMES(m1, m2, V1, V2, n1, n2, homogeneity="none", bias.adjust=TRUE,
       lavaan.output=FALSE)

## Output the fitted lavaan model
## It provides a likelihood ratio test to test the null hypothesis of
## homogeneity of variances.     
fit <- smdMES(m1, m2, V1, V2, n1, n2, homogeneity="cov", bias.adjust=TRUE,
              lavaan.output=TRUE)

lavaan::summary(fit)

lavaan::parameterestimates(fit)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ meta-analysis }
