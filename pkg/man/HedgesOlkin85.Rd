\name{HedgesOlkin85}
\alias{HedgesOlkin85}
\docType{data}
\title{Effects of open education reported by Hedges and Olkin (1985)
}
\description{
  Effects of open education on attitude toward school and on reading achievement reported by Hedges and Olkin (1985)
}
\usage{data(HedgesOlkin85)}
\format{
  The format is:
 num [1:4, 1:6] 1 2 3 4 0.458 0.363 0.162 0.294 0.1 0.241 ...
 - attr(*, "dimnames")=List of 2
  ..$ : NULL
  ..$ : chr [1:6] "study" "d_att" "d_ach" "var_att" ...
}
\details{
        The variables are:
  \describe{
	\item{study}{Study number}
	\item{d_att}{Standardized mean difference on \emph{attitude}}
	\item{d_ach}{Standardized mean difference on \emph{achievement}}
	\item{var_att}{Sampling variance of the effect size of \emph{attitude}}
	\item{cov_att_ach}{Sampling covariance between the effect sizes}
	\item{var_ach}{Sampling variance of the effect size of \emph{achievement}}
}
}
\source{
Hedges, L. V., & Olkin, I. (1985). \emph{Statistical methods for meta-analysis.} Orlando, FL: Academic Press. 
}
\references{
  Cheung, M. W. L. (2010). Fixed-effects meta-analyses as multiple-group structural equation models. \emph{Structural Equation Modeling}, \bold{17}, 481-509.
}
\examples{
data(HedgesOlkin85)
## maybe str(HedgesOlkin85) ; plot(HedgesOlkin85) ...
}
\keyword{datasets}
