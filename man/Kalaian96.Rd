\name{Kalaian96}
\alias{Kalaian96}
\docType{data}
\title{Multivariate effect sizes reported by Kalaian and Raudenbush (1996)
}
\description{
  This data set includes 47 multivariate effect sizes reported by
  Kalaian and Raudenbush (1996, Table 1). 
}
\usage{data(Kalaian96)}

\details{
A list of data with the following structure:
\describe{
  \item{Study}{Study name}
  \item{Year}{Year of publication}
  \item{n_e}{Sample size of the experimental group}
  \item{n_c}{Sample size of the control group}
  \item{dSAT_V}{Standardized mean difference of the Scholastic Aptitude
Test (SAT) on verbal}
  \item{dSAT_M}{Standardized mean difference of SAT on math}
  \item{var_V}{Sampling variance of \code{dSAT_V}}
  \item{cov_VM}{Sampling covariance of \code{dSAT_V} and \code{dSAT_M}
with a common correlation of 0.66}
  \item{var_M}{Sampling variance of \code{dSAT_M}}
  \item{Hr}{Hours of training}
  \item{ETS}{Educational Testing Service}
  \item{Study_type}{Either \code{Randomized}, \code{Matched} or
\code{Nonequivalent comparison}}
  \item{Home_work}{Home work}
}
}

\source{
  Kalaian, H. A., & Raudenbush, S. W. (1996). A multivariate mixed linear model for meta-analysis. \emph{Psychological Methods}, \emph{1}(3), 227-235. https://doi.org/10.1037/1082-989X.1.3.227
}

\examples{
\dontrun{
data(Kalaian96)
}
}
\keyword{datasets}
