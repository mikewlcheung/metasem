\name{uniR1}
\alias{uniR1}
\title{First Stage analysis of the univariate R (uniR) approach
}
\description{It conducts the first stage analysis of the uniR analysis
  by pooling elements of the correlation coefficients individually. 
}
\usage{
uniR1(Cor, n, ...)
}
\arguments{
  \item{Cor}{A list of correlation matrices
}
  \item{n}{A vector of sample sizes
}
  \item{\dots}{Further arguments which are currently ignored}
}
\value{An object of class \code{uniR1} of the original data, the sample
  sizes, the harmonic mean of sample sizes, the average correlation
  matrix, the standard errors of the correlation matrix, and the
  standard deviations (heterogeneity) of the correlation matrix. 
}
\details{
  This function implements the univariate r approach proposed by Viswesvaran
  and Ones (1995) to conduct meta-analytic structural equation modeling
  (MASEM). It uses Schmidt and Hunter's approach to combine correlation
  coefficients. It is included in this package for research interests. The two-stage
  structural equation modeling (TSSEM) approach is preferred (e.g.,
  Cheung, 2015; Cheung & Chan, 2005).
}
\references{
  Cheung, M. W.-L. (2015). \emph{Meta-analysis: A structural equation
	modeling approach}. Chichester, West Sussex: John Wiley & Sons, Inc.

  Cheung, M. W.-L., & Chan, W. (2005). Meta-analytic structural equation
  modeling: A two-stage approach. \emph{Psychological Methods},
  \bold{10}, 40-64.

  Schmidt, F. L., & Hunter, J. E. (2015). \emph{Methods of
  meta-analysis: Correcting error and bias in research findings (3rd
  ed.)}. Thousand Oaks, CA: Sage.

  Viswesvaran, C., & Ones, D. S. (1995). Theory testing: Combining
  psychometric meta-analysis and structural equations modeling.
  \emph{Personnel Psychology}, \bold{48}, 865-885. 
}
\author{Mike W.-L. Cheung <mikewlcheung@nus.edu.sg>
}
\seealso{ \code{\link[metaSEM]{uniR2}}, \code{\link[metaSEM]{Becker09}}
}
\keyword{uniR}
